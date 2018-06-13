#!/usr/bin/env python

from __future__ import division

import sys

import nibabel as nib
import numpy as np
from builtins import int
from dipy.viz import fvtk
from six import iteritems

try:
    import itertools.izip as zip
except ImportError:
    pass
from itertools import islice
from .varifolds import get_prototypes
from dipy.viz.colormap import create_colormap
from scipy.ndimage import generate_binary_structure
from scipy.ndimage.morphology import binary_dilation
from .streamlines_operations import streamlines_resample, streamlines_mapvolume, streamlines_sort

try:
    from psutil import virtual_memory
    from joblib import Parallel, delayed, cpu_count
    from fast_shift import fast_shift3d_parallel

    parallelize = True
except ImportError:
    from fast_shift import fast_shift3d

    parallelize = False


def combinator(*args):
    args = list(args)
    n = args.pop()
    cur = np.arange(n)
    cur = cur[:, None]

    while args:
        d = args.pop()
        cur = np.kron(np.ones((d, 1)), cur)
        front = np.arange(d).repeat(n)[:, None]
        cur = np.column_stack((front, cur))
        n *= d

    return cur


def batch_iterable(iterable, size):
    iterable = iter(iterable)
    batch = list(islice(iterable, size))
    while batch:
        yield batch
        batch = list(islice(iterable, size))


def compute_cone(mask, angle_dir):
    dims = mask.shape

    struct = generate_binary_structure(3, 3)
    borders = np.logical_xor(mask, binary_dilation(mask, struct))

    dir_x = np.cos(angle_dir[1]) * np.cos(angle_dir[0])
    dir_y = np.cos(angle_dir[1]) * np.sin(angle_dir[0])
    dir_z = np.sin(angle_dir[1])

    center = (int(dims[0] / 2), int(dims[1] / 2), int(dims[2] / 2))
    angle = combinator(dims[0], dims[1], dims[2]) - center
    angle = np.arccos(np.sum(angle * (dir_x, dir_y, dir_z), axis=1) / np.sqrt(np.sum((angle ** 2), axis=1)))
    angle[np.ravel_multi_index(center, dims)] = 0

    angle = np.reshape(angle, dims)

    k = np.pi
    angle[angle < (np.amax(angle) * (1 - k / (2 * np.pi)))] = 0
    angle = angle - np.min(angle[np.nonzero(angle)])
    angle = angle * 255 / np.amax(angle)

    dist = np.transpose(np.nonzero(borders)) - center
    if parallelize:
        num_cores = cpu_count()
        mem = virtual_memory()
        dt = angle.dtype
        if sys.maxsize == 2 ** 63 - 1:
            x = 64
        else:
            x = 32
        max_array_dim = angle.size * x / dt.itemsize
        num_array = int((mem.total / max_array_dim) * 0.9)
        cones = Parallel(n_jobs=num_cores)(
            delayed(fast_shift3d_parallel)(angle, i) for i in batch_iterable(dist, num_array))
        cone = np.maximum.reduce(cones)
    else:
        cone = np.zeros(dims, angle.dtype)
        for i in dist:
            cone = np.maximum(cone, fast_shift3d(angle, i))
    return cone


class SpatialFuzzySets:
    # LPI
    cones_direction = {'anterior_of': [np.pi / 2, 0], 'posterior_of': [-np.pi / 2, 0], 'left_of': [np.pi, 0],
                       'right_of': [0, 0], 'superior_of': [0, np.pi / 2], 'inferior_of': [0, -np.pi / 2]}

    def __init__(self, key, value, parcellation, affine, labels_dict, show):
        self.ID = key
        self.dict_conf = value
        self.affine = affine
        self.show = show
        self.parcellation = parcellation
        self.labels_dict = labels_dict

        orientation = nib.aff2axcodes(self.affine)
        if orientation[0] == 'R':
            SpatialFuzzySets.cones_direction['right_of'] = [np.pi, 0]
            SpatialFuzzySets.cones_direction['left_of'] = [0, 0]
        if orientation[1] == 'A':
            SpatialFuzzySets.cones_direction['anterior_of'] = [-np.pi / 2, 0]
            SpatialFuzzySets.cones_direction['posterior_of'] = [np.pi / 2, 0]
        if orientation[2] == 'S':
            SpatialFuzzySets.cones_direction['superior_of'] = [0, -np.pi / 2]
            SpatialFuzzySets.cones_direction['inferior_of'] = [0, np.pi / 2]

        self.cones = []
        for key, value in iteritems(self.dict_conf[0]):
            if value is not None:
                direction = SpatialFuzzySets.cones_direction[key]
                for v in value:
                    labels = []
                    for freesurfer_key in labels_dict:
                        if v.lower() in freesurfer_key.lower():
                            labels.append(labels_dict[freesurfer_key])

                    if labels:
                        mask = np.isin(parcellation, labels)
                        if np.count_nonzero(mask) != 0:
                            cone = compute_cone(mask, direction)
                            self.cones.append(cone)

    def __repr__(self):
        return "{}({},{},{},{},{},{})".format(self.__class__.__name__, self.ID, self.dict_conf, self.parcellation,
                                              self.affine, self.labels_dict, self.show)

    def __str__(self):
        return "{}({},{},{},{},{},{})".format(self.__class__.__name__, 'Bundle_ID', 'Bundle_config', 'Parcellation',
                                              'Affine', 'Labels', 'Show')

    def get_intersection(self):
        intersection = np.full(self.cones[0].shape, True)
        for cone in self.cones:
            intersection = np.logical_and(cone != 0, intersection)

        return intersection, self.affine

    def put_streamlines(self, streamlines):
        self.streamlines = streamlines

    @property
    def streamlines(self):
        return self._streamlines

    @streamlines.setter
    def streamlines(self, value):
        if not value:
            raise ValueError('No streamlines found after the initial segmentation. Incorrect configuration.')
        self._streamlines = streamlines_sort(value)

    def downsampling(self, perc):
        self.streamlines = streamlines_resample(self.streamlines, perc)

    def clustering(self, fa, fa_affine):
        kernel_dict = self.dict_conf[1]

        lambda_a = kernel_dict['lambdaA']
        lambda_b = kernel_dict['lambdaB']
        lambda_w = kernel_dict['lambdaW']
        lambda_m = kernel_dict['lambdaM']

        self.prototypes, self.prototypes_keys, self.labels = get_prototypes(self.streamlines, fa, fa_affine, lambda_a,
                                                                            lambda_b, lambda_w, lambda_m)

        if self.show:
            colormap = create_colormap(np.linspace(0, 1, len(self.prototypes)), name='jet')
            ren = fvtk.ren()
            fvtk.clear(ren)
            ren.SetBackground(0, 0, 0)
            fvtk.camera(ren, (0, 70, 45))
            fvtk.add(ren, fvtk.streamtube(self.streamlines, fvtk.colors.white, linewidth=0.1, opacity=0.4))
            fvtk.add(ren, fvtk.streamtube(self.prototypes, colormap, linewidth=0.4))
            fvtk.show(ren)

    def segmentation(self, threshold):

        first = True
        for cone in self.cones:
            mapping = streamlines_mapvolume(self.prototypes, cone, self.affine)
            mapping = [np.array(m) for m in mapping]
            mapping = np.asarray([np.mean(m) if 0 not in m else 0 for m in mapping])
            if first:
                fuzzy_centroids = mapping
                first = False
            else:
                fuzzy_centroids = np.minimum(fuzzy_centroids, mapping)

        if np.amax(fuzzy_centroids) != np.amin(fuzzy_centroids):
            fuzzy_centroids = fuzzy_centroids - np.amin(fuzzy_centroids)

        fuzzy_centroids = (fuzzy_centroids / np.amax(fuzzy_centroids))

        fuzzy_mapping = [0 if label == -1 else fuzzy_centroids[self.prototypes_keys.index(str(label))] for label in
                         self.labels]

        self.final = [e1 for (e1, e2) in zip(self.streamlines, fuzzy_mapping) if e2 > threshold]
        self.fuzzy_per_line = {
            'fuzzyness': [np.array(e2, dtype="f4") for (e1, e2) in zip(self.streamlines, fuzzy_mapping) if
                          e2 > threshold]}

        if self.show:
            colormap = create_colormap(np.asarray(fuzzy_mapping), name='plasma')
            ren = fvtk.ren()
            fvtk.clear(ren)
            ren.SetBackground(255, 255, 255)
            fvtk.camera(ren, (0, 70, -10))
            for i, line in enumerate(self.streamlines):
                fvtk.add(ren, fvtk.streamtube([line], colormap[i], linewidth=0.1))
            fvtk.show(ren)

            fvtk.clear(ren)
            for line in self.streamlines:
                flag = [np.array_equal(fin, line) for fin in self.final]
                if any(flag):
                    fvtk.add(ren, fvtk.streamtube([line], fvtk.colors.blue, linewidth=0.2))
                else:
                    fvtk.add(ren, fvtk.streamtube([line], fvtk.colors.red, linewidth=0.1, opacity=0.1))
            fvtk.show(ren)

    @property
    def final(self):
        if self._flag_f:
            raise ValueError('Threshold value too high: no fibers left.')
        return self._final

    @final.setter
    def final(self, value):
        if not value:
            self._flag_f = True
            self._final = self.streamlines
        else:
            self._flag_f = False
            self._final = value

    @property
    def fuzzy_per_line(self):
        return self._fuzzy_per_line

    @fuzzy_per_line.setter
    def fuzzy_per_line(self, value):
        self._fuzzy_per_line = value
