#!/usr/bin/env python

from __future__ import division

import numpy as np
from builtins import int
from dipy.tracking.streamline import set_number_of_points, values_from_volume
from nibabel.affines import apply_affine
from sklearn.cluster import DBSCAN


def streamlines_merge(*args, **kwargs):
    list_pack = list(args)
    streamlines_total = [np.asfarray(line) for sublist in list_pack for line in sublist]

    if kwargs['hdr']:
        kwargs['hdr']['nb_streamlines'] = len(streamlines_total)
        return streamlines_total, kwargs['hdr']
    else:
        return streamlines_total


def streamlines_resample(streamlines, perc):
    resampled = [set_number_of_points(s, max(int(len(s) * perc / 100.), 2)) for s in streamlines]

    return resampled


def streamlines_sort(streamlines):
    streamlines_sorted = []
    template = streamlines[0]
    for s in streamlines:
        dist_norm = np.linalg.norm(s[0] - template[0]) + np.linalg.norm(s[-1] - template[-1])
        dist_flipped = np.linalg.norm(s[-1] - template[0]) + np.linalg.norm(s[0] - template[-1])
        if dist_flipped < dist_norm:
            streamlines_sorted.append(s[::-1])
        else:
            streamlines_sorted.append(s)

    return streamlines_sorted


def streamlines_mapvolume(streamlines, volume, affine):
    inverse = np.linalg.inv(affine)
    streamlines = [apply_affine(inverse, s) for s in streamlines]
    mapping = values_from_volume(volume, streamlines)

    return mapping


def streamlines_clusters(sim):
    db_labels = DBSCAN(min_samples=2, eps=0.2, metric="precomputed", n_jobs=-1).fit_predict(sim / sim.max())

    return db_labels
