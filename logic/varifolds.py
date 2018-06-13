#!/usr/bin/env python

from __future__ import division

import sys

import numpy as np
from builtins import range
from numpy.linalg import norm
from psutil import virtual_memory
from six import itervalues

try:
    import itertools.izip as zip
except ImportError:
    pass
from itertools import count
from sklearn.metrics.pairwise import pairwise_distances as pd_sk
from .streamlines_operations import (streamlines_mapvolume, streamlines_clusters)


def pairwise_distances(s1, s2, operation):
    if operation == 'Fibers':
        return (s1 * s1).sum(1)[:, None] - 2.0 * np.dot(s1, s2.T) + (s2 * s2).sum(1)
    if operation == 'FA':
        return pd_sk(s1.reshape(-1, 1), s2.reshape(-1, 1), n_jobs=1) ** 2


def inner_varifolds(centa, tanga, faa, centb, tangb, fab, lambda_w, lambda_m):
    exp_dist = np.exp(-pairwise_distances(centa, centb, operation='Fibers') / lambda_w)

    exp_fa = np.exp(-pairwise_distances(faa, fab, operation='FA') / lambda_m)

    norm_tanga = np.sqrt(norm(tanga, axis=1)) + 1.0e-9
    tanga_normalised = tanga / norm_tanga[:, None]

    norm_tangb = np.sqrt(norm(tangb, axis=1)) + 1.0e-9
    tangb_normalised = tangb / norm_tangb[:, None]

    gram_matrix = np.power(tanga_normalised.dot(tangb_normalised.T), 2.0)

    return (exp_fa * exp_dist * gram_matrix).sum()


def inner_product(starta, enda, centa, tanga, faa, startb, endb, centb, tangb, fab, lambda_a, lambda_b, lambda_w,
                  lambda_m):
    exp_start = np.exp(- ((starta - startb) ** 2).sum() / lambda_a)

    exp_end = np.exp(- ((enda - endb) ** 2).sum() / lambda_b)

    var = inner_varifolds(centa, tanga, faa, centb, tangb, fab, lambda_w, lambda_m)

    return exp_start * exp_end * var


def varifolds_calc(streamlines, fa, affine, lambda_a=np.inf, lambda_b=np.inf, lambda_w=np.inf, lambda_m=np.inf):
    n_stream = len(streamlines)

    start, end, cent, tang = ([] for _ in range(4))

    for s in streamlines:
        start.append(np.array(s[0]))
        end.append(np.array(s[-1]))
        cent.append((s[:-1, :] + s[1:, :]) / 2.0)
        tang.append(np.diff(s, axis=0))

    mapping = streamlines_mapvolume(cent, fa, affine)
    mapping = [np.array(m) for m in mapping]

    inners = np.zeros((n_stream, n_stream))
    for i in range(n_stream):
        j = i
        while j < n_stream:
            inners[i][j] = inner_product(start[i], end[i], cent[i], tang[i], mapping[i], start[j], end[j], cent[j],
                                         tang[j], mapping[j], lambda_a, lambda_b, lambda_w, lambda_m)
            inners[j][i] = inners[i][j]
            j = j + 1

    distances = np.zeros((n_stream, n_stream))

    for i in range(n_stream):
        j = i + 1
        while j < n_stream:
            distances[i][j] = np.sqrt(inners[i][i] + inners[j][j] - 2.0 * inners[i][j])
            distances[j][i] = distances[i][j]
            j = j + 1

    return distances, inners


def prototypes_extraction(streamlines, labels, inners):
    proto_dict = {}
    diag_norm2 = np.array([np.sqrt(inners[i][i]) ** 2 for i in range(np.size(inners, 0))])
    argmax_inner = np.sum(inners, axis=0) / diag_norm2
    for label in set(labels):
        if label != -1:
            argmax_label = np.NINF
            for i, la, s in zip(count(), labels, streamlines):
                if la == label and argmax_inner[i] > argmax_label:
                    argmax_label = argmax_inner[i]
                    proto_dict[str(label)] = s

    return list(itervalues(proto_dict)), list(proto_dict)


def get_prototypes(streamlines, fa, affine, lambda_a, lambda_b, lambda_w, lambda_m):
    if sys.maxsize == 2 ** 63 - 1:
        x = 64
        dt = np.dtype(np.float64)
    else:
        x = 32
        dt = np.dtype(np.float32)
    array_dim = (0.5 * x * len(streamlines) ** 2) / dt.itemsize
    if array_dim < (virtual_memory().total * 0.9):
        fvar, inners = varifolds_calc(streamlines, fa, affine, lambda_a=lambda_a, lambda_b=lambda_b, lambda_w=lambda_w,
                                      lambda_m=lambda_m)
    else:
        # IMPLEMENT MATRIX APPROXIMATION
        raise MemoryError(
            'Too many fibers in output of the initial segmentation (fVar not able to resolve the distance matrix)')

    labels = streamlines_clusters(fvar)

    prototypes, prototypes_keys = prototypes_extraction(streamlines, labels, inners)

    return prototypes, prototypes_keys, labels
