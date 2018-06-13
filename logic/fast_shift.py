#!/usr/bin/env python

import numpy as np


def fast_shift3d(in_array, positions, padding=0):
    out_array = np.zeros_like(in_array)

    pos0 = positions[0]
    pos1 = positions[1]
    pos2 = positions[2]

    if pos0 >= 0 and pos1 >= 0 and pos2 >= 0:
        if pos0 == 0 and pos1 > 0 and pos2 > 0:
            out_array[:, pos1:, pos2:] = in_array[:, :-pos1, :-pos2]
            out_array[:, :pos1, :pos2] = padding
        elif pos0 == 0 and pos1 > 0 and pos2 == 0:
            out_array[:, pos1:, :] = in_array[:, :-pos1, :]
            out_array[:, :pos1, :] = padding
        elif pos0 == 0 and pos1 == 0 and pos2 > 0:
            out_array[:, :, pos2:] = in_array[:, :, :-pos2]
            out_array[:, :, :pos2] = padding
        elif pos0 == 0 and pos1 == 0 and pos2 == 0:
            out_array[:, :, :] = in_array[:, :, :]
            out_array[:, :, :] = padding
        elif pos0 > 0 and pos1 == 0 and pos2 > 0:
            out_array[pos0:, :, pos2:] = in_array[:-pos0, :, :-pos2]
            out_array[:pos0, :, :pos2] = padding
        elif pos0 > 0 and pos1 == 0 and pos2 == 0:
            out_array[pos0:, :, :] = in_array[:-pos0, :, :]
            out_array[:pos0, :, :] = padding
        elif pos0 > 0 and pos1 > 0 and pos2 == 0:
            out_array[pos0:, pos1:, :] = in_array[:-pos0, :-pos1, :]
            out_array[:pos0, :pos1, :] = padding
        else:
            out_array[pos0:, pos1:, pos2:] = in_array[:-pos0, :-pos1, :-pos2]
            out_array[:pos0, :pos1, :pos2] = padding
    elif pos0 >= 0 and pos1 >= 0 > pos2:
        if pos0 == 0 and pos1 != 0:
            out_array[:, pos1:, :pos2] = in_array[:, :-pos1, -pos2:]
            out_array[:, :pos1, pos2:] = padding
        elif pos0 != 0 and pos1 == 0:
            out_array[pos0:, :, :pos2] = in_array[:-pos0, :, -pos2:]
            out_array[:pos0, :, pos2:] = padding
        elif pos0 == 0 and pos1 == 0:
            out_array[:, :, :pos2] = in_array[:, :, -pos2:]
            out_array[:, :, pos2:] = padding
        else:
            out_array[pos0:, pos1:, :pos2] = in_array[:-pos0, :-pos1, -pos2:]
            out_array[:pos0, :pos1, pos2:] = padding
    elif pos0 >= 0 > pos1 and pos2 >= 0:
        if pos0 == 0 and pos2 != 0:
            out_array[:, :pos1, pos2:] = in_array[:, -pos1:, :-pos2]
            out_array[:, pos1:, :pos2] = padding
        elif pos0 != 0 and pos2 == 0:
            out_array[pos0:, :pos1, :] = in_array[:-pos0, -pos1:, :]
            out_array[:pos0, pos1:, :] = padding
        elif pos0 == 0 and pos2 == 0:
            out_array[:, :pos1, :] = in_array[:-pos0, -pos1:, :]
            out_array[:, pos1:, :] = padding
        else:
            out_array[pos0:, :pos1, pos2:] = in_array[:-pos0, -pos1:, :-pos2]
            out_array[:pos0, pos1:, :pos2] = padding
    elif pos0 >= 0 > pos1 and pos2 < 0:
        if pos0 == 0:
            out_array[:, :pos1, :pos2] = in_array[:, -pos1:, -pos2:]
            out_array[:, pos1:, pos2:] = padding
        else:
            out_array[pos0:, :pos1, :pos2] = in_array[:-pos0, -pos1:, -pos2:]
            out_array[:pos0, pos1:, pos2:] = padding
    elif pos0 < 0 <= pos1 and pos2 >= 0:
        if pos1 == 0 and pos2 != 0:
            out_array[:pos0, :, pos2:] = in_array[-pos0:, :, :-pos2]
            out_array[pos0:, :, :pos2] = padding
        elif pos1 != 0 and pos2 == 0:
            out_array[:pos0, pos1:, :] = in_array[-pos0:, :-pos1, :]
            out_array[pos0:, :pos1, :] = padding
        elif pos1 == 0 and pos2 == 0:
            out_array[:pos0, :, :] = in_array[-pos0:, :, :]
            out_array[pos0:, :, :] = padding
        else:
            out_array[:pos0, pos1:, pos2:] = in_array[-pos0:, :-pos1, :-pos2]
            out_array[pos0:, :pos1, :pos2] = padding
    elif pos0 < 0 <= pos1 and pos2 < 0:
        if pos1 == 0:
            out_array[:pos0, :, :pos2] = in_array[-pos0:, :, -pos2:]
            out_array[pos0:, :, pos2:] = padding
        else:
            out_array[:pos0, pos1:, :pos2] = in_array[-pos0:, :-pos1, -pos2:]
            out_array[pos0:, :pos1, pos2:] = padding
    elif pos0 < 0 and pos1 < 0 <= pos2:
        if pos2 == 0:
            out_array[:pos0, :pos1, :] = in_array[-pos0:, -pos1:, :]
            out_array[pos0:, pos1:, :] = padding
        else:
            out_array[:pos0, :pos1, pos2:] = in_array[-pos0:, -pos1:, :-pos2]
            out_array[pos0:, pos1:, :pos2] = padding
    else:
        out_array[:pos0, :pos1, :pos2] = in_array[-pos0:, -pos1:, -pos2:]
        out_array[pos0:, pos1:, pos2:] = padding

    return out_array


def fast_shift3d_parallel(in_array, positions, padding=0):
    cone = np.zeros(in_array.shape)

    for position in positions:
        pos0 = position[0]
        pos1 = position[1]
        pos2 = position[2]

        out_array = np.zeros_like(in_array)

        if pos0 >= 0 and pos1 >= 0 and pos2 >= 0:
            if pos0 == 0 and pos1 > 0 and pos2 > 0:
                out_array[:, pos1:, pos2:] = in_array[:, :-pos1, :-pos2]
                out_array[:, :pos1, :pos2] = padding
            elif pos0 == 0 and pos1 > 0 and pos2 == 0:
                out_array[:, pos1:, :] = in_array[:, :-pos1, :]
                out_array[:, :pos1, :] = padding
            elif pos0 == 0 and pos1 == 0 and pos2 > 0:
                out_array[:, :, pos2:] = in_array[:, :, :-pos2]
                out_array[:, :, :pos2] = padding
            elif pos0 == 0 and pos1 == 0 and pos2 == 0:
                out_array[:, :, :] = in_array[:, :, :]
                out_array[:, :, :] = padding
            elif pos0 > 0 and pos1 == 0 and pos2 > 0:
                out_array[pos0:, :, pos2:] = in_array[:-pos0, :, :-pos2]
                out_array[:pos0, :, :pos2] = padding
            elif pos0 > 0 and pos1 == 0 and pos2 == 0:
                out_array[pos0:, :, :] = in_array[:-pos0, :, :]
                out_array[:pos0, :, :] = padding
            elif pos0 > 0 and pos1 > 0 and pos2 == 0:
                out_array[pos0:, pos1:, :] = in_array[:-pos0, :-pos1, :]
                out_array[:pos0, :pos1, :] = padding
            else:
                out_array[pos0:, pos1:, pos2:] = in_array[:-pos0, :-pos1, :-pos2]
                out_array[:pos0, :pos1, :pos2] = padding
        elif pos0 >= 0 and pos1 >= 0 > pos2:
            if pos0 == 0 and pos1 != 0:
                out_array[:, pos1:, :pos2] = in_array[:, :-pos1, -pos2:]
                out_array[:, :pos1, pos2:] = padding
            elif pos0 != 0 and pos1 == 0:
                out_array[pos0:, :, :pos2] = in_array[:-pos0, :, -pos2:]
                out_array[:pos0, :, pos2:] = padding
            elif pos0 == 0 and pos1 == 0:
                out_array[:, :, :pos2] = in_array[:, :, -pos2:]
                out_array[:, :, pos2:] = padding
            else:
                out_array[pos0:, pos1:, :pos2] = in_array[:-pos0, :-pos1, -pos2:]
                out_array[:pos0, :pos1, pos2:] = padding
        elif pos0 >= 0 and pos1 < 0 <= pos2:
            if pos0 == 0 and pos2 != 0:
                out_array[:, :pos1, pos2:] = in_array[:, -pos1:, :-pos2]
                out_array[:, pos1:, :pos2] = padding
            elif pos0 != 0 and pos2 == 0:
                out_array[pos0:, :pos1, :] = in_array[:-pos0, -pos1:, :]
                out_array[:pos0, pos1:, :] = padding
            elif pos0 == 0 and pos2 == 0:
                out_array[:, :pos1, :] = in_array[:-pos0, -pos1:, :]
                out_array[:, pos1:, :] = padding
            else:
                out_array[pos0:, :pos1, pos2:] = in_array[:-pos0, -pos1:, :-pos2]
                out_array[:pos0, pos1:, :pos2] = padding
        elif pos0 >= 0 > pos1 and pos2 < 0:
            if pos0 == 0:
                out_array[:, :pos1, :pos2] = in_array[:, -pos1:, -pos2:]
                out_array[:, pos1:, pos2:] = padding
            else:
                out_array[pos0:, :pos1, :pos2] = in_array[:-pos0, -pos1:, -pos2:]
                out_array[:pos0, pos1:, pos2:] = padding
        elif pos0 < 0 <= pos1 and pos2 >= 0:
            if pos1 == 0 and pos2 != 0:
                out_array[:pos0, :, pos2:] = in_array[-pos0:, :, :-pos2]
                out_array[pos0:, :, :pos2] = padding
            elif pos1 != 0 and pos2 == 0:
                out_array[:pos0, pos1:, :] = in_array[-pos0:, :-pos1, :]
                out_array[pos0:, :pos1, :] = padding
            elif pos1 == 0 and pos2 == 0:
                out_array[:pos0, :, :] = in_array[-pos0:, :, :]
                out_array[pos0:, :, :] = padding
            else:
                out_array[:pos0, pos1:, pos2:] = in_array[-pos0:, :-pos1, :-pos2]
                out_array[pos0:, :pos1, :pos2] = padding
        elif pos0 < 0 <= pos1 and pos2 < 0:
            if pos1 == 0:
                out_array[:pos0, :, :pos2] = in_array[-pos0:, :, -pos2:]
                out_array[pos0:, :, pos2:] = padding
            else:
                out_array[:pos0, pos1:, :pos2] = in_array[-pos0:, :-pos1, -pos2:]
                out_array[pos0:, :pos1, pos2:] = padding
        elif pos0 < 0 and pos1 < 0 <= pos2:
            if pos2 == 0:
                out_array[:pos0, :pos1, :] = in_array[-pos0:, -pos1:, :]
                out_array[pos0:, pos1:, :] = padding
            else:
                out_array[:pos0, :pos1, pos2:] = in_array[-pos0:, -pos1:, :-pos2]
                out_array[pos0:, pos1:, :pos2] = padding
        else:
            out_array[:pos0, :pos1, :pos2] = in_array[-pos0:, -pos1:, -pos2:]
            out_array[pos0:, pos1:, pos2:] = padding
        np.maximum(cone, out_array, out=cone)
    return cone
