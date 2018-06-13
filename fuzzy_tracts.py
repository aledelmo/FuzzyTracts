#!/usr/bin/env python

from __future__ import division

import argparse
import json
import os.path
import sys
import time
import warnings

from six import iteritems

import logic

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

__author__ = 'Alessandro Delmonte'
__email__ = 'delmonte.ale92@gmail.com'


def main():
    in_file, use_wmql, r, fa, show = setup()

    with open(os.path.join(os.path.dirname(__file__), 'config.json')) as config_file:
        config = json.load(config_file)

    t = check_threshold(config['conf']['threshold'])
    atlas = check_nii(config['conf']['atlas_path'])
    legend = config['conf']['labels_legend']
    if use_wmql:
        query = check_qry(config['conf']['WMQL_query'])
    if r:
        r = check_resample(config['conf']['resampling'])
    if fa:
        fa = check_nii(config['conf']['FA_map'])

    folder = config['conf']['diffusion_folder']
    out_ext = check_format(config['conf']['out_format'])

    print ('Initializing Fuzzy Sets...')
    fuzzy_set_list = fuzzy_cones(config, atlas, legend, show)

    print ('Computing initial segmentation...')
    _, file_extension = os.path.splitext(in_file)
    if use_wmql:
        if file_extension == '.tck':
            if folder:
                in_file = logic.tck_to_trk(in_file, check_fld(folder))
            else:
                raise argparse.ArgumentTypeError(
                    "When using WMQL with a .tck file it is mandatory to specify the subject folder")

        header = logic.wmql(in_file, fuzzy_set_list, atlas, query)
    else:
        header = logic.fuzzy(in_file, fuzzy_set_list, file_extension)

    if r:
        print ('Resampling fibers...')
        for fuzzy_set in fuzzy_set_list:
            fuzzy_set.downsampling(r)

    if fa:
        print ('Loading FA map...')
        fa, affine = logic.load_nifti(fa)
    else:
        print ('Extracting FA map from diffusion data...')
        if folder:
            fa, affine = logic.dwi2fa(check_fld(folder))
        else:
            raise argparse.ArgumentTypeError(
                'Impossible to extract diffusion information (no FA map or subject folder given.)')

    print ('Fibers clustering and prototypes selection...')
    for fuzzy_set in fuzzy_set_list:
        fuzzy_set.clustering(fa, affine)

    print ('Final segmentation...')
    for fuzzy_set in fuzzy_set_list:
        fuzzy_set.segmentation(t)

    print ('Saving results...')
    for fuzzy_set in fuzzy_set_list:
        filename = fuzzy_set.ID + '_' + str(t) + '_FuzzyTracts'
        if out_ext == 'tck':
            logic.save_tck(filename + '.' + out_ext, fuzzy_set.final, header, fuzzy_set.fuzzy_per_line)
        if out_ext == 'trk':
            logic.save_trk(filename + '.' + out_ext, fuzzy_set.final, header, fuzzy_set.fuzzy_per_line)
        if out_ext == 'vtk' or out_ext == 'xml' or out_ext == 'vtp':
            logic.save_vtk(filename + '.' + out_ext, fuzzy_set.final)


def setup():
    parser = argparse.ArgumentParser()
    parser.add_argument('Input_Tractogram', help='Name of the input tractography file', type=check_ext)
    parser.add_argument('-w', '--WMQL', help='Initialization using WMQL', action='store_true')
    parser.add_argument('-fa', '--FA', help='Precomputed FA map', action='store_true')
    parser.add_argument('-r', '--resample', help='Downsampling streamlines (might improve computational time)',
                        action='store_true')
    parser.add_argument('-s', '--show', help='Use the visualization tools.', action='store_true')

    args = parser.parse_args()

    return args.Input_Tractogram, args.WMQL, args.resample, args.FA, args.show


def check_ext(value):
    filename, file_extension = os.path.splitext(value)
    if file_extension in ('.vtk', '.xml', '.vtp', '.tck', '.trk'):
        return value
    else:
        raise argparse.ArgumentTypeError("Invalid file extension (file format supported: tck, trk, vtk): %s" % value)


def check_threshold(value):
    try:
        t = float(value)
        if 0 < t <= 1:
            return t
        else:
            raise argparse.ArgumentTypeError("Invalid threshold (must be between 0 and 1): %r" % value)
    except ValueError:
        raise argparse.ArgumentTypeError("Invalid threshold value: %r" % value)


def check_nii(value):
    if value.split('.')[1] in ('nii', 'nii.gz'):
        return os.path.abspath(value)
    else:
        raise argparse.ArgumentTypeError("Invalid file extension (support for .nii and .nii.gz): %s" % value)


def check_qry(value):
    filename, file_extension = os.path.splitext(value)
    if file_extension in ('.qry'):
        return os.path.abspath(value)
    else:
        raise argparse.ArgumentTypeError("Invalid query format (.qry): %s" % value)


def check_resample(value):
    try:
        t = float(value)
        if 0 < t < 100:
            return t
        else:
            raise argparse.ArgumentTypeError("Invalid threshold (must be between 0%% and 100%%): %r" % value)
    except ValueError:
        raise argparse.ArgumentTypeError("Invalid threshold value: %r" % value)


def check_fld(value):
    if os.path.isdir(value):
        return value
    else:
        raise argparse.ArgumentTypeError("Invalid diffusion folder: %r", value)


def check_format(value):
    if value.lower() in ('vtk', 'xml', 'vtp', 'tck', 'trk'):
        return value.lower()
    else:
        raise argparse.ArgumentTypeError("Invalid file extension (file format supported: tck, trk, vtk): %r" % value)


def fuzzy_cones(config, atlas, legend, show):
    parcellation, affine = logic.load_nifti(atlas)
    labels_dict = logic.load_label_legend(os.path.abspath(legend))

    sets = []
    for key, value in iteritems(config):
        if key != 'conf':
            set_fuzzy = logic.SpatialFuzzySets(key, value, parcellation, affine, labels_dict, show)
            sets.append(set_fuzzy)

    return sets


if __name__ == "__main__":
    start = time.time()
    main()
    print ('Total runtime: %.2f m.' % ((time.time() - start) / 60))
    sys.exit(0)
