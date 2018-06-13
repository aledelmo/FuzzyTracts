#!/usr/bin/env python

import glob
import os
import shutil
import tempfile

import numpy as np
from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, traits

from .input_output import read_tck, read_trk, read_vtk
from .streamlines_operations import streamlines_merge, streamlines_mapvolume
from .utilis import pipe


class TractQuerierInputSpec(CommandLineInputSpec):
    input_query = traits.File(desc='Input Query file', exists=False, mandatory=True, argstr='-q %s',
                              copy_file=False)
    input_atlas = traits.File(desc="Input Atlas volume", exists=False, mandatory=True, argstr="-a %s", copy_file=False)
    input_tractography = traits.File(desc="Input Tractography", exists=False, mandatory=True, argstr="-t %s",
                                     copy_file=False)
    out_prefix = traits.Str('query', des="prefix for the results", mandatory=False, argstr="-o %s", usedefault=True)


class TractQuerierOutputSpec(TraitedSpec):
    output_queries = traits.List(exists=True, desc='resulting query files')


class TractQuerier(CommandLine):
    _cmd = 'tract_querier'
    input_spec = TractQuerierInputSpec
    output_spec = TractQuerierOutputSpec

    def _format_arg(self, name, spec, value):
        return super(TractQuerier, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['output_queries'] = glob.glob(os.path.join(os.getcwd(),
                                                           self.inputs.out_prefix +
                                                           '*.trk')
                                              )
        return outputs


def wmql(in_file, object_list, atlas, query):
    folder = tempfile.mkdtemp()
    cmd = 'tract_querier -a ' + atlas + ' -t ' + in_file + ' -q ' + query + ' -o ' + \
          os.path.join(folder, 'temp.trk')
    pipe(cmd)
    # current = os.getcwd()
    # os.chdir(folder)
    # tract_querier = TractQuerier()
    # tract_querier.inputs.input_query = query
    # tract_querier.inputs.input_atlas = atlas
    # tract_querier.inputs.input_tractography = in_file
    # tract_querier.inputs.out_prefix = 'temp'
    # tract_querier.cmdline
    # os.chdir(current)

    directory = os.listdir(folder)
    for cones_set in object_list:
        streamlines = []
        for in_file in directory:
            in_file = os.path.join(folder, in_file)
            if cones_set.ID.lower() in in_file.lower():
                new_streamlines, new_header = read_trk(in_file)
                streamlines, header = streamlines_merge(streamlines, new_streamlines, hdr=new_header)
        cones_set.put_streamlines(streamlines)

    shutil.rmtree(folder)

    return header


def fuzzy(in_file, object_list, file_extension):
    if file_extension == '.tck':
        streamlines, header = read_tck(in_file)
    if file_extension == '.trk':
        streamlines, header = read_trk(in_file)
    if file_extension == '.vtk' or file_extension == '.xml' or file_extension == '.vtp':
        streamlines = read_vtk(in_file)[0]

    for cones_set in object_list:
        mask, affine = cones_set.get_intersection()
        mapping = streamlines_mapvolume(streamlines, mask, affine)
        mapping = [np.array(m) for m in mapping]
        mapping = np.asarray([np.amax(m) for m in mapping])
        segmentation = [streamlines[i] for i, m in enumerate(mapping) if m != 0]
        cones_set.put_streamlines(segmentation)

    if file_extension == '.tck' or file_extension == '.trk':
        return header
    else:
        return None
