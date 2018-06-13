#!/usr/bin/env python

import numpy as np
import vtk
from builtins import range
from six import iteritems

try:
    import itertools.izip as zip
except ImportError:
    pass
from vtk.util import numpy_support as ns
from nibabel.streamlines.tractogram import Tractogram
from nibabel.streamlines.tck import TckFile as tck
from nibabel.streamlines.trk import TrkFile as tv


def read_tck(filename):
    tck_object = tck.load(filename)
    streamlines = tck_object.streamlines
    header = tck_object.header

    return streamlines, header


def read_trk(filename):
    trk_object = tv.load(filename)
    streamlines = trk_object.streamlines
    header = trk_object.header

    return streamlines, header


def read_vtk(filename):
    if filename.endswith('xml') or filename.endswith('vtp'):
        polydata_reader = vtk.vtkXMLPolyDataReader()
    else:
        polydata_reader = vtk.vtkPolyDataReader()

    polydata_reader.SetFileName(filename)
    polydata_reader.Update()

    polydata = polydata_reader.GetOutput()

    return vtkpolydata_to_tracts(polydata)


def vtkpolydata_to_tracts(polydata):
    result = {'lines': ns.vtk_to_numpy(polydata.GetLines().GetData()),
              'points': ns.vtk_to_numpy(polydata.GetPoints().GetData()), 'numberOfLines': polydata.GetNumberOfLines()}

    data = {}
    if polydata.GetPointData().GetScalars():
        data['ActiveScalars'] = polydata.GetPointData().GetScalars().GetName()
    if polydata.GetPointData().GetVectors():
        data['ActiveVectors'] = polydata.GetPointData().GetVectors().GetName()
    if polydata.GetPointData().GetTensors():
        data['ActiveTensors'] = polydata.GetPointData().GetTensors().GetName()

    for i in range(polydata.GetPointData().GetNumberOfArrays()):
        array = polydata.GetPointData().GetArray(i)
        np_array = ns.vtk_to_numpy(array)
        if np_array.ndim == 1:
            np_array = np_array.reshape(len(np_array), 1)
        data[polydata.GetPointData().GetArrayName(i)] = np_array

    result['pointData'] = data

    tracts, data = vtkpolydata_dictionary_to_tracts_and_data(result)
    return tracts, data


def vtkpolydata_dictionary_to_tracts_and_data(dictionary):
    dictionary_keys = {'lines', 'points', 'numberOfLines'}
    if not dictionary_keys.issubset(dictionary):
        raise ValueError("Dictionary must have the keys lines and points" + repr(dictionary))

    tract_data = {}
    tracts = []

    lines = np.asarray(dictionary['lines']).squeeze()
    points = dictionary['points']

    actual_line_index = 0
    number_of_tracts = dictionary['numberOfLines']
    original_lines = []
    for l in range(number_of_tracts):
        tracts.append(points[lines[actual_line_index + 1:actual_line_index + lines[actual_line_index] + 1]])
        original_lines.append(
            np.array(lines[actual_line_index + 1:actual_line_index + lines[actual_line_index] + 1], copy=True))
        actual_line_index += lines[actual_line_index] + 1

    if 'pointData' in dictionary:
        point_data_keys = [it[0] for it in iteritems(dictionary['pointData']) if isinstance(it[1], np.ndarray)]

        for k in point_data_keys:
            array_data = dictionary['pointData'][k]
            if k not in tract_data:
                tract_data[k] = [array_data[f] for f in original_lines]
            else:
                np.vstack(tract_data[k])
                tract_data[k].extend([array_data[f] for f in original_lines[-number_of_tracts:]])

    return tracts, tract_data


def save_vtk(filename, tracts, lines_indices=None):
    lengths = [len(p) for p in tracts]
    line_starts = ns.numpy.r_[0, ns.numpy.cumsum(lengths)]
    if lines_indices is None:
        lines_indices = [ns.numpy.arange(length) + line_start for length, line_start in zip(lengths, line_starts)]

    ids = ns.numpy.hstack([ns.numpy.r_[c[0], c[1]] for c in zip(lengths, lines_indices)])
    vtk_ids = ns.numpy_to_vtkIdTypeArray(ids, deep=True)

    cell_array = vtk.vtkCellArray()
    cell_array.SetCells(len(tracts), vtk_ids)
    points = ns.numpy.vstack(tracts).astype(ns.get_vtk_to_numpy_typemap()[vtk.VTK_DOUBLE])
    points_array = ns.numpy_to_vtk(points, deep=True)

    poly_data = vtk.vtkPolyData()
    vtk_points = vtk.vtkPoints()
    vtk_points.SetData(points_array)
    poly_data.SetPoints(vtk_points)
    poly_data.SetLines(cell_array)

    poly_data.BuildCells()

    if filename.endswith('.xml') or filename.endswith('.vtp'):
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetDataModeToBinary()
    else:
        writer = vtk.vtkPolyDataWriter()
        writer.SetFileTypeToBinary()

    writer.SetFileName(filename)
    if hasattr(vtk, 'VTK_MAJOR_VERSION') and vtk.VTK_MAJOR_VERSION > 5:
        writer.SetInputData(poly_data)
    else:
        writer.SetInput(poly_data)
    writer.Write()


def save_tck(filename, tracts, header, fuzzy):
    tractogram = Tractogram(tracts, affine_to_rasmm=np.eye(4), data_per_streamline=fuzzy)
    if not header:
        header = tck.create_empty_header()
    tck_obj = tck(tractogram, header)
    tck_obj.save(filename)


def save_trk(filename, tracts, header, fuzzy):
    tractogram = Tractogram(tracts, affine_to_rasmm=np.eye(4), data_per_streamline=fuzzy)
    if not header:
        header = tv.create_empty_header()
    trk_obj = tv(tractogram, header)
    trk_obj.save(filename)
