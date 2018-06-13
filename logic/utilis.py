#!/usr/bin/env python

import os.path
import shlex
import shutil
import tempfile
from subprocess import Popen, PIPE, call

import nibabel as nib
import nipype.interfaces.mrtrix as mrt
import nipype.interfaces.mrtrix3 as mrt3
import nipype.pipeline.engine as pe
import numpy as np
from builtins import int
from dipy.core.gradients import gradient_table
from dipy.io import read_bvals_bvecs
from dipy.reconst.dti import fractional_anisotropy, TensorModel


def load_nifti(fname):
    img = nib.load(fname)
    data = img.get_data()
    affine = img.get_affine()

    return data, affine


def load_label_legend(filepath):
    labels_dict = {}
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            labels_dict[line.split()[1]] = int(line.split()[0])

    return labels_dict


def dwi2fa(folder):
    try:
        tensor = pe.Node(interface=mrt3.FitTensor(args='-quiet'), name='tensor')
        data_path = os.path.join(os.path.abspath(folder), 'T1w/Diffusion/data.nii.gz')
        tensor.inputs.in_file = data_path
        mask_path = os.path.join(os.path.abspath(folder), 'T1w/Diffusion/nodif_brain_mask.nii.gz')
        tensor.inputs.in_mask = mask_path
        fbval = os.path.join(os.path.abspath(folder), 'T1w/Diffusion/bvals')
        fbvec = os.path.join(os.path.abspath(folder), 'T1w/Diffusion/bvecs')
        tensor.inputs.grad_fsl = (fbvec, fbval)

        folder_tmp = tempfile.mkdtemp()

        fa = pe.Node(
            interface=mrt3.TensorMetrics(out_fa=os.path.join(os.path.abspath(folder_tmp), 'FA.nii'), args='-quiet'),
            name='FA')

        workflow = pe.Workflow(name='FA_workflow')
        workflow.base_dir = folder_tmp
        workflow.connect(tensor, 'out_file', fa, 'in_file')
        workflow.run()

        fa, affine = load_nifti(os.path.join(os.path.abspath(folder_tmp), 'FA.nii'))

        shutil.rmtree(folder_tmp)

    except:
        data, affine = load_nifti(os.path.join(os.path.abspath(folder), 'T1w/Diffusion/data.nii.gz'))

        fbval = os.path.join(os.path.abspath(folder), 'T1w/Diffusion/bvals')
        fbvec = os.path.join(os.path.abspath(folder), 'T1w/Diffusion/bvecs')
        bvals, bvecs = read_bvals_bvecs(fbval, fbvec)
        gtab = gradient_table(bvals, bvecs)

        tenmodel = TensorModel(gtab)
        tenfit = tenmodel.fit(data)
        fa = fractional_anisotropy(tenfit.evals)
        fa[np.isnan(fa)] = 0

    return fa, affine


def pipe(cmd, verbose=False, topipe=False):
    if verbose:
        print ('Executing command: ' + str(cmd) + '\n')
    if topipe:
        p = Popen(cmd, shell=True, stdin=None, stdout=PIPE, stderr=PIPE)
        [stdout, stderr] = p.communicate()
        return stdout, stderr
    else:
        cmd = shlex.split(cmd)
        call(cmd, shell=False, stdin=None, stdout=None, stderr=None)


def tck_to_trk(filename, folder):
    name, extension = os.path.splitext(filename)

    tck2trk = mrt.MRTrix2TrackVis()
    tck2trk.inputs.in_file = filename
    tck2trk.inputs.image_file = os.path.join(os.path.abspath(folder), 'T1w/T1w_acpc_dc_restore_1.25.nii.gz')
    tck2trk.inputs.out_filename = name + '.trk'
    tck2trk.run()

    return name + '.trk'
