#!/usr/bin/env python

import os.path

try:
	from setuptools import setup
except ImportError:
	from distutils.core import setup


def read(fname):
	return open(os.path.join(os.path.dirname(__file__), fname)).read()


if __name__ == "__main__":
	setup(
		name='Fuzzy Tracts',
		version='1.0',
		packages=['logic'],
		scripts=['fuzzy_tracts.py'],
		url='https://hal.archives-ouvertes.fr/hal-01744267/document',
		license=read('license..rst'),
		author='Alessandro Delmonte',
		author_email='delmonte.ale92@gmail.com',
		description='Fuzzy Tracts - A Tractogram Segmentation Tool',
		requires=['dipy == 0.12.0', 'scikit-learn == 0.19.0', 'joblib == 0.11', 'future == 0.16.0', 'nipype == 0.13.1',
		          'vtk == 7.1.1'],
		long_description=read('README.md'),
		classifiers=[
			'Intended Audience :: Science/Research',
			'Programming Language :: Python',
			'Topic :: Scientific/Engineering',
			'Operating System :: Unix'
		],
	)
