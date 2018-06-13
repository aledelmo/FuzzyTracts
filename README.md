# Fuzzy Tracts

Fuzzy Tracts is a tractograms segmentation tool that combines spatial information with a fuzzy sets approach!

  - Define the desired bundles.
  - Choose your parameters (or use ours!)
  - Set the segmentation threshold that best fits you.
  - Done !!!

For more information visit: <https://hal.archives-ouvertes.fr/hal-01744267/document> <br />
Fuzzy Tracts is open source with a [public repository][gitlab] on GitLab.

Data from:
[![HCP](https://wiki.humanconnectome.org/download/attachments/589826/global.logo?version=2&modificationDate=1326402585274&api=v2)](http://www.humanconnectomeproject.org/)

## Installation and Usage.

Enter in your freshly downloaded Fuzzy Tracts folder:
```sh
$ cd fuzzy_tracts
```

Fuzzy tracts requires some dependencies to run smoothly. To install everything needed run from a terminal:
```sh
$ pip install -r requirements.txt
```

Mandatory Python Dependencies: **Dipy**, **scikit-learn**, **scipy** <br />
Highly suggested Python Dependencies: **joblib** <br />
Python3 Compatibility Dependencies: **future** <br />
Optional Python Dependencies: **Nipype**, **VTK** <br />
Optional Software Dependencies:
    **MRTrix3**: <http://www.mrtrix.org/>
    **WMQL**: <http://tract-querier.readthedocs.io/en/latest/#>

Configure the tool with your favourite text editor:
```sh
$ vim config.json
```

Launch the tool:
```sh
$ python fuzzy_tracts.py input_tractogram
```


## Updates

### New Features !!!

  - Improved perfomance.
  - Added native support for tck and vtk.

### Todos:

 - fVar matrix approximation
 - Adding IFOF.
 


## Who's behind this!

Fuzzy Tracts is currently being developed at [LTCI] lab in Télécom ParisTech.

| Contributor | Personal Site |
| ------ | ------ |
| Alessandro Delmonte | Work in progress =) |
| Isabelle Bloch | https://perso.telecom-paristech.fr/bloch/ |
| Dominique Hasboun | ...|
| Corentin Mercier | https://perso.telecom-paristech.fr/comercier/ |
| Johan Pallud | http://www.ch-sainte-anne.fr/Offres-de-soins/Neuro-Sainte-Anne/Neuro-oncologie/Johan-PALLUD |
| Pietro Gori | https://perso.telecom-paristech.fr/pgori/ |


### Development.

The tool is constantly updated with new functionalities.

Want to contribute? <br />
Please contact Alessandro (alessandro.delmonte@telecom-paristech.fr) or Pietro (pietro.gori@telecom-paristech.fr)



License.
----

Apache License 2.0

[//]: #
   [gitlab]: <https://gitlab.telecom-paristech.fr/equipe-images/biomed/Segmentation_Tractography>
   [LTCI]: <https://ltci.telecom-paristech.fr/>