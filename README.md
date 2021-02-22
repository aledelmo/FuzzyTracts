# Fuzzy Tracts

Streamlines selection with unsupervised clustering and fuzzy functional varifolds.

![acs](https://i.imgur.com/GCzmyiO.jpg)

Polylines segmentation combining spatial fuzzy logic and clustering techniques. Streamlines coherence evaluated using
spatial fuzzy sets based on positional a priori. Similarity matrix between streamlines built using the 
following varifold formulation:

![varifold](https://i.imgur.com/sSEo2Bc.png)

DBSCAN for streamlines clustering and fibers prototypes selection. Prototypes selected/discarded based
on their fuzzy score.

Method developed for the recognition of neural connections from MRI scans. Data from [Human Connectome Project](http://www.humanconnectomeproject.org/)

If you use this work please cite:
>**[White Matter Multi-resolution Segmentation using Fuzzy Theory](https://hal.archives-ouvertes.fr/hal-01983010/document)** *Alessandro Delmonte, Corentin Mercier, Johan Pallud, Isabelle Bloch and Pietro Gori.* ISBI 2019 - IEEE International Symposium on
Biomedical Imaging, April 2019, Venice, IT.
> 
>**[Segmentation of White Matter Tractograms Using Fuzzy Spatial Relations](https://hal.archives-ouvertes.fr/hal-01744267/document)** *Alessandro Delmonte, Isabelle Bloch, Dominique Hasboun, Corentin Mercier, Johan Pallud and Pietro Gori.* OHBM 2018 -
> Organization for Human Brain Mapping, June 2018, Singapore, SG.

## Installation and Usage.

The project may use external tools. Check [WMQL](http://tract-querier.readthedocs.io/en/latest/#) and [MRTrix3](http://www.mrtrix.org/) guides for installation.

How to run:
```shell
$ cd FuzzyTracts
$ pip install -r requirements.txt
$ vim config.json
$ python fuzzy_tracts.py <input_tractogram>
```

Run the provided example with
```shell

$ python fuzzy_tracts.py test_data/input_test.vtk -fa -s -r
```

## Contacts

For any inquiries please contact: 
[Alessandro Delmonte](https://aledelmo.github.io) @ [alessandro.delmonte@institutimagine.org](mailto:alessandro.delmonte@institutimagine.org)

## License

This project is licensed under the [Apache License 2.0](LICENSE) - see the [LICENSE](LICENSE) file for
details