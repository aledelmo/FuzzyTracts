# Fuzzy Tracts

Streamlines selection with unsupervised clustering and fuzzy functional varifolds.

![acs](https://i.imgur.com/GCzmyiO.jpg)

Polylines segmentation combining spatial fuzzy logic and clustering techniques. Streamlines coherence is evaluated using
spatial fuzzy sets based on positional a priori. Similarity matrix between streamlines is built using the 
following varifold formulation:

![equation](http://www.sciweavers.org/tex2img.php?eq=%24%5Clangle%20C_X%2CC_Y%5Crangle%3De%5E%7B-%5Cfrac%7B%7C%7Cf_a-t_a%20%7C%7C%5E2%7D%7B%5Clambda%5E2_a%7D%7D%2Ae%5E%7B-%5Cfrac%7B%7C%7Cf_b-t_b%20%7C%7C%5E2%7D%7B%5Clambda%5E2_b%7D%7D%2A%5Csum_%7Bi%7D%5E%7BN%7D%5Csum_%7Bi%7D%5E%7BN%7D%5Be%5E%7B-%5Cfrac%7B%7C%7Cx_i-y_i%20%7C%7C%5E2%7D%7B%5Clambda%5E2_w%7D%7D%2Ae%5E%7B-%5Cfrac%7B%7C%7CFA_i-FA_j%20%7C%7C%5E2%7D%7B%5Clambda%5E2_m%7D%7D%2A%28%5Cfrac%7B%5Calpha%5ET_i%2A%5Cbeta_j%7D%7B%7C%5Calpha_i%7C%2A%7C%5Cbeta_j%7C%7D%29%5E2%2A%7C%5Calpha_i%7C%2A%7C%5Cbeta_j%7C%5D%24&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0[/img])

DBSCAN is used to cluster streamlines and to select fibers prototypes. Prototypes are finally selected/discarded based
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