# FAST FUSION OF HYPERSPECTRAL AND MULTISPECTRAL IMAGES : A TUCKER APPROXIMATION APPROACH

Copyright (c) 2022 Clemence Prevost, Pierre Chainais, Remy Boyer <br>
Contact: ```clemence.prevost@univ-lille.fr```

This software reproduces the results from the following:
```
@unpublished{prevost:hal-03617759,
  TITLE = {{FAST FUSION OF HYPERSPECTRAL AND MULTISPECTRAL IMAGES : A TUCKER APPROXIMATION APPROACH}},
  AUTHOR = {Pr{\'e}vost, Cl{\'e}mence and Chainais, Pierre and Boyer, Remy},
  URL = {https://hal.archives-ouvertes.fr/hal-03617759},
  NOTE = {working paper or preprint},
  YEAR = {2022},
  MONTH = Mar,
  KEYWORDS = {hyperspectral super-resolution ; data fusion ; low-rank tensor factorizations ; recovery ; least-squares problem},
  PDF = {https://hal.archives-ouvertes.fr/hal-03617759/file/icip.pdf},
  HAL_ID = {hal-03617759},
  HAL_VERSION = {v1},
}
```

<br><br>
[Link to the project](https://github.com/cprevost4/enhanced_scott)

## Content

 - /metrics : contains the metrics used to assess the quality of the reconstruction
 - /src : contains helpful files to run the demos

## Minimal requirements

 In order to run the demo files, you will need to:
 - Download and install Tensorlab 3.0: https://www.tensorlab.net
 - Download the dataset from the AVIRIS website: https://aviris.jpl.nasa.gov/dataportal/
 
  Please quote the corresponding paper if you decide to use these codes.

 ## How it works
 
 This software reproduces the figures and tables contained in the paper. You can play with the two datasets.
 
 ### Generate coupled tensor model
 
 Real datasets are used. First, the HSI and MSI are generated following Wald's protocol. Then, white Gaussian noise is added to the observations.

 ### Run algorithms
 
 In ```fusion_isabella.m``` and ```fusion_lockwood.m```, we showcase the performance of:
  - The proposed algorithm, with [2,2] and [4,4] block-patterns
  - SCOTT, for comparison

### Plot results 
The metrics and computation time are then displayed in a table.
Slices of the reference and reconstructions are plotted in a figure.
The color scale is generated according to the image's spectral support.

<img src="img/illu.png?raw=true"/>

## Available demos

They are available in the ```/demos``` folder.
The table below summarized what does what

| Name                       | Content                                           |
|----------------------------|---------------------------------------------------|
| ```fusion_isabella.m```    | Simulations for Isabella Lake dataset             |
| ```fusion_lockwood.m```    | Simulations for Lockwood dataset                  |
| ```choice_ranks.m```       | plots R-SNR as a function of the multilinear ranks|


