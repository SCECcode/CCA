# CCA

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)


## Description

#CCA

CCA-06 model is a 3D velocity model for Central California.  CCA-06 was 
produced through tomographic inversion, using the USGS Bay Area (v8.3.0),
CVM-S4.26, and Lin-Thurber models as starting models. 6 iterations were 
performed on a 500m-resolution mesh, down to a minimum Vs of 900 m/s,
to generate the final model.  This model includes an optional Ely-Jordan GTL.

## Installation

This package is intended to be installed as part of the UCVM framework,
version 25.x or higher.

## Library

The library ./lib/libcca.a may be statically linked into any
user application. Also, if your system supports dynamic linking,
you will also have a ./lib/libcca.so file that can be used
for dynamic linking. The header file defining the API is located
in ./include/cca.h.

## Contact the authors

If you would like to contact the authors regarding this software,
please e-mail software@scec.org. Note this e-mail address should
be used for questions regarding the software itself (e.g. how
do I link the library properly?). Questions regarding the model's
science (e.g. on what paper is the CCA based?) should be directed
to the model's authors, located in the AUTHORS file.

