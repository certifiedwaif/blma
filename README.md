Bayesian Linear Model Averaging
===============================
 [![Travis-CI Build Status](https://travis-ci.org/certifiedwaif/blma.svg?branch=master)](https://travis-ci.org/certifiedwaif/blma) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/certifiedwaif/blma?branch=master&svg=true)](https://ci.appveyor.com/project/certifiedwaif/blma)
 Linux and Mac installation instructions
---------------------------------------
If you have gsl installed,
 ```r
library(devtools)
devtools::install_devtools("certifiedwaif/blma")
```
 should install and build the package. The package relies on the gsl library
being installed. On Windows, gsl will be downloaded as part of the package
installation. On Mac OS X, if you are using Homebrew, gsl can be installed by
executing the command
 ```shell
brew install gsl
```
 If you are using Linux, you should consult your distribution's documentation on
how to install gsl.
 Windows installation instructions
---------------------------------
 Install the blma package via the R commands:
 ```r
library(devtools)
install_github("certifiedwaif/blma")
```
 Special thanks to Jerone Ooms [@opencpu](http://twitter.com/opencpu) for his
help with getting the package building on Windows automatically, and to Maëlle
Salmon [@ma_salmon](http://twitter.com/ma_salmon) for CCing him in the first
place and introducing me to Travis-CI and RAppVeyor.
