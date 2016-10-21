#!/usr/bin/env bash

codedir=`pwd`
cidir=`mktemp -d`

# pull code for CI suite
cd $cidir
git clone git@gitlab.mpcdf.mpg.de:connectomics/matlab-ci.git .
matlab -nosplash -nodisplay -nojvm -r "matlabCI('/home/amotta/Desktop/config.mat', '$codedir'); exit;"

# reset
cd $codedir
rm -rf $cidir
