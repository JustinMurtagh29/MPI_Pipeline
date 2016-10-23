#!/usr/bin/env bash

echo "> Updating .matlab-ci file..."
echo "apiProjectId: $CI_PROJECT_ID" >> .matlab-ci
echo "apiPrivateToken: $PRIVATE_TOKEN" >> .matlab-ci

codedir=`pwd`
cidir=`mktemp -d`

# pull code for CI suite
cd $cidir
git clone git@gitlab.mpcdf.mpg.de:connectomics/matlab-ci.git .
matlab -nosplash -nodisplay -nojvm -r "matlabCI('$codedir'); exit;"

# reset
cd $codedir
rm -rf $cidir
