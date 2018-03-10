#!/usr/bin/env bash

echo "> Updating .matlab-ci file..."
echo "apiProjectId: $CI_PROJECT_ID" >> .matlab-ci
echo "apiPrivateToken: $PRIVATE_TOKEN" >> .matlab-ci

codedir=`pwd`
cidir=`mktemp -d`

# pull code for CI suite
cd $cidir
git clone git@gitlab.mpcdf.mpg.de:connectomics/matlab-ci.git .

# select MATLAB version and run
export PATH=$PATH:/usr/local/matlab/r2017b/bin
matlab -nosplash -nodisplay -nojvm -r "runMatlabCI('$codedir');"

# capture exit code
errcode=$?

# reset
cd $codedir
rm -rf $cidir

# repeat error
exit $errcode
