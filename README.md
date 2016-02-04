Short introduction on how to use this repositorium on gaba cluster
==================================================================

Preparation:
============

1. Get a github account by signing up on the website github.com
2. Get a account at the RZG and make sure you can connect to the gaba cluster at the Rechenzentrum Garching (see [here](https://wiki.hest.brain.mpg.de/doku.php?id=knowledge:organization:it:connecting_to_garching) for a guide)
3. Send an email to heiko.wissler@brain.mpg.de with your github and RZG username to be added to the respective teams/repos
4. Wait for a positive reply from Heiko

Usage:
========================

Replace YOUR_DIRECTORY with a directory that you can write to
You will be prompted for a username & password here, use the RZG one
```
module load git/2.1.1
module load matlab/R2014b
git clone https://gitlab.mpcdf.mpg.de/mberning/pipeline.git YOUR_DIRECTORY
cd YOUR_DIRECTORY
```

You will be prompted for a username & password here, use github account
```
git submodule init
git submodule update
matlab -nosplash -nodesktop
```

Now open/run exampleUsage.m once you have edited to to your needs.
