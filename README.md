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

Connect to GABA over the intranet with
```
ssh gaba.opt.rzg.mpg.de
```

Replace YOUR_DIRECTORY with a directory that you can write to.
You will be prompted for a username & password here, use the RZG one.
```
module load git/2.1.1
module load matlab/R2014b
git clone https://gitlab.mpcdf.mpg.de/mberning/pipeline.git YOUR_DIRECTORY
cd YOUR_DIRECTORY
```

You will be prompted for a username & password here. Use your GitHub account.
```
git submodule init
git submodule update
```

Now edit configuration.m to your needs with e.g.:`nano exampleUsage.m`. Save the file
by pressing Control+X and then 'Y' for "yes".

Finally, start matlab with `matlab -nosplash -nodesktop` and:
1. Set configuration for pipeline you just edited by running: `run configuration.m;`
2. Look at a preview of the segmentation using (second argument is region to be put into movie in voxel coordiantes): 
```
makeSegmentationPreviewMovie(p, [1001 1720; 1001 2280; 1001 1100])
```
This will print a file name of a segmentation movie you can use to judge quality (e.g. over vs. undersegmentation).
3. If you are not satisfied with the results, edit configuration.m again and repeat steps 1 & 2
4. Once you are satisfied with the results run `runPipeline(p)`

