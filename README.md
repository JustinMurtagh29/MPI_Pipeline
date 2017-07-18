Introduction to the Pipline Repository
======================================


Preparation
-----------
1. Get a account at the Max Planck Computing and Data Facility (MPCDF). Make sure you can connect to the GABA cluster at the Rechenzentrum Garching (see [here](https://wiki.hest.brain.mpg.de/doku.php?id=knowledge:organization:it:connecting_to_garching) for a guide)
2. Send an email to heiko.wissler@brain.mpg.de with your GitLab username to be added to the respective teams/repos
3. Wait for a positive reply from Heiko


Connecting to GABA
------------------

Connect to GABA over the intranet with
```
ssh gaba.opt.rzg.mpg.de
```


Setting up the code
-------------------

Replace YOUR_DIRECTORY with a directory that you can write to.
You will be prompted for a username & password here, use the RZG one.
```
module load git
module load matlab/R2015b
git clone git@gitlab.mpcdf.mpg.de:connectomics/pipeline.git YOUR_DIRECTORY
cd YOUR_DIRECTORY
```

Download all the dependencies with the following command. You will be prompted for a username & password here. Use your GitLab account.
```
git submodule update --init --recursive
```


Configuration
-------------

Next, you need to create the `configuration.m` file:
```
cp configuration-sample.m configuration.m
```

The `configuration.m` file is special in that it is ignored by git and excluded from versioning. This way you can be sure that your `configuration.m` will never be overwritten by a `git pull`.

Now edit configuration.m to your needs with `nano configuration.m`. Save the file
by pressing Control+X and then 'Y' for "yes". Finally, press enter.


Running the pipeline
--------------------

Finally, start MATLAB with `matlab -nosplash -nodesktop`. Please make sure you are **in the pipline directory** when starting MATLAB. Then it's time for the kick-off:

1. Run the startup script
```
run startup.m;
```
2. Set configuration for pipeline you just edited by running:
```
run configuration.m;
```
3. Look at a preview of the segmentation using: 
```
makeSegmentationPreviewMovie(p, [1001, 1001, 1001, 720, 1280, 100])
```
The second input argument is the bounding box in webKNOSSOS format, i.e., `[x_min, y_min, z_min, x_width, y_width, z_width]`.   
Alternatively, you can specify the bounding box in MATLAB format, i.e., `[x_min, x_max; y_min, y_max; z_min, z_max]`.   
This will print a file name of a segmentation movie you can use to judge quality (e.g. over vs. undersegmentation).
4. If you are not satisfied with the results, edit configuration.m again and repeat steps 1 & 2
5. Once you are satisfied with the results run:
```
runPipeline(p)
```

Updating
--------

To update your pipeline repository to the latest version, just switch into the pipline directory and run these two lines of code:
```
git pull origin master
git submodule update --init --recursive
```
