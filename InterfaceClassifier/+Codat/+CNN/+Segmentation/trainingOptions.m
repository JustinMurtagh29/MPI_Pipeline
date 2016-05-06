%% Training options for local debugging.
options.display = 10;
options.tr_size = [10, 10, 10];
options.val_size = [100, 100, 100];
options.val_iter = 0;
options.val_fwd_alg = 'fft2';
options.gpuDev = false;
options.max_iter = 1000;
options.augment = false;
options.snapshot = 10000;
options.snapshot_name = 'myCnet';
options.save_imp = false;
options.class_ratio = [];
options.skip_single_class_cubes = false;
options.data_pre = func2str(@Codat.CNN.Segmentation.loadData);
options.lr_policy = 'step';
options.step_size = 100000;
options.gamma = 0.1;
options.rot_inv = false;

%% training stacks
load('E:\workspace\data\cortexTrainingData\cortexTrainingDataParameter.mat')
stacksT = {stacks(1:10).stackFile};
stacksV = {stacks(11).stackFile};