%% This is an old version comparing Cortex and Retina Intracellular stain
a = readKnossosRoi('Z:\CortexConnectomics\shared\2012-09-28_ex145_07x2\mag1\', '2012-09-28_ex145_07x2_mag1', [4500 5500; 3000 4000; 1200 2200]);
b = readKnossosRoi('\\P1-374\Musik\e_k0563\k0563_mag1', '100527_k0563_mag1', [1800 2800; 2300 3300; 2300 3300]);

c = a(250:750, 250:750, 250:750);
d = b(250:750, 250:750, 250:750);

e = readKnossosRoi('Z:\CortexConnectomics\shared\2012-09-28_ex145_07x2\mag1\', '2012-09-28_ex145_07x2_mag1', [1000 1500; 1000 1500; 1000 1500]);

figure;
subplot(3,1,1);
hist(single(c(:)), 60);
xlim([0 255]);
title('07x2 centered');
set(gca, 'YScale', 'log');
subplot(3,1,2);
hist(single(d(:)), 60);
xlim([0 255]);
title('e_k0563 centered');
set(gca, 'YScale', 'log');
subplot(3,1,3);
hist(single(e(:)), 60);
xlim([0 255]);
title('07x2 upper left');
set(gca, 'YScale', 'log');
