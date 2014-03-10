function b = detectBloodVessels()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
bbox = [1 10240; 1 7168; 1 3422];
root = '/zdata/manuel/data/cortex/2012-09-28_ex145_07x2/mag1/';
prefix = '2012-09-28_ex145_07x2_mag1';
xTunnel = 1800;

display('Reading data from KNOSSOS hierachy');
tic;
raw = readKnossosRoi(root, prefix, bbox);
toc;

display('Thresholding & Filling:');
tic;
% generate logical arrrays
a = raw > 160 | raw < 10;
b = zeros(size(a), 'uint8');
for i=1:size(a,3)
        a(:,:,i) = bwareaopen(a(:,:, i), 1000, 8);
        temp = padarray(a(:,:,i), [1 1], 1);
        idx = 1;
        while temp(xTunnel,idx) == 1
                temp(xTunnel,idx) = 0;
                idx = idx + 1;
        end
        temp = imfill(temp, 'holes');
        temp([1 end],:) = [];
        temp(:,[1 end]) = [];
        b(:,:,i) = temp;
        display(['Progress : ' num2str(i/size(a,3)*100, '%.2f') '%']);
end

end

