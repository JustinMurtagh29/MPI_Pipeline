function [ raw, target, targetWeights ] = loadDataSTH( path )
%LOADDATASTH Calculate ST and Hessian as additional input channels.
% INPUT path: Path to training data file.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

m = load(path);
raw = (single(m.raw) - 122)./22;
target = single(m.target);
targetWeights = ~(m.target == 0);
sigma = 12./[11.24,11.24,28];
filterSiz = ceil(2*sigma);
H = Hessian( raw, sigma, filterSiz );
H = H{3};
H = H./2.2;
ST = StructureTensor(raw,sigma,filterSiz,sigma,filterSiz);
ST = ST{3};
raw = cat(4,raw,H,ST);

end
