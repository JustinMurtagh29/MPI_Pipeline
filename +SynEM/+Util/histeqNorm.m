function raw = histeqNorm( raw, raw_ref )
%HISTEQNORM Raw data normalization using histogram equalization.
% INPUT raw: uint8
%           The data that is normalized.
%       raw_ref: uint8
%           The reference raw data to which the histogram of raw is
%           adapted.
% OUTPUT raw: uint8
%           The transformed raw data.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

N = histcounts(raw_ref(:), -0.5:255.5);
raw = reshape(histeq(raw(:), N), size(raw));

end

