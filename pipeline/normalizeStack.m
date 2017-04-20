function out = normalizeStack( in, meanValue, stdValue )
%out = normalizeStack( in, meanValue, stdValue ) Subtracts mean and devide by std

% This is necessary as auxiliary methods should always be backwards compatible
% This sets values as used for all pipelines on 07x2 prior to now
if nargin == 1
    meanValue = 122;
    stdValue = 22;
end

if ~isfloat(in)
    in = single(in);
end

out = in-meanValue;
out = out./stdValue; 

end

