function y = softmax( x )
%SOFTMAX Softmax activation function.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

num = exp(x);
denom = sum(num,4);
y = bsxfun(@rdivide,num,denom);

end

