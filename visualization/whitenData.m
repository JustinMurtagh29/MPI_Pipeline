function white = whitenData(black, epsilon)

mu = mean(black); 
black = bsxfun(@minus, black, mu);
A = black'*black;
[V,D,~] = svd(A);
M = sqrt(size(black,1)-1)*V*sqrtm(inv(D + eye(size(D))*epsilon))*V';
white = black*M;

end

