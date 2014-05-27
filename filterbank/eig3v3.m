function [varargout] = eig3v3(A)
% function D = eig3(A)
%
% Compute in one shot the eigen-values of multiples (3 x 3) matrices
%
% INPUT:
%   A: (n x 3 x 3) array
% OUTPUT:
%   D: (n x 3). EIG3 returns in D(k,:) three eigen-values of A(:,:,k)
%   
% See also: CardanRoots, eig2, eig
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% History:
%     Original 20-May-2010
%
% function [V,D] = EIG3(A)
%
% Compute in one shot the eigen-values of multiples (3 x 3) matrices
%
% INPUT:
%   A: (n x 3 x 3) array
% OUTPUT:
%   D: (n x 3). EIG3 returns in D(k,i) the ith eigen-value of A(:,:,k)
%   V: (n x 3 x 3). EIG3 returns in V(k,:,i)  a normalized ith eigen-vector of A(:,:,k)
%   
if size(A,2) ~= 3 || size(A,3) ~= 3
    error('A must be [nx3x3] array');
end


P3 = 1;
% Trace
P2 = -(A(:,1)+A(:,5)+A(:,9));
% Principal minors
M11 = A(:,5).*A(:,9) - A(:,8).*A(:,6);
M22 = A(:,9).*A(:,1) - A(:,3).*A(:,7);
M33 = A(:,1).*A(:,5) - A(:,4).*A(:,2);
P1 = (M11 + M22 + M33);
% Determinant
P0 = - A(:,1).*M11 ...
    + A(:,4).*(A(:,2).*A(:,9)-A(:,8).*A(:,3)) ...
    - A(:,7).*(A(:,2).*A(:,6)-A(:,5).*A(:,3));

D = CardanRoots(P3, P2, P1, P0);
D = sortrows(D.',1).';
    
    
if nargout<2
    varargout={D};
else
    I = permute(eye(3),[3,1,2]);
    C1 = mprod(A-bsxfun(@times,I,D(:,2)),A-bsxfun(@times,I,D(:,1)));
    C2 = mprod(A-bsxfun(@times,I,D(:,3)),A-bsxfun(@times,I,D(:,1)));
    
    n=size(A,1);
    V=zeros(size(A),class(A));
    
    idx=any(C1,2);
    idx(:,:,1)=idx(:,:,1)&~(idx(:,:,2)|idx(:,:,3));
    idx(:,:,2)=idx(:,:,2)&~idx(:,:,3);
undefined=~any(idx,3);
	idx(undefined,:,3)=1;
    idx=repmat(idx,[1,3,1]);
    V(:,:,3)=reshape(C1(idx),[n,3,1]);
	V(undefined,:,3)=randn(sum(undefined),3,class(A));
    
    idx=any(C2,2);
    idx(:,:,1)=idx(:,:,1)&~(idx(:,:,2)|idx(:,:,3));
    idx(:,:,2)=idx(:,:,2)&~idx(:,:,3);
undefined=~any(idx,3);
        idx(undefined,:,3)=1;

    idx=repmat(idx,[1,3,1]);    
    V(:,:,2)=reshape(C2(idx),[n,3,1]);
 V(undefined,:,2)=cross(V(undefined,:,3),randn(sum(undefined),3,class(A)));


    V(:,:,1)=cross(V(:,:,3).',V(:,:,2).').';
    V=bsxfun(@rdivide,V,sqrt(sum(real(V).^2+imag(V).^2,2)));
    D2=zeros(size(A),class(A));
    D2(:,1,1)=D(:,1);
    D2(:,2,2)=D(:,2);
    D2(:,3,3)=D(:,3);    
    varargout={V,D2};
end

end
