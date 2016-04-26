function [ D ] = eig3S( a, b, c, d, e, f )
%EIG3S Calculate the eigenvalues of multipe 3x3 real symmetric matrices.
% INPUT All inputs should be [Nx1] vectors. The i-th entry of each input
%       corresponds to the entries of a real symmetric matrix of the form
%       A(i) = [a(i), b(i), c(i);
%               b(i), d(i), e(i);
%               c(i), e(i), f(i)]
% OUTPUT D: [Nx3] array containing the eigenvalues of A(i). Eigenvalues are
%           sorted in descending order, i.e.
%           D(:,1) >= D(:,2) >= D(:,3)
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

p1 = b.^2 + c.^2 + e.^2;
q = (a + d + f)./3; %trace of A
p2 = (a - q).^2 + (d - q).^2 + (f - q).^2 + 2*p1;
p = sqrt(p2/6);

%determinant det(1/p(B-q*Id))/2
r = 0.5./p.^3.*((a-q).*(d-q).*(f-q) + 2.*b.*c.*e -  ...
    e.^2.*(a-q) - b.^2.*(f-q) - c.^2.*(d-q));

phi = zeros(size(r),'like',a);
idx1 = r >= 1;
idx2 = r <= -1;
phi( idx2) = pi/3;
phi( ~(idx1 & idx2)) = real(acos(r(~(idx1 & idx2)))/3);

%the eigenvalues satisfy D(1,:) >= D(2,:) >= D(3,:)
D = zeros([length(a), 3], 'like', a);
D(:,1) = q + 2.*p.*cos(phi);
D(:,3) = q + 2.*p.*cos(phi + 2*pi/3);
D(:,2) = 3.*q - D(:,1) - D(:,3);

end

%pseudocode from wikipedia
%(reference: Smith, Oliver K. (April 1961), "Eigenvalues of a symmetric
% 3 × 3 matrix.", Communications of the ACM 4 (4): 168,
% doi:10.1145/355578.366316):
%
% % % Given a real symmetric 3x3 matrix A, compute the eigenvalues
% % 
% % p1 = A(1,2)^2 + A(1,3)^2 + A(2,3)^2
% % if (p1 == 0) 
% %    % A is diagonal.
% %    eig1 = A(1,1)
% %    eig2 = A(2,2)
% %    eig3 = A(3,3)
% % else
% %    q = trace(A)/3
% %    p2 = (A(1,1) - q)^2 + (A(2,2) - q)^2 + (A(3,3) - q)^2 + 2 * p1
% %    p = sqrt(p2 / 6)
% %    B = (1 / p) * (A - q * I)       % I is the identity matrix
% %    r = det(B) / 2
% % 
% %    % In exact arithmetic for a symmetric matrix  -1 <= r <= 1
% %    % but computation error can leave it slightly outside this range.
% %    if (r <= -1) 
% %       phi = pi / 3
% %    elseif (r >= 1)
% %       phi = 0
% %    else
% %       phi = acos(r) / 3
% %    end
% % 
% %    % the eigenvalues satisfy eig3 <= eig2 <= eig1
% %    eig1 = q + 2 * p * cos(phi)
% %    eig3 = q + 2 * p * cos(phi + (2*pi/3))
% %    eig2 = 3 * q - eig1 - eig3     % since trace(A) = eig1 + eig2 + eig3
% % end