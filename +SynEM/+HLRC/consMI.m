function [mi, mi_ref] = consMI( raw, raw_ref )
%CONSMI Mutual information between consecutive images.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

mi= zeros(size(raw, 3) - 1, 1);
for i = 1:(size(raw, 3) - 1)
    mi(i) = SynEM.HLRC.MI_GG(raw(:,:,i), raw(:,:,i+1));
end

if exist('raw_ref', 'var') && ~isempty(raw_ref)
    mi_ref = zeros(size(raw_ref, 3) - 1, 1);
    for i = 1:(size(raw_ref, 3) - 1)
        mi_ref(i) = SynEM.HLRC.MI_GG(raw_ref(:,:,i), raw_ref(:,:,i+1));
    end
else
    mi_ref = [];
end

end

