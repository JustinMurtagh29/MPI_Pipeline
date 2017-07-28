function [ rp, areaT ] = areaDepRP( group, scores, area )
%AREADEPRP Calculate RP for large synapses only.
% INPUT group: [Nx1] int
%           Target group vector.
%       scores: [Nx1] float
%           Prediction scores.
%       area: [Nx1] float
%           Interface area.
% OUTPUT rp: [Nx2] float
%           Optimal RP values considering all interfaces above an area
%               thresold.
%        areaT: [Nx1] float
%           The area threshold for the respective row in rp.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

areaT = 150:50:max(area);
rp = zeros(length(areaT), 2);

for i = 1:length(areaT)
    curGroup = group(area >= areaT(i));
    if ~any(curGroup)
        break;
    end
    curRP = SynEM.Eval.interfaceRP(curGroup, scores(area >= areaT(i)));
    [~,idx] = min(pdist2(curRP,[1, 1]));
    rp(i,:) = curRP(idx,:);
end

%delete large area thresholds
rp(i:end,:) = [];
areaT(i:end) = [];

end

