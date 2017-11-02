load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');

for idx = 1 : 11*17*13
    idx
    results{idx} = load(['/tmpscratch/kboerg/borders2/borders_' num2str(idx)],'findings');
    [X,Y,Z] = ind2sub([11,17,13],idx);
    fR = reshape(results{idx}.findings,ceil((diff(p.local(idx).bboxSmall,[],2)'+1)./[32,32,16]));
    globalD( ...
    (X-1)*16+1:(X-1)*16+size(fR,1), ...
    (Y-1)*16+1:(Y-1)*16+size(fR,2), ...
    (Z-1)*16+1:(Z-1)*16+size(fR,3)) = fR;
end

axonid = -28;
load('/gaba/u/mberning/results/pipeline/20170217_ROI/aggloState/dendrites_09.mat');
types = load('/gaba/scratch/sahilloo/forKevin/idxPostIn09.mat');

for i_x = 1 : size(globalD,1)
    i_x
    for i_y = 1 : size(globalD,2)
        for i_z = 1 : size(globalD,3)
            globalAxon(i_x,i_y,i_z) = 0;
            if ~isempty(globalD{i_x,i_y,i_z})
                globalAxon(i_x,i_y,i_z)=sum(globalD{i_x,i_y,i_z}(any(ismember(globalD{i_x,i_y,i_z}(:,1:2),axonid),2),3));
                targetAll(i_x,i_y,i_z)=sum(sum(ismember(globalD{i_x,i_y,i_z}(:,1:2),fIBD).*repmat(globalD{i_x,i_y,i_z}(:,3),1,2)));
                targetSmooth(i_x,i_y,i_z)=sum(sum(ismember(globalD{i_x,i_y,i_z}(:,1:2),double(types.idxAD)).*repmat(globalD{i_x,i_y,i_z}(:,3),1,2)));
                targetAD(i_x,i_y,i_z)=sum(sum(ismember(globalD{i_x,i_y,i_z}(:,1:2),double(types.idxSmooth)).*repmat(globalD{i_x,i_y,i_z}(:,3),1,2)));
            end
        end
    end
end

distanceT = bwdist(globalAxon>0);
for i_dist = 0 : 60
    selector = @(x)sum(sum(sum(x(distanceT>=i_dist&distanceT<i_dist+1))));
    graph_AD(i_dist+1) = selector(targetAD)/selector(targetAll);
    graph_Smooth(i_dist+1) = selector(targetSmooth)/selector(targetAll);
end    
    