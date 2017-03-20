function [ segIdsAgglomerate,edgesAgglomerate ] = getAgglomerate( p,meta,mode,islocal )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% get seg ids of agglomeration and transform meta data to agglomeration meta data
if nargin < 3
    mode = 'automatic';
end
    
if nargin < 4
    islocal = 1;
end
switch mode
    case 'manual'
        [fname,pname] = uigetfile('.nml','Please give merger mode tracing of %s.');
        skel = skeleton(fullfile(pname,fname));
        if islocal && isempty(strfind(p.seg.root,'G:'))
            p.seg.root = strcat('G:',p.seg.root);
        end
        numNodes = cellfun(@(x) size(x,1),skel.nodes);
        allNodes = cat(1,skel.nodes{:});
        segIdsAgglomerate = Seg.Global.getSegIds(p,allNodes(:,1:3));
        segIdsAgglomerate = mat2cell(segIdsAgglomerate,numNodes);
        edgesAgglomerate = skel.edges;
%         remainingEdges = cat(1,skel.edges{:});
    case 'automatic'
        if ~exist(fullfile(p.saveFolder,'agglomeration','initialAgglo.mat'),'file')
            alignEM.agglomerate(p,0)
        end
        segIdsAgglomerate = load(fullfile(p.saveFolder,'agglomeration','initialAgglo.mat'),'initialPartition','remainingEdges');
%         remainingEdges = segIdsAgglomerate.remainingEdges{1};
%         remainingEdges = cellfun(@(x) intersect(x,segIdsAgglomerate.remainingEdges{1}),segIdsAgglomerate.initialPartition{1},'uni',0);
        edgesAgglomerate = cell(numel(segIdsAgglomerate.initialPartition{1}),1);
        for n = 1:numel(segIdsAgglomerate.initialPartition{1})
            inthisbox = intersect(meta.segIds,segIdsAgglomerate.initialPartition{1}{n}); % get only segIds of agglomerate that are exist in the meta (i.e. are in the bbox)
            if ~isempty(inthisbox)
                im = ismember(segIdsAgglomerate.remainingEdges{1},inthisbox);
                edgesAgglomerate{n} = segIdsAgglomerate.remainingEdges{1}(all(im,2),:);
                IdsNotConsidered = setdiff(segIdsAgglomerate.remainingEdges{1}(im),edgesAgglomerate{n});  % find Ids that are in the box but lost their edge
                if ~isempty(IdsNotConsidered)
                    if isempty(edgesAgglomerate{n})
                        IdsNotConsidered = setdiff(IdsNotConsidered,inthisbox(1));
                        if ~isempty(IdsNotConsidered)
                            edgesAgglomerate{n} = cat(1,edgesAgglomerate{n},cat(2,ones(numel(IdsNotConsidered),1,'uint32')*inthisbox(1),IdsNotConsidered)); % connect these to first entry of agglomerate
                        end
                    else
                        edgesAgglomerate{n} = cat(1,edgesAgglomerate{n},cat(2,ones(numel(IdsNotConsidered),1,'uint32')*edgesAgglomerate{n}(1),IdsNotConsidered)); % connect these to first entry of agglomerate
                    end
                end
                segIdsAgglomerate.remainingEdges{1}(all(im,2),:) = []; % delete edges to speed up intersect of next iteration
                
            end
            segIdsAgglomerate.initialPartition{1}{n} = inthisbox;
        end
        segIdsAgglomerate = segIdsAgglomerate.initialPartition{1};
        edgesAgglomerate = edgesAgglomerate(~cellfun(@isempty,segIdsAgglomerate));
        segIdsAgglomerate = segIdsAgglomerate(~cellfun(@isempty,segIdsAgglomerate));
end
% notused = cellfun(@(x) setdiff(unique(x),meta.segIds),segIdsAgglomerate,'uni',0); % delete edges of segIDs that are not in the box
% [~,i1] = intersect(edgesAgglomerate(:,1),cat(1,notused{:}));
% [~,i2] = intersect(edgesAgglomerate(:,2),cat(1,notused{:}));
% edgesAgglomerate(unique([i1;i2]),:) = [];
% segIdsAgglomerate = cellfun(@(x) intersect(meta.segIds,unique(x)),segIdsAgglomerate,'uni',0); % get only unique segIds per agglomerate and only those that are in bbox
% edgesAgglomerate = edgesAgglomerate(~cellfun(@isempty,segIdsAgglomerate));
% segIdsAgglomerate = segIdsAgglomerate(~cellfun(@isempty,segIdsAgglomerate)); % delete segs not in bbox

end

