% write out the agglos/components to nml
function write( agglosToWrite, sM, path, distthr )

%vertcat(agglos(aggloInd), somaAgglos(78), cut);
if nargin < 4
    distthr = 3000; % maximum distance between nodes that can be connected
    fprintf('distance threshold set to %d \n', distthr);
end

tic;
connectEM.generateSkeletonFromNodes(path,...
    cellfun(@(x) sM.point(:,x)', agglosToWrite,'uni',0), ...
    arrayfun(@(x) strcat('skelNod_',num2str(x,'%.2i')), 1:numel(agglosToWrite),'uni',0),[],[],[],distthr); 
fprintf('created %d skeletons from filtered agglos.\n', length(agglosToWrite)); toc;

end

