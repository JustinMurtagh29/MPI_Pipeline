function cells = getSomaNodes( p, graph, meta, rp, gb, probThreshold, sizeThreshold, somaIDs )
%GETSOMANODESPAR getSomaNodes in parallel
% INPUT     p, graph, meta, rp, gb: connectomics stuff 
%               default data that is used.
%           probThreshold: int
%               Only edges with prob > probThreshold will be taken.
%           sizeThreshold: int
%               Merge probability for the corresponding edges.
%           somaIDs: Nx1
%               IDs of the soma (corresponding to Christians nuclei)
% OUTPUT    
%           cells: cell containing:
%           nodes: [Nx3] int
%               meta.point of ids that are in the agglo.
%           name: str
%               The edges that connect the ids.
%           segIds: [Nx1]
%               seg ids that are in the agglo.
%           The edges that connect the ids.
% Author: Robin Hesse

    %%start Job 
    inputCell = cell(size(somaIDs,2),1);
    for i=1:size(inputCell,1)
       inputCell{i} = somaIDs(i);
    end
    params = {p, graph, meta, rp, gb, probThreshold, sizeThreshold};
    job = Cluster.startJob(@Soma.getSomaNodes, inputCell,...
        'sharedInputs',params,'name','somaAgglo','numOutputs',3,'cluster',...
        sprintf('-p -200 -l s_rt=99:59:30 -l h_rt=100:00:00 -l h_vmem=42G'));
    Cluster.waitForJob(job)
    cells = Cluster.fetchTaskOutputs(job,1:numel(inputCell));
    
end

