function runPipeline(p, pipelineStep)

    %check whether pipelineStep is specified else start with classification
    if ~exist('pipelineStep','var') || isempty(pipelineStep)
        pipelineStep = PipelineStep.Classification;
    end

    % Runs CNN based forward pass for region defined in p.bbox,p.raw and saves as p.class	
    % Because CNN is translation invariant, saved as KNOSSOS hierachy again 
    % Pass bounding box as well as tileSize will be added in bigFwdPass otherwise (legacy stuff)
    % This uses CNN subfolder in code repository    
    if pipelineStep <= PipelineStep.Classification
        job = bigFwdPass(p, p.bbox);
        Cluster.waitForJob(job);
    end
    
    % If you set p.myelin.isUsed = true in configuration.m, the result of a myelin detection is
    % used to ensure a voxel detected as myelin will not be in one segment with voxel not detected
    % as myelin
    if p.myelin.isUsed && pipelineStep <= PipelineStep.MyelinFix
    
        %define new classification prefix (NOTE: This will not be saved, but temporary files anyway)
        newPrefix = [p.class.prefix, '_mod'];
    
        % run myelin detection
        job = Myelin.runFix(p, newPrefix);
        Cluster.waitForJob(job);
    
        % use modified classification in subsequent steps
        p.class.prefix = newPrefix;
        clear newPrefix;
    end

    % Runs watershed based segmentation for region defined in p.bbox on p.raw and saves as p.local(X).tempSegFile
    % Because watershed segmentation has FOV effects (no translation invariance), processed with large
    % overlap and later joined together (see correspondences)
    % Uses segmentation subfolder in code repository
    if pipelineStep <= PipelineStep.Segmentation
        job = miniSegmentation(p);
        Cluster.waitForJob(job);
    end
    
    % Find correspondences between tiles processed (segmented) in paralell
    % Uses correspondences subfolder in code repository
    if pipelineStep <= PipelineStep.Correspondence
        job = correspondenceFinder(p);
        Cluster.waitForJob(job);
    end
    
    % Transfer segmentation from p.local(X).tempSegFile to p.local(X).segFile and drop overlaps 
    % Also in correspondences subfolder
    % Added routine to renumber CC of segments after cutting to non-overlapping region
    % This also involves renumbering all correspondences if a segment is renumbered in this step
    if pipelineStep <= PipelineStep.OverlapRemoval
        removeOverlaps(p);
    end
    
    % Make segmentation IDs unique over dataset (were only unique in each
    % tile before), this will be called global IDs from time to time
    % will work on correspondences and segmentation
    % see globalization subfolder 
    % globalization is run in two parts: globalizeSegmentation and globalizeCorrespondences
    if pipelineStep <= PipelineStep.GlobalSegmentID
        job = globalizeSegmentation(p);
        Cluster.waitForJob(job);
    end 
    
    % Build segment meta data
    if pipelineStep <= PipelineStep.BuildSegmentMetaData
        job = buildSegmentMetaData(p);
        Cluster.waitForJob(job);
    end
    
    if pipelineStep <= PipelineStep.GlobalCorrespondences
        job = globalizeCorrespondences(p);
        Cluster.waitForJob(job);
    end
        
    % Create resolution pyramid for the segmentation
    display('Downsampling segmentation:');
    tic;
    thisBBox = [1 1 1; (ceil(p.bbox(:,2)./1024).*1024)']';
    % This weird command line stuff is necessary to deference symbolic links
    [~, thisRoot] = system(['readlink -f ' p.seg.root ' < /dev/null']);
    thisRoot = [strrep(thisRoot, sprintf('\n'), '') filesep];
    createResolutionPyramid(thisRoot, p.seg.prefix, thisBBox, strrep(thisRoot, '/1/', ''), true);
    toc;

    % Construct graph on globalized version of segmentation
    % This will create p.local(:).edgeFile & borderFile
    % See graphConstruction subfolder
    if pipelineStep <= PipelineStep.GraphConstruction
        job = graphConstruction(p);
        Cluster.waitForJob(job);
    end
    
    % Run synapse detection as presented in Staffler et. al, 2017
    % see: http://biorxiv.org/content/early/2017/01/22/099994
    if pipelineStep <= PipelineStep.SynapseDetection
        [p, job] = SynEM.Seg.pipelineRun(p);
        % Save parameter file to new 
        Util.save([p.saveFolder 'allParameterWithSynapses.mat'], p);
        % Wait for completion of job
        Cluster.waitForJob(job);
    end
        
    % Calculate raw features on smaller (wrt SynEM) borders (down to 10 voxel)
    if pipelineStep <= PipelineStep.RawFeatures
        job = connectEM.calculateRawFeatures(p);
        Cluster.waitForJob(job);
    end
    
    % Calculate features on CNN output as well
    if pipelineStep <= PipelineStep.ClassFeatures
        job = connectEM.calculateClassFeatures(p);
        Cluster.waitForJob(job);
    end
    
    % Calculate neurite continuity predeictions
    if pipelineStep <= PipelineStep.NeuriteContinuityPrediction
        job = connectEM.predictDataset(p);
        Cluster.waitForJob(job);
    end
        
    %Save the global SVG data
    if pipelineStep <= PipelineStep.SaveGlobalSvgData
        job = collectSvgDataOnCluster(p);
        Cluster.waitForJob(job);
    end
    
    %Create graph struct 
    if pipelineStep <= PipelineStep.GlobalGraphStruct
        job = collectGraphStructOnCluster(p);
        Cluster.waitForJob(job);
    end
end

% Some comments that one might want to run in addition:
% if pipelineStep <= PipelineStep.HeuristicLookup
%   connectEM.getHeuristicResult(p)         % -> lookups result from heuritic (nuclei, vessel, myelin) detections
% end
% if pipelineStep <= PipelineStep.Agglomeration
%   connectEM.agglomerate                   % -> Generate CC of graph at a threshold chosen in script
% end