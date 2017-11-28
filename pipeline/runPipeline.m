function p = runPipeline(p, startStep, endStep,runlocal)
%RUNPIPELINE
% INPUT p: struct
%           Segmentation parameter struct.
%           (see also setParameterSettings)
%       startStep: (Optional) PipelineStep object
%           The step where the pipeline is started.
%           (Default: PipelineStep.Classification)
%       endStep: (Optional) PipelineStep object
%           The step where the pipeline is stopped.
%           (Default: PipelineStep.GraphConstruction)
%       runlocal: (Optional) boolean to switch if
%           some main functions of the pipeline should
%           be run on the interactive node or on a worker (default)

    % store run info in hidden folder of pipeline
    runInfo = Util.runInfo(); %#ok<NASGU>
    runInfoFolder = fullfile(p.saveFolder, '.runInfo');
    if ~exist(runInfoFolder, 'dir')
        mkdir(runInfoFolder);
    end
    outFile = fullfile(runInfoFolder, datestr(now, 30));
    Util.log('Storing run info at %s', outFile);
    save(outFile, 'runInfo');

    %check whether pipelineStep is specified else start with classification
    if ~exist('startStep','var') || isempty(startStep)
        startStep = PipelineStep.Classification;
    end
    if ~exist('endStep','var') || isempty(endStep)
        endStep = PipelineStep.GlobalGraphStruct;
    end
    if ~exist('runlocal','var') || isempty(runlocal)
        runlocal = 0;
    end

    % Runs CNN based forward pass for region defined in p.bbox,p.raw and saves as p.class
    % Because CNN is translation invariant, saved as KNOSSOS hierachy again
    % Pass bounding box as well as tileSize will be added in bigFwdPass otherwise (legacy stuff)
    % This uses CNN subfolder in code repository
    if startStep <= PipelineStep.Classification && ...
       endStep >= PipelineStep.Classification
        job = bigFwdPass(p, p.bbox);
        Cluster.waitForJob(job);
    end

    % If you set p.myelin.isUsed = true in configuration.m, the result of a myelin detection is
    % used to ensure a voxel detected as myelin will not be in one segment with voxel not detected
    % as myelin
    if isfield(p,'myelin') && p.myelin.isUsed && startStep <= PipelineStep.MyelinFix && ...
        endStep >= PipelineStep.MyelinFix

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

    if startStep <= PipelineStep.Segmentation && ...
       endStep >= PipelineStep.Segmentation
        job = miniSegmentation(p);
        Cluster.waitForJob(job);
    end

    % Transfer segmentation from p.local(X).tempSegFile to p.local(X).segFile and drop overlaps
    % Also in correspondences subfolder
    % Added routine to renumber CC of segments after cutting to non-overlapping region
    if startStep <= PipelineStep.OverlapRemoval && ...
       endStep >= PipelineStep.OverlapRemoval
        removeOverlaps(p);
    end

    % Make segmentation IDs unique over dataset (were only unique in each
    % tile before), this will be called global IDs from time to time
    % will work on segmentation see globalization subfolder
    if startStep <= PipelineStep.GlobalSegmentID && ...
       endStep >= PipelineStep.GlobalSegmentID
        job = globalizeSegmentation(p);
        Cluster.waitForJob(job);
    end

    % Find correspondences between tiles processed (segmented) in parallel
    % Uses correspondences subfolder in code repository
    if startStep <= PipelineStep.Correspondence && ...
       endStep >= PipelineStep.Correspondence
        job = correspondenceFinderGlobal(p);
        Cluster.waitForJob(job);
    end

    % Build segment meta data
    if startStep <= PipelineStep.BuildSegmentMetaData && ...
       endStep >= PipelineStep.BuildSegmentMetaData
        job = buildSegmentMetaData(p);
        Cluster.waitForJob(job);
    end

    if startStep <= PipelineStep.SegmentationPyramid && ...
       endStep >= PipelineStep.SegmentationPyramid
        % Create resolution pyramid for the segmentation
        tic; fprintf('Downsampling segmentation... ');
        thisBBox = [1, 1, 1; (ceil(p.bbox(:, 2) ./ 1024) .* 1024)']';
        createResolutionPyramid(p.seg, thisBBox, [], true);
        fprintf('done!\n'); toc;
    end

    if startStep <= PipelineStep.CompressSegmentation && ...
       endStep >= PipelineStep.CompressSegmentation
        % Compress segmentation (at all resolutions)
        compressSegmentation(p.seg);
    end

    % Construct graph on globalized version of segmentation
    % This will create p.local(:).edgeFile & borderFile
    % See graphConstruction subfolder
    if startStep <= PipelineStep.GraphConstruction && ...
       endStep >= PipelineStep.GraphConstruction
        job = graphConstruction(p);
        Cluster.waitForJob(job);
    end

    % Run synapse detection as presented in Staffler et. al, 2017
    % see: http://biorxiv.org/content/early/2017/01/22/099994
    if startStep <= PipelineStep.SynapseDetection && ...
       endStep >= PipelineStep.SynapseDetection
        classifierPath = fullfile(fileparts(fileparts( ...
            mfilename('fullpath'))), ...
            '+SynEM', 'data', 'SynEMPaperClassifier.mat');
        [p, job] = SynEM.Seg.pipelineRun(p, classifierPath);
        % Save parameter file to new
        Util.save([p.saveFolder 'allParameterWithSynapses.mat'], p);
        disp('SynEM modified parameter file. Load allParameterWithSynapses.mat for the next steps.')
        % Wait for completion of job
        Cluster.waitForJob(job);
    end

    % Calculate raw features on smaller (wrt SynEM) borders (down to 10 voxel)
    if startStep <= PipelineStep.RawFeatures && ...
       endStep >= PipelineStep.RawFeatures
        job = connectEM.calculateFeatures(p, 'Raw');
        Cluster.waitForJob(job);
    end

    % Calculate features on CNN output as well
    if startStep <= PipelineStep.ClassFeatures && ...
       endStep >= PipelineStep.ClassFeatures
        job = connectEM.calculateFeatures(p, 'Class');
        Cluster.waitForJob(job);
    end

    % Calculate neurite continuity predictions
    if startStep <= PipelineStep.NeuriteContinuityPrediction && ...
       endStep >= PipelineStep.NeuriteContinuityPrediction
        job = connectEM.predictDataset(p);
        Cluster.waitForJob(job);
    end

    %Save the global SVG data
    if startStep <= PipelineStep.SaveGlobalSvgData && ...
            endStep >= PipelineStep.SaveGlobalSvgData
        if runlocal
            Seg.Global.saveGlobalSvgData(p,[],[],1);
        else
            job = collectSvgDataOnCluster(p);
            Cluster.waitForJob(job);
        end
    end

    %Create graph struct
    if startStep <= PipelineStep.GlobalGraphStruct && ...
            endStep >= PipelineStep.GlobalGraphStruct
        if runlocal
            collectGlobalGraphStruct(p);
        else
            job = collectGraphStructOnCluster(p);
            Cluster.waitForJob(job);
        end
    end
end

% Some comments that one might want to run in addition:
% if pipelineStep <= PipelineStep.HeuristicLookup
%  connectEM.getHeuristicsResult(p)         % -> lookups result from heuritic (nuclei, vessel, myelin) detections
% end
% if pipelineStep <= PipelineStep.Agglomeration
%  connectEM.agglomerate                   % -> Generate CC of graph at a threshold chosen in script
% end
