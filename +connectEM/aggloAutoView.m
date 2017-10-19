function aggloAutoView(aggloFile,show,skelOutput)
% 
if ~exist('show','var') || isempty(show)
    show = 'wc';
end
if ~exist('skelOutput','var') || isempty(skelOutput)
    skelOutput = 1;
end

if ispc()
    load('G:/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
    p.saveFolder = strcat('G:',p.saveFolder);
    p.seg.root = strcat('G:',p.seg.root);
    aggloFolder = fullfile(p.saveFolder, 'aggloState');
    outputFolder = 'G:/tmpscratch/mbeining/autoView/';
else
    load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
    aggloFolder = fullfile(p.saveFolder, 'aggloState');
    outputFolder = '/tmpscratch/mbeining/autoView/';
end

if ~isempty(strfind(aggloFile,'dendrites'))
    agglos = load(fullfile(aggloFolder,[aggloFile,'.mat']),'dendrites');
    agglos = agglos.dendrites;
else
    agglos = load(fullfile(aggloFolder,[aggloFile,'.mat']),'axons');
    agglos = agglos.axons;
end

offset = [ 128, 128, 128 ].* [ 11.239999771118164, 11.239999771118164, 28 ];
bbox = [5504,8448,3328] .* [ 11.239999771118164, 11.239999771118164, 28 ];
smallbbox = bbox - 4000;
bbox = bsxfun(@plus,[zeros(3,1), bbox'],offset');
smallbbox = bsxfun(@plus,[zeros(3,1), smallbbox'],offset'+2000);


% % code for random dendrite ids created with the April version
% load(fullfile(aggloFolder,[dendritesApril.mat']),'dendrites','indBigDends');
% numSeg = arrayfun(@(x) size(x.nodes,1),bigD);
% bigD = dendrites(indBigDends & numSeg>50);
% s = RandStream('mt19937ar','Seed',0);
% bigD = bigD(randperm(s,numel(bigD),50));
% segIdsForRandDend = arrayfun(@(x) x.nodes(1,4),bigD);
% save(fullfile(aggloFolder,'randDendIds50.mat'),'segIdsForRandDend')
load(fullfile(aggloFolder,'randDendIds50.mat'),'segIdsForRandDend')

% load somata including merged somata
somaAgglos = connectEM.getCenterSomaAgglos;
somaSegIds = cell2mat(somaAgglos);
% remove duplicate segIds
[~,ic] = unique(somaSegIds);
duplicates = somaSegIds(setdiff(1:numel(somaSegIds),ic));
% somaDuplicateIds = cellfun(@(x) any(intersect(x,duplicates)),somaAgglos);
somaAgglos = cellfun(@(x) setdiff(x,duplicates),somaAgglos,'uni',0);
somaSegIds = cell2mat(somaAgglos);

somaLUT(somaSegIds) = repelem(1:numel(somaAgglos),cellfun(@numel,somaAgglos));

aggloSegIds = cell2mat(arrayfun(@(x) x.nodes(:,4),agglos,'uni',0));
aggloLUT = zeros(1,max(aggloSegIds));
aggloLUT(aggloSegIds)  = repelem(1:numel(agglos),arrayfun(@(x) numel(x.nodes(:,4)),agglos));


switch show
    case {'cells','wc'}
        
        % check which segId of soma is found in which dendrite agglo
        [ismem,ind] = ismember(somaSegIds,aggloSegIds);
        
        % create list with soma id in first row and dendrite id in second row
        dendOfSoma = unique([somaLUT(somaSegIds(ismem))',aggloLUT(aggloSegIds(ind(ismem)))'],'rows');
        
        % get each dend id which contains most of the seg ids of each soma
        aggloSomaId = accumarray(somaLUT(somaSegIds(ismem))',aggloLUT(aggloSegIds(ind(ismem)))',[],@mode);
        
        
        % check if dendrite overlaps with more than one soma!
        [ismem,ind] = ismember(aggloSegIds,somaSegIds);
        dendOfSoma2 = unique([somaLUT(somaSegIds(ind(ismem)))',aggloLUT(aggloSegIds(ismem))'],'rows');
        
        overlapMatrix = false(numel(somaAgglos),numel(agglos));
        ind = sub2ind(size(overlapMatrix),cat(1,dendOfSoma(:,1),dendOfSoma2(:,1)),cat(1,dendOfSoma(:,2),dendOfSoma2(:,2)));
        overlapMatrix(ind) = true;
%         in2Soma = find(sum(overlapMatrix,1)>1);
        any(sum(overlapMatrix,2))
        fprintf('%d agglos have overlap with multiple soma or vice versa\n',size(unique(cat(1,dendOfSoma,dendOfSoma2),'rows')) )
        
        % loop through somata, get all agglos belonging to each soma and write
        % them as nml
        raster = [4,3];
        fig = figure('Units','centimeters','Position',[0 -2 23 30]);
        startN = 1;
        if exist(fullfile(aggloFolder,[aggloFile,'_endings.mat']),'file')
            aggloEndings = load(fullfile(aggloFolder,[aggloFile,'_endings.mat']));
        else
            aggloEndings = [];
        end
        % delete(sprintf('SomaAgglo_%s',aggloFile))
        for n = 1:numel(somaAgglos)
            if skelOutput
                if ~isempty(aggloEndings)
                    idxComments = ismember(agglos(aggloSomaId(n)).nodes(:,4),aggloEndings{aggloSomaId(n)});
                    if any(idxComments)
                        agglos(aggloSomaId(n)).comments = repmat({''},size(agglos(aggloSomaId(n)).nodes,1),1);
                        agglos(aggloSomaId(n)).comments(idxComments) = repmat({'ending'},sum(idxComments),1);
                    end
                end
                skel = connectEM.generateSkeletonFromAggloNew(agglos(aggloSomaId(n)), {sprintf('SomaAgglo_%02d',n)} , outputFolder, max(aggloSegIds),[],sprintf('SomaAgglo_%s_%02d.nml',aggloFile,n));
            else
                skel = Superagglos.toSkel(agglos(aggloSomaId(n)));
            end
            %     theseDends = agglos(dendOfSoma(dendOfSoma(:,1)==n,2));
            %     skel = connectEM.generateSkeletonFromAggloNew(theseDends, arrayfun(@(x) sprintf('SomaAgglo_%02d_%03d',n,x),1:numel(theseDends),'uni',0) , outputFolder, max(dendriteSegIds),[],sprintf('SomaAggloAll_%s_%02d.nml',aggloFile,n));
            if ~ishandle(fig)
                fig = figure('Units','centimeters','Position',[0 -2 23 30]);
                startN = n;
            end
            subplot(raster(1),raster(2),rem(n-1,prod(raster))+1)
            skel.plot;
            Visualization.plotBbox(bbox);
            axis equal
            xlim([0 70000])
            ylim([0 100000])
            zlim([0 100000])
            if rem(n,prod(raster))==0 || n == numel(somaAgglos)
                pause(0.5)
                drawnow
                %         export_fig(sprintf('SomaAgglo_%s',aggloFile),'-pdf','-append')
                tprint(sprintf('SomaAgglo_%s_%02d-%02d',aggloFile,startN,n),'-pdf-R')
                pause(0.5)
                close(fig);
            end
        end
        
    case 'singles'
        
        
        %% find and plot the 50 pseudo random agglos
        [~,ind] = ismember(segIdsForRandDend,aggloSegIds);
        randDends = NaN(numel(segIdsForRandDend),1);
        randDends(ind~=0) = aggloLUT(aggloSegIds(ind(ind~=0)));
        raster = [4,3];
        fig = figure('Units','centimeters','Position',[0 -2 23 30]);
        startN = 1;
        for n = 1:numel(randDends)
            if ~ isnan(randDends(n))
                theseDends = agglos(randDends(n));
                if skelOutput
                    skel = connectEM.generateSkeletonFromAggloNew(theseDends,{sprintf('RandDend_%s_%02d_aggloId%d',aggloFile,n,randDends(n))} , outputFolder, max(aggloSegIds),[],sprintf('RandDend_%02d.nml',n));
                else
                    skel = Superagglos.toSkel(theseDends);
                end
                if ~ishandle(fig)
                    fig = figure('Units','centimeters','Position',[0 -2 23 30]);
                    startN = n;
                end
                subplot(raster(1),raster(2),rem(n-1,prod(raster))+1)
                
                Visualization.plotBbox(bbox,[0 1 0]);
                [~,inBBox] = skel.intersectsBbox(smallbbox,[],1);
                skel.plotValue(~inBBox);
                axis equal
                view(45,45)
                xlim([0 70000])
                ylim([0 100000])
                zlim([0 100000])
            end
            if (rem(n,prod(raster))==0 || n == numel(randDends))
                drawnow
                pause(0.5)
                print(sprintf('RandDend_%s_%02d-%02d',aggloFile,startN,n),'-dpdf')
                pause(1)
                close(fig);
            end
        end
end