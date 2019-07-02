%train data
trainFiles = {'C:\Users\loombas\Downloads\+TypeEM_+Data_tracings_ex145-07x2-roi2017_dense-v3_region-1.nml',...
              'C:\Users\loombas\Downloads\+TypeEM_+Data_tracings_ex145-07x2-roi2017_dense-v3_region-2.nml',...
              'C:\Users\loombas\Downloads\+TypeEM_+Data_tracings_ex145-07x2-roi2017_dense-v3_region-3.nml',...
              'C:\Users\loombas\Downloads\+TypeEM_+Data_tracings_ex145-07x2-roi2017_dense-v2_region-4.nml'};

skels = cellfun(@(x) skeleton(x), trainFiles);

numAxon = [];
numDend = [];
numGlia = [];
numAxonN = [];
numDendN = [];
numGliaN = [];

for i=1:numel(skels)
    [numA, numAN] = doThis(skels(i),'Axon');
    [numD, numDN] = doThis(skels(i),'Dendrite');
    [numG, numGN] = doThis(skels(i),'Astrocyte');
    
    numAxon = vertcat(numAxon,numA);
    numAxonN = vertcat(numAxonN,numAN);
    numDend = vertcat(numDend,numD);
    numDendN = vertcat(numDendN,numGN);
    numGlia = vertcat(numGlia,numG);
    numGliaN = vertcat(numGliaN,numGN);
    clear numA numAN numD numDN numG numGN
end


%%
figure
hold on
bar(cat(2,numGlia, numAxon, numDend))
bar(cat(2,numGliaN, numAxonN, numDendN),'FaceColor','none','EdgeColor','g')
legend({'+ glia','+ axon','+ dendrite','- glia','- axon','- dendrite' })
xlabel('Test box Id')
ylabel('Numder of labels')
set(gca,'xtick',[1,2,3,4])

function [numAxonTemp, numAxonNTemp] = doThis(skel,name)
    % find trees with name for class axon
    skelAxon = skel.keepTreeWithName( name, 'exact' );
    % find nodes per class
    numAxonTemp = size(skelAxon.getNodes,1);
    % remove merger comment nodes
    toDel = skelAxon.getNodesWithComment('Merger','','',true);
    if ~isempty(toDel)
        toDel = vertcat(toDel{:});
        numAxonTemp = numAxonTemp - size(toDel,1);
    end
    skelAxonN = skel.deleteTreeWithName(name, 'exact' );
    % find nodes per class
    numAxonNTemp = size(skelAxonN.getNodes,1);
    % remove merger comment nodes
    toDel = skelAxonN.getNodesWithComment('Merger','','',true);
    if ~isempty(toDel)
        toDel = vertcat(toDel{:}); 
        numAxonNTemp = numAxonNTemp - size(toDel,1);
    end    
end