% Code to generate boutons from the synapses.mat using a given threshold

function boutons = boutonGeneration(p)

cubeIndices = 1:numel(p.local);
synScoreThr = 0;
boutonIDs = [];
boutonsCoMs = [];
for i =1:length(cubeIndices)

	m=load(p.local(i).synapseFile);
	scores = single(m.scores);
	m=load(p.local(i).edgeFile)
	edges = m.edges;
	boutonIDs = vertcat(boutonIDs,edges(scores>0,:));
end

%load com.mat to get CoMs of boutonIDs
load([p.saveFolder 'globalCoMList.mat']);

boutonCoMs = comList(boutonIDs,:);

boutons.boutonIDs = boutonIDs;
boutons.boutonCoMs = boutonCoMs;
save([p.saveFolder 'Boutons.mat'],'boutons');

end
