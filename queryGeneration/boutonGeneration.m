% Code to generate boutons from the synapses.mat using a given threshold

function boutons = boutonGeneration(p)

cubeIndices = 1:numel(p.local);
synScoreThr = 0;
boutonIDs = [];
boutonsCoMs = [];
for i =1:length(cubeIndices)

	m=load([p.local(i).saveFolder 'synapses.mat']); % change to p.local(i).synapseFile
	scores = single(m.scores);
        m = load([p.local(i).probFile]);
        prob = m.prob;
        scores(prob>0.5,:)=-Inf;
	m=load(p.local(i).edgeFile);
	edges = m.edges;
	boutonIDs = vertcat(boutonIDs,edges(scores>0));
end
boutonIDs = unique(boutonIDs);
%load com.mat to get CoMs of boutonIDs
load([p.saveFolder 'globalCoMList.mat']);

boutonCoMs = globalCoMList(boutonIDs,:);

boutons.boutonIDs = boutonIDs;
boutons.boutonCoMs = boutonCoMs;
save([p.saveFolder 'Boutons.mat'],'boutons');

end
