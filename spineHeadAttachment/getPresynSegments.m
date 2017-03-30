function presynSegments = getPresynSegments(graph, documentID, windowspath)
if ~exist('windowspath', 'var')
    windowspath = false;
end
if windowspath
    addpath('D:\Git\st07x2\globalAnalysis\connectome\')
    addpath('D:\Git\st07x2\globalAnalysis\dendriteAnnotation\')
    c = classSmoothConnectome;
    synapses = load('D:\Git\st07x2\globalAnalysis\connectome\LDKMB_syns\globalSynScores.mat');
    synapses.edges = c.id(load('D:\Git\st07x2\globalAnalysis\connectome\LDKMB_syns\globalEdges.mat')).r.edges;
else    
    addpath('/gaba/u/kboerg/Git/st07x2/globalAnalysis/connectome/')
    addpath('/gaba/u/kboerg/Git/st07x2/globalAnalysis/dendriteAnnotation/')
    c = classSmoothConnectome;
    synapses = load('/gaba/u/mberning/results/pipeline/20141007T094904/globalSynScores.mat');
    synapses.edges = c.id(load('/gaba/u/mberning/results/pipeline/20141007T094904/globalEdges.mat')).r.edges;
end
presyns = synapses.edges(synapses.synScores>0);
presynSegments = getXSegments(graph, presyns, 0.97, 49, documentID, Inf, []);