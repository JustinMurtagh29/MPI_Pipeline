function result = analyzeGraph( graph )

% Size measures
result.order = length(graph.edges);
result.size = sum(graph.edges(:) ~= 0);
% Degree Distibution
neighbours = full(sum(graph.edges~=0,2));
result.degreeMean = mean(neighbours);
result.degreeHist = histc(neighbours, unique(neighbours));
% Connected Components
diagEdges = graph.edges|speye(size(graph.edges));
[~,p,~,r] = dmperm(diagEdges~=0);
result.comp_sizes = diff(r);
result.cc = numel(result.comp_sizes);
result.comps = zeros(1,result.order); 
result.comps(r(1:result.cc)) = ones(1,result.cc);
result.comps = cumsum(result.comps);
result.comps(p) = result.comps;
% Node Strength
result.nodeStrengthP = sum( graph.edges.*(graph.edges>0));
result.nodeStrengthN = sum(-graph.edges.*(graph.edges<0));
result.nodeStrengthPsum = sum(result.nodeStrengthP);
result.nodeStrengthNsum = sum(result.nodeStrengthN);
% Clustering Coefficent
nonZero=sum(graph.edges~=0);
cyc3=diag((graph.edges.^(1/3))^3);
nonZero(cyc3==0)=inf;
result.clustering=cyc3'./(nonZero.*(nonZero-1));
% Characteristic path length
D = graph.edges + 1.8;
result.charLength = sum(sum(D(D~=Inf)))/length(nonzeros(D~=Inf));
result.eccentricity = max(D.*(D~=Inf),[],2);
result.radius = min(result.eccentricity);
result.diameter = max(result.eccentricity);
n = size(D,1);
D = 1./D;
D(1:n+1:end) = 0;
result.efficiency = sum(D(:))/(n*(n-1));

end

