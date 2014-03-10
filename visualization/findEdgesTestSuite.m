load('R:\forRonja\forRonja.mat', 'seg');
%segTest = seg(257:868,257:868,129:384);
segTest = seg(1:50,1:50,1:50);

%%
profile on;
edges1 = findedgesversuch(segTest, 26);
edges2 = findEdges(segTest);
profile viewer;
profile off;

%%
testSizes = [];
for i=1:10