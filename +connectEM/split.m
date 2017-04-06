
% Try to find and split merged agglomerates on current state
load('/home/mberning/Desktop/20170326T220813_agglomeration/initialAgglo.mat');

for i=1:length(initialPartition)
    [value(i), goodWeight(i), badWeight(i), nrEdgesGood(i), nrEdgesBad(i)] = connectEM.aggloObjective2(graph, initialPartition{i});
    display(['Component ' num2str(i, '%.4i')]);
end

%% Plot 

figure;
subplot(2,2,1);
plot(cat(1, goodWeight, abs(badWeight))');
legend('good weight', 'bad weight');
xlabel('agglomerate ID');
ylabel('value (see legend)');
subplot(2,2,2);
plot(cat(1, nrEdgesGood, nrEdgesBad)');
legend('#edges good weight', '#edges bad weight');
xlabel('agglomerate ID');
ylabel('value (see legend)');
subplot(2,2,3);
plot(cat(1, goodWeight./nrEdgesGood, badWeight./nrEdgesBad)');
legend('average good weight', 'average bad weight');
xlabel('agglomerate ID');
ylabel('value (see legend)');
subplot(2,2,4);
plot(cat(1, nrEdgesGood./(nrEdgesGood+nrEdgesBad), nrEdgesBad./(nrEdgesGood+nrEdgesBad))');
legend('fraction edges good weight', 'fraction edges bad weight');
xlabel('agglomerate ID');
ylabel('value (see legend)');
