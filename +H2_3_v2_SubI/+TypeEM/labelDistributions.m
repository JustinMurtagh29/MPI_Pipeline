% get label distributions in training and test data

numGlia = sum(gtTrain.label(:,1)>0);
numAxon = sum(gtTrain.label(:,2)>0);
numDend = sum(gtTrain.label(:,3)>0);
trainStatNum = [numGlia, numAxon, numDend];
trainStatFrac = trainStatNum./size(gtTrain.label,1);

numGlia = sum(gtTest.label(:,1)>0);
numAxon = sum(gtTest.label(:,2)>0);
numDend = sum(gtTest.label(:,3)>0);
testStatNum = [numGlia, numAxon, numDend];
testStatFrac = testStatNum./size(gtTest.label,1);

figure;
subplot(1,2,1)
x = [1,2];
y = [trainStatNum;testStatNum];
barh(x,y,'grouped');
set(gca,'yticklabels',{'train','test'})
subplot(1,2,2)
x = [1,2];
y = [trainStatFrac;testStatFrac];
barh(x,y,'grouped');
set(gca,'yticklabels',{'train','test'})
legend({'glia','axon','dendrite'})
saveas(gcf,fullfile(param.saveFolder,'typeEM','train_test_statistics.png'))
close all



