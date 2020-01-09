% get label distributions in training and test data
numSH = sum(gtTrain.label(:,1)>0);
numNSH = sum(gtTrain.label(:,1)<0);
trainStatNum = [numSH, numNSH];
trainStatFrac = trainStatNum./size(gtTrain.label,1);

numSH = sum(gtTest.label(:,1)>0);
numNSH = sum(gtTest.label(:,1)<0);
testStatNum = [numSH, numNSH];
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
legend({'SH','nonSH'})
saveas(gcf,fullfile(param.saveFolder,'typeEM','spine','train_test_statistics_dense.png'))
close all



