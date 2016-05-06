function f = plotParamSearchSegResult( result )
%PLOTPARAMSEARCHSEGRESULT
% INPUT result: See Codat.Pipeline.parameterSearchSeg
% OUTPUT f: Figure handle.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%iterate over node thresholds
legendNames = cell(2*length(result.smCurveTrain),1);
for i = 1:length(result.smCurveTrain)
    resultTrain = result.smCurveTrain{i}{:,:};
    resultTest = result.smCurveTest{i}{:,:};
    scatter(resultTrain(:,1),resultTrain(:,2),'filled');
    hold on
    scatter(resultTest(:,1),resultTest(:,2),'filled');
    legendNames{2*i-1} = ['training node threshold ' num2str(i)];
    legendNames{2*i} = ['test node threshold ' num2str(i)];
end
xlabel('average path length between merger [um]');
ylabel('average path length between splits [um]');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
legend(legendNames{:});
grid on
set(gca,'FontSize',24);
set(gca,'YLim',[1e-1,1e4]);

end