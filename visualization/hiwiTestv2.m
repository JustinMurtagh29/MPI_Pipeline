%% Some first results from Hiwi testing
load('D:\sync\activeTraining\DBreadout20131002.mat');

%% Convert timestamps to MATLAB format (-> transfer to getMissionData.m on cluster)
for i=1:length(data)
    data(i).timestampM = changeTimestampsToMatlabFormat(data(i).timestamp);
end

%% Print timestamps (sanity check)
display(datestr([data(:).timestampM]));


%% Exlude all data points before certain timestamp
startTime = '18 Sep 2013 16:56:04';
startTime = datenum(startTime, 'dd mmm yyyy HH:MM:SS');
selectedData = data([data(:).timestampM] > startTime);

%% Split data according to game
dataS.b4b1 =  selectedData(strcmp('b4b_1', {selectedData(:).game}));
dataS.b4b2 = selectedData(strcmp('b4b_2', {selectedData(:).game}));
dataS.b4b3 = selectedData(strcmp('b4b_3', {selectedData(:).game}));
dataS.tracer = selectedData(strcmp('tracer', {selectedData(:).game}));

%% Level deployment according to Yates
nrMissions.b4b1 = 918;
nrMissions.b4b2 = 90;
nrMissions.b4b3 = 918;
nrMissions.tracer = 541;

%% Mission coverage
% Generate control data and counts of missions
games = fieldnames(dataS);
for i=1:length(games)
    missionNumbers.(games{i}) = single([dataS.(games{i}).missionId]);
    controlData.(games{i}) = randi(nrMissions.(games{i}), length(missionNumbers.(games{i})), 1);
    counts.(games{i}) = histc(single(missionNumbers.(games{i})), 1:1000);
    countsC.(games{i}) = histc(single(controlData.(games{i})), 1:1000);
end

for i=1:length(games)
    figure;
    h = stem(1:1000, [counts.(games{i})' countsC.(games{i})]);
    set(h(1), 'MarkerSize', 3, 'LineStyle', ':', 'LineWidth', .5, 'Color', 'b', 'MarkerEdgeColor', 'b', 'Marker', 'd');
    set(h(2), 'MarkerSize', 3, 'LineStyle', ':', 'LineWidth', .5, 'Color', 'r', 'MarkerEdgeColor', 'r', 'Marker', 's');
    title([games{i} ': Randomness of data: ' num2str(numel(unique(missionNumbers.(games{i})))) ' answered missions; ' num2str(numel(unique(controlData.(games{i})))) ' in control data']);
    xlabel('Mission ID');
    ylabel('Counts');
    legend('Data from missions', 'Random values from Matlab RNG');
    xlim([-0.5 1000.5]);
end

for i=1:length(games)
    figure;
    h = bar(0:7, cat(1,hist(counts.(games{i}),0:7), hist(countsC.(games{i}), 0:7))', 'grouped');
    title([games{i} ': Distibution: Answer count vs control sample']);
    xlabel('Answers to missions');
    ylabel('Counts');
    legend('Data from missions', 'Random values from Matlab RNG');
end
