function missionData = getMissionData( batchID )
% Load data from section BatchID passed as input argument and supllies
% a subset of the data in a clean format as first input argument, but
% also if a second output argument is requested the whole structure
root = '/zdata/manuel/results/missions/2012-09-28_ex145_07x2/';
allData = readJson([root 'batch' num2str(batchID) '.json']);
missionData = struct();
% Linearize all entriess
idx = 1;
for i=1:length(allData)
    for j=1:length(allData{i}.renderedStacks)
        for k=1:length(allData{i}.renderedStacks{j}.solutions)
            temp = [allData{i}.renderedStacks{j}.solutions{k}.solution{:}];
            if ~isempty(temp)
                % For unique identification of problem 
                missionData(idx).dataset = allData{i}.dataSetName;
                missionData(idx).batchId = allData{i}.batchId;
                missionData(idx).missionId = allData{i}.missionId;
                % Which ID's was the problem rendered for
                missionData(idx).startId = allData{i}.start.id;
                missionData(idx).endId = allData{i}.end.id;
                % Additional information from missions.json
                missionData(idx).errorCenter = [allData{i}.errorCenter{:}];
                missionData(idx).start = allData{i}.start;
                missionData(idx).isControl = allData{i}.isControl;
                missionData(idx).difficulty = allData{i}.difficulty;
                % Solutions proposed by the active classifier
                missionData(idx).possibleEnds = [allData{i}.possibleEnds{:}];
                % Render & play status in Mongo DB
                missionData(idx).missionStatus = allData{i}.missionStatus;
                % Level & code version used for rendering
                missionData(idx).levelName = allData{i}.renderedStacks{j}.levelId.name;
                missionData(idx).levelCodeVersion = allData{i}.renderedStacks{j}.levelId.version;
                missionData(idx).levelParameter = allData{i}.renderedStacks{j}.paraInfo;
                % User, Game, timestamp and provided solution(s) for each
                missionData(idx).game = allData{i}.renderedStacks{j}.solutions{k}.game;
                missionData(idx).user = allData{i}.renderedStacks{j}.solutions{k}.userId.id;
                missionData(idx).timestamp = allData{i}.renderedStacks{j}.solutions{k}.timestamp;
                missionData(idx).id = [temp(:).segmentId];
                missionData(idx).label = [temp(:).solution];
                idx = idx + 1;
            end
        end
    end
end

end
