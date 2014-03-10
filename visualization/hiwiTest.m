%% Some first results from Hiwi testing
tic;
batch1 = getMissionData(1);
toc
% Exclude missions before some date
dateStart = '18 Sep 2013 16:56:04';
timestampStart = datenum(date, 'dd mmm yyyy HH:MM:SS');
% Create lists from structure
levelList = strcat({batch1(:).levelName}', cellstr(num2str(cat(1,batch1(:).levelCodeVersion), '%.5i')))';
userList = {batch1(:).user};
missionList = [batch1(:).missionId];
timestampList = [batch1(:).missionId];
% Find unique elements of the list
levelUnique = unique(levelList);
userUnique = unique(userList);
missionUnique = unique(missionList);



%%
% Group accordingly
results = cell(length(uniqueMissions), length(uniqueUsers), length(uniqueLevels));
for i=1:length(levelString)
    idxLevel = find(strcmp(levelString{i}, uniqueLevels));
    idxUser = find(strcmp(userString{i}, uniqueUsers));
    idxMission = find(missions(i) == uniqueMissions);
    results{idxMission, idxUser, idxLevel} = [results{idxMission,idxUser,idxLevel} solutions{i}];
end

%% Exclude everything before a certain time stamp
missions = batch1(:)

%% Plot
clc;

for i=1:length(levelUnique)
    for j=1:length(missionUnique)
        if sum(~cellfun(@isempty,results(j,:,i))) > 2
            figure('Units', 'centimeters', 'Position', [0 0 29.7 21], ...
                'PaperType', 'A4', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 29.7 21], 'PaperOrientation', 'portrait');
            hold on;
            votes = [];
            for k=1:size(results,2)
                if ~isempty(results{j,k,i})
                    votes = [votes [results{j,k,i}; repmat(k,1,size(results{j,k,i},2))]];
                end
            end
            [uniqueIDs, m, n] = unique(votes(1,:));
            for l=1:length(uniqueIDs)
                subplot(length(uniqueIDs), 1, l);
                label1 = sum(n==l & votes(2,:) == 1);
                label0 = sum(n==l & votes(2,:) == 0);
                labelM1 = sum(n==l & votes(2,:) == -1);
                plot([-1 0 1], [labelM1 label0 label1], 'xb', 'MarkerSize', 10, 'LineWidth', 3);
                set(gca, 'Xtick', -1:1);
                set(gca, 'XTickLabel', {'not connected' 'no answer' 'connected'});
                xlim([-1.5 1.5]);
                ylim([0 max([label1 label0 labelM1])+1]);
                ylabel('vote counts');
                title(['Statistics for mission: ' num2str(uniqueMissions(j)) ' level: ' uniqueLevels{i} ' ID: ' num2str(uniqueIDs(l))]);
            end
            set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
            print(gcf, ['C:\Users\mberning\Desktop\hiwiEval\mission' num2str(uniqueMissions(j), '%.3i') 'level' uniqueLevels{i} '.pdf'], '-dpdf', '-r300');
            close all;
        end
    end
end

%% Mission coverage
clc;
figure('Units', 'centimeters', 'Position', [0 0 29.7 21], ...
'PaperType', 'A4', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 29.7 21], 'PaperOrientation', 'portrait');
hold on;
solutionPresent = ~cellfun(@isempty, results);
for i=1:size(solutionPresent,3);
    subplot(9,1,i);
    h1 = imagesc(solutionPresent(:,:,i)');
    if i== 5
        ylabel('User ID');
    end
    if i == 9
        xlabel('Mission ID');
    end
    title(['level: ' uniqueLevels{i}]);
    set(gca,'XTickLabel','')
    set(gca,'YTickLabel','')
end
set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
print(gcf, ['C:\Users\mberning\Desktop\hiwiEval\missionCoverage.pdf'], '-dpdf', '-r300');
close all;