%function generateQueriesFromDirections
% zz = load('/gaba/scratch/kboerg/directionalityTest1.mat');
% directionsMB.result.latent{1006}(:, 1)
% directionsMB.result.pca{1000}(:, :, 1)
% directionsMB.result.neighbours{1000}(1)
% directionsMB.result.prob{1000}(1)
% directionsMB.result.borderIdx{1000}(1)
% directionsMB.result.scores{1000}(1)
direction_col = [];
startpoint_col = [];
memory_col = []
for aggloidx = 1:10000

    todo = directionsMB.result.latent{aggloidx}(:, 1) > 0.7;
    todo = todo & abs(directionsMB.result.scores{aggloidx}) > 0.9;
    if max(pdist(bsxfun(@times, double(borderMeta.borderCoM(directionsMB.result.borderIdx{aggloidx}(todo), :)), [11.24, 11.24, 28]))) < 5000
        continue;
    end
    todo = todo & directionsMB.result.scores{aggloidx} .* squeeze(sign(directionsMB.result.pca{aggloidx}(3, 1, :))) > 0.9;
    todoF = find(todo);
    for idx = 1 : length(todoF)
        direction = 100 * sign(directionsMB.result.scores{aggloidx}(todoF(idx))) * directionsMB.result.pca{aggloidx}(:, 1, todoF(idx))' ./ [11.24, 11.24, 28]; %direction is normalized to 100nm
        startpoint = round(double(borderMeta.borderCoM(directionsMB.result.borderIdx{aggloidx}(todoF(idx)),:)) - 2 * direction);
        if ~ismember(Seg.Global.getSegIds(p, startpoint), directionsMB.axons{1}{aggloidx})
            continue
        end
        aggloidx
        direction_col(end + 1, :) = direction;
        startpoint_col(end + 1, :) = startpoint;
        memory_col(end+1, :) = [aggloidx, idx];
        break
    end
end
projectpath = '/gaba/scratch/kboerg/annotationZips2467361239877643360testFlightQuery20170514_nmls';
mkdir(projectpath)
unzip([projectpath '.zip'], [projectpath filesep]);
liste = dir([projectpath '/*.nml']);
for idx = 1 : length(liste)
    skel = skeleton([projectpath filesep liste(idx).name]);
    assert(length(skel.nodes) == 1);
    find(ismember(startpoint_col, skel.nodes{1}(1, 1 : 3), 'rows'))
end
