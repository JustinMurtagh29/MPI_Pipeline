
    batchFolder = '/tmpscratch/scchr/AxonEndings/LeftEndings/';

    % Write first batch to tasks on webKnossos
    % This needs the anaconda/2 module loaded before starting matlab to work (on gaba)
    responses = cell(0);
    for i=1:100
        batch = load(fullfile(batchFolder, ['batch' num2str(i, '%.4i') '.mat']));
        tic;
        for j=1:numel(batch.q.pos)
            for k=1:numel(batch.q.pos{j})
                command = sprintf('python ../wk-paper/RESCOPaaS/createTaskWithScript.py %u %u %u %.2f %.2f %.2f 1 CS_MB_L4_AxonLeftQueries', ...
                    batch.q.pos{j}{k}(:), batch.q.angles{j}{k}(:));
                [~,response] = system(command);
                responses{i,j,k} = response;
                pause(.1);
            end
        end
        toc;
    end
    save(fullfile(batchFolder, 'responses.mat'));

