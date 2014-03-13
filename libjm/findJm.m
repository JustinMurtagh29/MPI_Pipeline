function jm = findJm()

jm(1) = findResource('scheduler', 'type', 'jobmanager', 'LookupURL', 'fermat01', 'name', 'fermat-job-manager');
jm(2) = findResource('scheduler', 'type', 'jobmanager', 'LookupURL', 'fermat01', 'name', 'fermat-gpu-manager');

end

