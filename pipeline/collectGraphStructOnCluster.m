function job = collectGraphStructOnCluster(p)
    job = startCPU(@collectGlobalGraphStruct, {{p}}, 'globalGraphStruct');
end
