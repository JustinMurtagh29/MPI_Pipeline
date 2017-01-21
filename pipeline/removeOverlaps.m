function job = removeOverlaps(p)

    % Transfer local cubes to segFile
    for i=1:numel(p.local)
        inputCell{i} = {p.local(i).tempSegFile, p.local(i).segFile, p.local(i).bboxSmall, p.local(i).bboxBig};
    end
    functionH = @rewriteWithoutOverlaps;
    job = startCPU(functionH, inputCell, 'overlapRemoval', 12, 10);

end

