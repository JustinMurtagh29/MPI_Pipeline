function data = findFilesCortex( parentDir )

temp = dir(parentDir);
temp = temp(3:end);
for i=1:length(temp)
    temp2 = dir([parentDir temp(i).name '\parameter*']);
    if ~isempty(temp2)
        data{i} = [parentDir temp(i).name '\' temp2(end).name];
    end
end

end

