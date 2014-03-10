function data = findFiles( parentDir )

temp = dir(parentDir);
temp = temp(3:end);
for i=1:length(temp)
    temp2 = dir([parentDir temp(i).name '\']);
    temp2 = temp2(3:end);
    for j=1:length(temp2)
        temp3 = dir([parentDir temp(i).name '\' temp2(j).name '\*.mat']);
        if ~isempty(temp3)
            data(i,j).directory = [parentDir temp(i).name '\' temp2(j).name '\'];
            data(i,j).files = {temp3(:).name};
        end
    end
end

end

