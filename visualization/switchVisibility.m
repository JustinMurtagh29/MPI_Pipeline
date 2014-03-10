function switchVisibility( k, cellNr )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for i=1:length(cellNr)
    if strcmp(get(k{cellNr(i)}, 'Visible'), 'on')
        set(k{cellNr(i)}, 'Visible', 'off');
    else
        set(k{cellNr(i)}, 'Visible', 'on');
    end
end

end

