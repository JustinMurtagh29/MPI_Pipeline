function newmap = changem(map, newcode, oldcode)

    [~, idx] = ismember(map, oldcode);
    newmap = newcode(idx);
    if all(size(map)~=size(newmap))  % fix only-1-edge problem
        newmap = newmap';
    end
end
