function newmap = changem(map, newcode, oldcode)

    [~, idx] = ismember(map, oldcode);
    newmap = newcode(idx);

end
