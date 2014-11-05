function h = cubeAverage(siz)
    h{1} = ones(1,siz,1)/prod([siz]);
    h{2} = ones(siz,1,1)/prod([siz]);
    h{3} = ones(1,1,siz)/prod([siz]);
end