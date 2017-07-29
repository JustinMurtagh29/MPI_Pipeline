function [ idxis,sub2idxis ] = bucketsort( subs )
%BUCKETSORT sorts a set of numbers from 1 to n & returns a cell idxis of size n+1
%containing the indexes all(subs(idxis{val(i)})==i)=true for all i that are
%in subs
    idxis=accumarray(subs,1:numel(subs),[],@(x){x});
    empties=cellfun(@(idxi)isempty(idxi),idxis);
    sub2idxis=cumsum(~empties);
    sub2idxis(empties)=nan;
    idxis=idxis(~empties);
end