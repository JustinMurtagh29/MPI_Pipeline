function [seg,newEdges,borders] = correctLeaves(seg,leaves,newEdges,borders)

for i=1:size(leaves,1)
    %Change values of leaves to the values of their parent nodes
    seg(ismember(seg,leaves{i,2}))= leaves{i,1};

    %Find linear indices of newedges with leaves as nodes (use modulo to
    %derive the rownumbers, as intersect gives just the linear indices of
    %the nx2 matrix newedges
    [~,temp]=ismember(newEdges,leaves{i,2});
    temp=find((temp(:,1)+temp(:,2))~=0);
    
    %Set all borderpixel of leaves to the value of the parent node
    indices = borders(temp);
    seg(cat(2,indices(:).PixelIdxList))=leaves{i,1};
    
    %Delete borders and newedges of leave
    borders(temp)=[];
    newEdges(temp,:)=[];
end

end
