% check if agglo (in points representation) candidate has more than 
%just the two contact points at the x bounds
function b = moreContactPoints(agglo, foundPts, tolx, toly, tolz, xbounds, ybounds, zbounds)

b = 0;
for i=1:length(agglo)
    p = agglo(i, :);
    x = agglo(i, 1);
    y = agglo(i, 2);
    z = agglo(i, 3);
    % if border is touched at both places within a agglo
    if y <= (ybounds(1) + toly) || y >= (ybounds(2) - toly)
       b = 1;
       break;
    end
    if z <= (zbounds(1) + tolz) || z >= (zbounds(2) - tolz)
       b = 1;
       break;
    end
    % check for multiple in x
    % not working properly
%     if x <= (xbounds(1) + tolx) && (pdist2(p, foundPts(1,:),'euclidean') > 1000)
%         %disp(length(agglo));
%         b = 1;
%         break;
%     elseif x >= (xbounds(2) - tolx) && (pdist2(p, foundPts(2,:),'euclidean') > 1000)
%         %disp(length(agglo));
%         b = 1;
%         break;
%     end
    
end


end