function plot3D( volumeLabel )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
contourListOneCell = volumeLabel.contourList(volumeLabel.contourList(:,3) == 1, 2);
coordsOneCell = cell2mat(volumeLabel.contours(1,contourListOneCell)');
plot3(coordsOneCell(:,1), coordsOneCell(:,2), coordsOneCell(:,3), 'x')

end

