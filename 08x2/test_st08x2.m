addpath auxiliaryMethods/cubes

parameter.bbox = [257+512 9984-2048; 769+1024+512 7936-1024-512; 129+128 7552-128]; % this should be aligned with KNOSSOS cubes and be divisble by tileSize
parameter.tileSize =  [512; 512; 256]; % Size of local segmentation and local graph construction
parameter.tileBorder = [-256 256; -256 256; -128 128]; % border of local segmentation included for gloablization and large size due to games
parameter.tiles = (parameter.bbox(:,2) - parameter.bbox(:,1) + 1) ./ parameter.tileSize;
%roi=readKnossosRoi('Z:\Data\R-Drive\2012-11-23_ex144_st08x2\color\16\','2012-11-23_ex144_st08x2_mag16',[1 128*5;1 128*4;1 128*4]);
%implay(roi)
%roi2=zeros(size(roi),'uint8')+255;
roi2=roi/2+128;
for i=1:parameter.tiles(1)
    for j=1:parameter.tiles(2)
        for k=1:parameter.tiles(3)
            parameter.local(i,j,k).bboxSmall = round((k-28)/2.9)*[0 0;128 128;0 0]+[parameter.bbox(:,1) + [i-1; j-1; k-1] .* parameter.tileSize parameter.bbox(:,1) + [i; j; k] .* parameter.tileSize - [1; 1; 1]];
            parameter.local(i,j,k).bboxBig = parameter.local(i,j,k).bboxSmall + parameter.tileBorder;
            z=round(parameter.local(i,j,k).bboxBig/16);
            roi2(z(1,1)+1:z(1,2),z(2,1)+1:z(2,2),z(3,1)+1:z(3,2))=...
                roi(z(1,1)+1:z(1,2),z(2,1)+1:z(2,2),z(3,1)+1:z(3,2));
        end
    end
end

implay(roi2(:,:,1:4:end))
