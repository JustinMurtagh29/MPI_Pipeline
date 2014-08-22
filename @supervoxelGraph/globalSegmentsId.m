function globalSegId(p)

clear coords    
%construct coord array for bBox Border area
coords{1} = [257,768;257,768;129,384];
coords = repmat(coords,[size(p.local,1),size(p.local,2),size(p.local,3)]);
coords(1,:,:) = cellfun(@(x)([x(1,1)+p.tileBorder(1,1) x(1,2); x(2,1) x(2,2); x(3,1) x(3,2)]),coords(1,:,:), 'UniformOutput', false);
coords(end,:,:) = cellfun(@(x)([x(1,1) x(1,2)+p.tileBorder(1,2); x(2,1) x(2,2); x(3,1) x(3,2)]),coords(end,:,:), 'UniformOutput', false);
coords(:,1,:) = cellfun(@(x)([x(1,1) x(1,2); x(2,1)+p.tileBorder(2,1) x(2,2); x(3,1) x(3,2)]),coords(:,1,:), 'UniformOutput', false);
coords(:,end,:) = cellfun(@(x)([x(1,1) x(1,2); x(2,1) x(2,2)+p.tileBorder(2,2); x(3,1) x(3,2)]),coords(:,end,:), 'UniformOutput', false);
coords(:,:,1) = cellfun(@(x)([x(1,1) x(1,2); x(2,1) x(2,2); x(3,1)+p.tileBorder(3,1) x(3,2)]),coords(:,:,1), 'UniformOutput', false);
coords(:,:,end) = cellfun(@(x)([x(1,1) x(1,2); x(2,1) x(2,2); x(3,1) x(3,2)+p.tileBorder(3,2)]),coords(:,:,end), 'UniformOutput', false);

%%  globalize segmentaiton Ids

numEl = 0;
numElTotal = [0,0,0];

	for i=1:size(p.local,1)
		for j=1:size(p.local,2)
			for k= 1:size(p.local,3)

				load(p.local(i,j,k).segFile)
				seg = seg(coords{i,j,k}(1,1):coords{i,j,k}(1,2), coords{i,j,k}(2,1):coords{i,j,k}(2,2), coords{i,j,k}(3,1):coords{i,j,k}(3,2));
			    seg = uint32(seg + numEl);				
			    numEl = numEl + length(unique(seg(:)));
                numElTotal(i,j,k) = numEl;
                
			    writeKnossosRoi(p.seg.segFolder, p.seg.prefix , [p.local(i,j,k).bboxBig(:,1) + coords{i,j,k}(:,1)]', seg, 'uint32')
				
			end
		end
	end
numElTotal = numElTotal - numElTotal(1,1,1);
save([p.seg.segFolder 'numEl.mat'], 'numElTotal') 

end

