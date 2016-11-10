function vessel = detectVesselsSingleImage(raw)

    temp = bwareaopen(raw > 162 | raw < 50, 1000, 4);
	temp = imclose(temp, ones(5,5));
    % Remove still small sections afterwards (to remove remaining detection in bright myelinated axons
    temp = bwareaopen(temp, 2000, 4);
	% Do we need to change location of hole drilling? -> checked, should be fine now
    temp = padarray(temp, [1 1], 1);
    idx = 1;
	while temp(1000,idx) == 1
		temp(1000,idx) = 0;
		idx = idx + 1;
	end
	temp = imfill(temp, 'holes');
	temp([1 end], :) = [];
	temp(:, [1 end]) = [];
    vessel = temp;

end

