function vessel = detectVesselSingleImage(raw)

    temp = bwareaopen(raw > 162 | raw < 50, 1000, 4);
	temp = imclose(temp, ones(5,5));	
	% Do we need to change location of hole drilling? -> test
    temp = padarray(temp, [1 1], 1);
    idx = 1;
	while temp(1800,idx) == 1
		temp(1800,idx) = 0;
		idx = idx + 1;
	end
	temp = imfill(temp, 'holes');
	temp([1 end], :) = [];
	temp(:, [1 end]) = [];
    vessel = temp;

end

