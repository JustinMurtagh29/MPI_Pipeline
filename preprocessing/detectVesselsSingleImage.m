function vessel = detectVesselsSingleImage(raw)

    vessel = bwareaopen(raw > 162 | raw < 50, 1000, 4);
	vessel = imclose(vessel, ones(5,5));
	% Do we need to change location of hole drilling? -> checked, should be fine now
    vessel = padarray(vessel, [1 1], 1);
    holeLoc = round(size(vessel,1)/2);
    idx = 1;
	while vessel(holeLoc,idx) == 1
		vessel(holeLoc,idx) = 0;
		idx = idx + 1;
	end
	vessel = imfill(vessel, 'holes');
    vessel = vessel(2:end-1,2:end-1);
end

