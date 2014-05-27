function weights = featureDesign(input, border)
%Comment

weights = zeros(size(border,1),5*length(input));
for i=1:size(border,1)
	for j=1:length(input)
		temp = real(input{j}(cell2mat(border)));
  		weights(i,(j-1)*7+1) = quantile(temp,0);
		weights(i,(j-1)*7+2) = quantile(temp,0.25);
        weights(i,(j-1)*7+3) = quantile(temp,0.5);
        weights(i,(j-1)*7+4) = quantile(temp,0.75);
        weights(i,(j-1)*7+5) = quantile(temp,1);
		weights(i,(j-1)*7+6) = mean(temp);
		weights(i,(j-1)*7+7) = std(temp);
	end
end

end

