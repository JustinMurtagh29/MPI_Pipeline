function [weights] = featureDesign(input, border)
%Comment

weights = [];
for j = 1:size(border,2)
	b = border(j).PixelIdxList;
	temp = real(input(b));
  	weights1(1) = quantile(temp,0);
	weights1(2) = quantile(temp,0.25);
	weights1(3) = quantile(temp,0.5);
     	weights1(4) = quantile(temp,0.75);
       	weights1(5) = quantile(temp,1);
	weights1(6) = mean(temp);
	weights1(7) = std(temp);
	weights = [weights; weights1];
end
end

