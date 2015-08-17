function out = normalizeStack( in )
%out = normalizeStack( in ) Subtract mean and devide by std
%for st08x2
out = in-162;
out = out/22;


%out = in-122;
%out = out./22; 
%out = in-132; % 01.04.2013, quick hack for new dataset (add to parameter construct?)
%out = out./43;
%out = in-155; % 26.03.2015, another quick hack, for ED stack
%out = out./57;
end

