function timestampMATLAB = changeTimestampsToMatlabFormat( timestampDB )
% Transform between timestamp formats

% Unix uses 1 Jan 1970 00:00:00 as reference
timestampMATLAB = datenum([1970 1 1 0 0 double(timestampDB/1000)]);

end
