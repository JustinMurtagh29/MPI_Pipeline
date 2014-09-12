function weights = normalizeDataForGPold(weights,a,normFile)
   minValues = min(weights,[],1);
   maxValues = max(weights,[],1);
   compFactor = 1./(maxValues-minValues);
   save(normFile, 'minValues', 'maxValues', 'compFactor');
   weights = bsxfun(@minus,weights,minValues);
   weights = bsxfun(@times,weights,compFactor);
end
