function X_norm = normalizeTrainingDataForGP(X, renormalize, normFile)

percentiles = [.1 .9];

if renormalize
    X_lower_cutoff = quantile(X,percentiles(1),1);
    X_upper_cutoff = quantile(X,percentiles(2),1);
else
    load(normFile);
end

X = bsxfun(@minus, X, X_lower_cutoff);
X_norm = bsxfun(@times, X, 1./(X_upper_cutoff-X_lower_cutoff));

if renormalize
    save(normFile, 'X_lower_cutoff', 'X_upper_cutoff');
end

end

