function w = gaussianRegression( w, sigmaSquared, X, Y)

    % regression for w
    A = inv(1/sigmaSquared .* (X*X') + inv(w.Sigma));
    w.Mu = w.Mu + 1/sigmaSquared .* A * X * (Y - w.Mu'*X)';
    w.Sigma = A;
    
end
