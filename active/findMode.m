function [ output_args ] = findMode( K, y, lik )

[alpha, sW, L, nlZ, dnlZ] = feval('approxLA', K, lik, x, y);
while Psi_old - Psi_new > tol && it < maxIterations
    it = it+1;
    Psi_old = Psi_new;
    sW = sqrt(W);
    L = chol(eye(n)+sW*sW'.*K);
    b = W.*(f-m) + dlp;
fminsearch(
end

end

function [alpha, sW, L, nlZ, dnlZ] = approxLA(K, lik, x, y)

persistent best_alpha best_nlZ
tol = 1e-9;
n = size(x,1);

if any(size(best_alpha) ~= [n,1])   % find a good starting point for alpha and f
  f = zeros(n,1); alpha = f;                                     % start at zero
  [lp,dlp,d2lp] = feval(lik,y,f,'deriv');   W=-d2lp;
  Psi_new = lp; best_nlZ = Inf; 
else
  alpha = best_alpha; f = K*alpha;                             % try best so far
  [lp,dlp,d2lp] = feval(lik,y,f,'deriv');   W=-d2lp;
  Psi_new = -alpha'*f/2 + lp;         
  if Psi_new < -n*log(2)                                 % if zero is better ..
    f = zeros(n,1); alpha = f;                                      % .. go back
    [lp,dlp,d2lp] = feval(lik,y,f,'deriv'); W=-d2lp; 
    Psi_new = -alpha'*f/2 + lp;
  end
end
Psi_old = -Inf;                                    % make sure while loop starts

while Psi_new - Psi_old > tol                        % begin Newton's iterations
  Psi_old = Psi_new; alpha_old = alpha; 
  sW = sqrt(W);                     
  L = chol(eye(n)+sW*sW'.*K);                            % L'*L=B=eye(n)+sW*K*sW
  b = W.*f+dlp;
  alpha = b - sW.*solve_chol(L,sW.*(K*b));
  f = K*alpha;
  [lp,dlp,d2lp,d3lp] = feval(lik,y,f,'deriv'); W=-d2lp;

  Psi_new = -alpha'*f/2 + lp;
  i = 0;
  while i < 10 && Psi_new < Psi_old               % if objective didn't increase
    alpha = (alpha_old+alpha)/2;                      % reduce step size by half
    f = K*alpha;
    [lp,dlp,d2lp,d3lp] = feval(lik,y,f,'deriv'); W=-d2lp;
    Psi_new = -alpha'*f/2 + lp;
    i = i+1;
  end
end

sW = sqrt(W);                                                    % recalculate L
L  = chol(eye(n)+sW*sW'.*K);                             % L'*L=B=eye(n)+sW*K*sW
nlZ = alpha'*f/2 - lp + sum(log(diag(L)));      % approx neg log marg likelihood
    
if nlZ < best_nlZ
  best_alpha = alpha; best_nlZ = nlZ;
end

end

function [lik, logLik, sumLogLik, d1sumLogLik, d2sumLogLik, d3sumLogLik] = cumGaussian(y, f)
% Function to calculate likelihood, log likelihood and derivatives
    yf = y.*f;
    lik = (1+erf(yf/sqrt(2)))/2;
    if nargout > 1
        logLik = log(lik);
        if nargout > 2
            sumLogLik   = sum(logLik);
            if nargout > 3
                n_p = exp(-yf.^2/2)/sqrt(2*pi)./lik;
                d1sumLogLik = y.*n_p;
                if nargout > 4
                    d2sumLogLik = -n_p.^2 - yf.*n_p;
                    if nargout > 5
                        d3sumLogLik = 2*y.*n_p.^3 +3*f.*n_p.^2 +y.*(f.^2-1).*n_p; 
                    end
                end
            end
        end
    end
end

function [lik, logLik, sumLogLik, d1sumLogLik, d2sumLogLik, d3sumLogLik] = logisitic(y, f)
	yf = y.*f;
	lik = 1./(1+exp(-yf));
    if nargout > 1
        logLik = log(lik);
        if nargout > 2
            s    = -yf;
            ps   = max(0,s);
            sumLogLik = -sum(ps+log(exp(-ps)+exp(s-ps)));
            if nargout > 3
                s    = min(0,f); 
                p    = exp(s)./(exp(s)+exp(s-f));
                d1sumLogLik = (y+1)/2-p;
                if nargout>4
                    d2sumLogLik = -exp(2*s-f)./(exp(s)+exp(s-f)).^2;
                    if nargout > 5
                        d3sumLogLik = 2*d2sumLogLik.*(0.5-p);
                    end
                end
            end
        end
    end
end

