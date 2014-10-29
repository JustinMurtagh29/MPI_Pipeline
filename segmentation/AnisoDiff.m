function diff = AnisoDiff( diff, kRange, niter, option)
%AnisoDiff Anisotropic Diffusion Filter with 6 neighbouring voxels
%   

[sz1,sz2,sz3] = size(diff);

    for i=1:niter
            
              diffl = zeros(sz1+2, sz2+2, sz3+2 );
              diffl(2:sz1+1, 2:sz2+1, 2:sz3+1) = diff;

              deltaN = diffl(1:sz1,1:sz2,3:sz3+2)   - diff;
              deltaS = diffl(1:sz1,1:sz2,1:sz3) - diff; 
              deltaE = diffl(3:sz1+2,1:sz2,1:sz3) - diff; 
              deltaW = diffl(1:sz1,1:sz2,1:sz3)   - diff;
              deltaV = diffl(1:sz1,1:sz2,1:sz3)   - diff;
              deltaH = diffl(1:sz1,3:sz2+2,1:sz3)   - diff;

              if option == 1
                cN = exp(-(deltaN/kRange).^2);   
                cS = exp(-(deltaS/kRange).^2);
                cE = exp(-(deltaE/kRange).^2);
                cW = exp(-(deltaW/kRange).^2);
                cV = exp(-(deltaV/kRange).^2);
                cH = exp(-(deltaH/kRange).^2);   

              elseif option == 2
                cN = 1./(1 + (deltaN/kRange).^2);
                cS = 1./(1 + (deltaS/kRange).^2);
                cE = 1./(1 + (deltaE/kRange).^2);
                cW = 1./(1 + (deltaW/kRange).^2);
                cV = 1./(1 + (deltaV/kRange).^2);
                cH = 1./(1 + (deltaH/kRange).^2);

               elseif option == 3
                cN = (0.5*(1-(deltaN/kRange).^2).^2);
                cS = (0.5*(1-(deltaS/kRange).^2).^2);
                cE = (0.5*(1-(deltaE/kRange).^2).^2);
                cW = (0.5*(1-(deltaW/kRange).^2).^2);
                cV = (0.5*(1-(deltaV/kRange).^2).^2);
                cH = (0.5*(1-(deltaH/kRange).^2).^2);
              end

            diff = diff + 0.16*(cN.*deltaN + cS.*deltaS + cE.*deltaE + cW.*deltaW + cV.*deltaV + cH.*deltaH);
    end  
end