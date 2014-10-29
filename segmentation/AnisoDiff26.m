function diff = AnisoDiff26( diff, kRange, niter, option)
%AnisoDiff Anisotropic Diffusion Filter with 26 neighbouring voxels
%   not memory efficient

[sz1,sz2,sz3] = size(diff);

    for i = 1:niter

      diffl = zeros(sz1+2, sz2+2, sz3+2 );
      diffl(2:sz1+1, 2:sz2+1, 2:sz3+1) = diff;

      deltaN = diffl(2:sz1+1,2:sz2+1,3:sz3+2)   - diff;                        % I_zz
      deltaNV = diffl(2:sz1+1,3:sz2+2,3:sz3+2)   - diff;                       % I_zy
      deltaNH = diffl(2:sz1+1,3:sz2+2,3:sz3+2)   - diff;                       % I_z(-y)
      deltaNE = diffl(3:sz1+2,2:sz2+1,3:sz3+2)   - diff;                       % I_zx
      deltaNW = diffl(1:sz1,2:sz2+1,3:sz3+2)   - diff;                         % I_z(-x)

      deltaNVW = diffl(3:sz1+2,2:sz2+1,3:sz3+2)   - diff;                      % I_zyx
      deltaNVE = diffl(3:sz1+2,1:sz2,3:sz3+2)   - diff;                        % I_zy(-x)
      deltaNHE = diffl(3:sz1+2,3:sz2+2,3:sz3+2)   - diff;                      % I_z(-y)x
      deltaNHW = diffl(1:sz1,3:sz2+2,3:sz3+2)   - diff;                        % I_z(-y)(-x)

      deltaSVW = diffl(1:sz1,1:sz2,1:sz3)   - diff;                            % I_zyx
      deltaSVE = diffl(3:sz1+2,1:sz2,1:sz3)   - diff;                          % I_zy(-x)
      deltaSHE = diffl(3:sz1+2,3:sz2+2,1:sz3)   - diff;                        % I_z(-y)x
      deltaSHW = diffl(1:sz1,3:sz2+2,1:sz3)   - diff;                          % I_z(-y)(-x)

      deltaS = diffl(2:sz1+1,2:sz2+1,1:sz3) - diff;                            % I_(-z)(-z)
      deltaSV = diffl(2:sz1+1,1:sz2,1:sz3)   - diff;                           % I_(-z)y
      deltaSH = diffl(2:sz1+1,3:sz2+2,1:sz3)   - diff;                         % I_(-z)(-y)
      deltaSE = diffl(3:sz1+2,2:sz2+1,1:sz3)   - diff;                         % I_(-z)x
      deltaSW = diffl(1:sz1,2:sz2+1,1:sz3)   - diff;                           % I_(-z)(-x)

      deltaE = diffl(3:sz1+2,2:sz2+1,2:sz3+1) - diff;                          % I_xx
      deltaEH = diffl(3:sz1+2,3:sz2+2,2:sz3+1) - diff;                         % I_x(-y)
      deltaEV = diffl(3:sz1+2,1:sz2,2:sz3+1) - diff;                           % I_xy

      deltaW = diffl(1:sz1,2:sz2+1,2:sz3+1)   - diff;                          % I_(-x)(-x)
      deltaWH = diffl(1:sz1,3:sz2+2,2:sz3+1)   - diff;                         % I_(-x)(-y)
      deltaWV = diffl(1:sz1,1:sz2,2:sz3+1)   - diff;                           % I_(-x)(y)

      deltaV = diffl(2:sz1+1,1:sz2,2:sz3+1)   - diff;                          % I_yy
      deltaH = diffl(2:sz1+1,3:sz2+2,2:sz3+1)   - diff;                        % I_(-y)(-y)

      if option == 1
        cN = exp(-(deltaN/kRange).^2);   
        cS = exp(-(deltaS/kRange).^2);
        cE = exp(-(deltaE/kRange).^2);
        cW = exp(-(deltaW/kRange).^2);
        cV = exp(-(deltaV/kRange).^2);
        cH = exp(-(deltaH/kRange).^2);
        cNH = exp(-(deltaNH/kRange).^2); 
        cNV = exp(-(deltaNV/kRange).^2); 
        cNE = exp(-(deltaNE/kRange).^2); 
        cNW = exp(-(deltaNW/kRange).^2);
        cSH = exp(-(deltaSH/kRange).^2); 
        cSV = exp(-(deltaSV/kRange).^2); 
        cSE = exp(-(deltaSE/kRange).^2); 
        cSW = exp(-(deltaSW/kRange).^2);
        cWH = exp(-(deltaWH/kRange).^2); 
        cWV = exp(-(deltaWV/kRange).^2); 
        cEH = exp(-(deltaEH/kRange).^2); 
        cEV = exp(-(deltaEV/kRange).^2);
        cNVE = exp(-(deltaNVE/kRange).^2); 
        cNVW = exp(-(deltaNVW/kRange).^2); 
        cNHE = exp(-(deltaNHE/kRange).^2); 
        cNHW = exp(-(deltaNHW/kRange).^2); 
        cSVE = exp(-(deltaSVE/kRange).^2); 
        cSVW = exp(-(deltaSVW/kRange).^2); 
        cSHE = exp(-(deltaSHE/kRange).^2); 
        cSHW = exp(-(deltaSHW/kRange).^2);

      elseif option == 2
        cN = 1./(1 + (deltaN/kRange).^2);
        cS = 1./(1 + (deltaS/kRange).^2);
        cE = 1./(1 + (deltaE/kRange).^2);
        cW = 1./(1 + (deltaW/kRange).^2);
        cV = 1./(1 + (deltaV/kRange).^2);
        cH = 1./(1 + (deltaH/kRange).^2);
        cNH = 1./(1 + (deltaNH/kRange).^2);
        cNV = 1./(1 + (deltaNV/kRange).^2);
        cNE = 1./(1 + (deltaNE/kRange).^2);
        cNW = 1./(1 + (deltaNW/kRange).^2);
        cSH = 1./(1 + (deltaSH/kRange).^2);
        cSV = 1./(1 + (deltaSV/kRange).^2);
        cSE = 1./(1 + (deltaSE/kRange).^2);
        cSW = 1./(1 + (deltaSW/kRange).^2);
        cWH = 1./(1 + (deltaWH/kRange).^2);
        cWV = 1./(1 + (deltaWV/kRange).^2);
        cEH = 1./(1 + (deltaEH/kRange).^2);
        cEV = 1./(1 + (deltaEV/kRange).^2);
        cNVE = 1./(1 + (deltaNVE/kRange).^2);
        cNVW = 1./(1 + (deltaNVW/kRange).^2);
        cNHE = 1./(1 + (deltaNHE/kRange).^2);
        cNHW = 1./(1 + (deltaNHW/kRange).^2);
        cSVE = 1./(1 + (deltaSVE/kRange).^2);
        cSVW = 1./(1 + (deltaSVW/kRange).^2);
        cSHE = 1./(1 + (deltaSHE/kRange).^2);
        cSHW = 1./(1 + (deltaSHW/kRange).^2);

       elseif option == 3
           cN = (0.5*(1-(deltaN/kRange).^2).^2);
           cS = (0.5*(1-(deltaS/kRange).^2).^2);
           cE = (0.5*(1-(deltaE/kRange).^2).^2);
           cW = (0.5*(1-(deltaW/kRange).^2).^2);
           cV = (0.5*(1-(deltaV/kRange).^2).^2);
           cH = (0.5*(1-(deltaH/kRange).^2).^2);
           cNH = (0.5*(1-(deltaNH/kRange).^2).^2);
           cNV = (0.5*(1-(deltaNV/kRange).^2).^2);
           cNE = (0.5*(1-(deltaNE/kRange).^2).^2);
           cNW = (0.5*(1-(deltaNW/kRange).^2).^2);
           cSH = (0.5*(1-(deltaSH/kRange).^2).^2);
           cSV = (0.5*(1-(deltaSV/kRange).^2).^2);
           cSE = (0.5*(1-(deltaSE/kRange).^2).^2);
           cSW = (0.5*(1-(deltaSW/kRange).^2).^2);
           cWH = (0.5*(1-(deltaWH/kRange).^2).^2);
           cWV = (0.5*(1-(deltaWV/kRange).^2).^2);
           cEH = (0.5*(1-(deltaEH/kRange).^2).^2);
           cEV = (0.5*(1-(deltaEV/kRange).^2).^2);
           cNVE = (0.5*(1-(deltaNVE/kRange).^2).^2);
           cNVW = (0.5*(1-(deltaNVW/kRange).^2).^2);
           cNHE = (0.5*(1-(deltaNHE/kRange).^2).^2);
           cNHW = (0.5*(1-(deltaNHW/kRange).^2).^2);
           cSVE = (0.5*(1-(deltaSVE/kRange).^2).^2);
           cSVW = (0.5*(1-(deltaSVW/kRange).^2).^2);
           cSHE = (0.5*(1-(deltaSHE/kRange).^2).^2);
           cSHW = (0.5*(1-(deltaSHW/kRange).^2).^2);
      end

        diff = diff + 0.0384*(cN.*deltaN + cS.*deltaS + cE.*deltaE + cW.*deltaW + cV.*deltaV + cH.*deltaH ...
         + cNV.*deltaNV+ cNH.*deltaNH+ cNE.*deltaNE+ cNW.*deltaNW ...
         + cSV.*deltaSV+ cSH.*deltaSH+ cSE.*deltaSE+ cSW.*deltaSW ...
         + cWH.*deltaWH+ cWV.*deltaWV+ cEH.*deltaEH + cEV.*deltaEV ...
         + cNVE.*deltaNVE+ cNVW.*deltaNVW+ cNHE.*deltaNHE + cNHW.*deltaNHW ...
         + cSVE.*deltaSVE+ cSVW.*deltaSVW+ cSHE.*deltaSHE + cSHW.*deltaSHW);

    end
end