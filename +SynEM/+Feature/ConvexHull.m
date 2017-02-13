classdef ConvexHull < SynEM.Feature.ShapeFeature
    %CONVEXHULL Convex hull.
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
    end
    
    methods
        function obj = ConvexHull(numFeat)
            %INPUT numFeat: (Optional) Number of features calculated by
            %           Volume.
            %           (Default: 1)
            
            obj.name = 'ConvHull';
            if exist('numFeat','var') && ~isempty(numFeat)
                obj.numFeatures = numFeat;
            else
                obj.numFeatures = 1;
            end
        end
        
        function feat = calculate(obj, vol, selFeat)
            if ~exist('selFeat','var') || isempty(selFeat)
                selFeat = true(obj.numFeatures,1);
            end
            vol = vol(:,selFeat);
            feat = cellfun(@obj.convHullVol3D, vol);
        end
    end
    methods (Static) 
        function V = convHullVol3D(X)
            %Convex hull for 3d object. If points are coplanar then a 2D
            % convex hull is calculated. Colinear is currently not handled
            % and causes an error.
            % INPUT X: [Nx3] integer vector of shape coordinates.
            
            try
                [~,V] = convhull(X);
            catch
                try
                    %points are coplanar
                    if length(unique(X(:,3))) == 1
                        [~,V] = convhull(X(:,[1, 2]));
                    elseif length(unique(X(:,1))) == 1
                        [~,V] = convhull(X(:,[2, 3]));
                    elseif length(unique(X(:,2))) == 1
                        [~,V] = convhull(X(:,[1 3]));
                    else
                        V = 0;
                    end
                catch
                    %points are colinear
                    V = 0;
                end
            end
        end
    end
    
end

