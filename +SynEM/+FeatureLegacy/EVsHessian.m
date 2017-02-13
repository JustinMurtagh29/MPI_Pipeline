classdef EVsHessian < SynEM.Feature.TextureFeature
    %EVSHESSIAN Eigenvalues for hessian sorted by increasing absolute
    % absolute value.
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

    properties
        sigma
        filterSiz
    end

    methods
        function obj = EVsHessian(sigma, filterSiz)
            obj.name = 'EVsHessian';
            obj.sigma = sigma;
            obj.numChannels = 3;
            obj.filterSiz = filterSiz;
            obj.border = 2.*filterSiz;
        end

        function fm = calc(obj, raw)
            
            %calculate hessian
            sigma = obj.sigma; %#ok<*PROPLC>
            filterSiz = obj.filterSiz;
            sigma = num2cell(sigma);
            coords = arrayfun(@(siz)(-siz:siz).',filterSiz, ...
                'UniformOutput',false);
            coords = cellfun(@(coords,dim) ...
                permute(coords,[2:dim,1,(dim+1):3]), ...
                coords,num2cell(1:3),'UniformOutput',false);
            gauss = cellfun(@(coords,sigma)exp(-(coords./sigma).^2./2), ...
                coords,sigma,'UniformOutput',false);
            gauss = cellfun(@(gauss)gauss./sum(gauss),gauss, ...
                'UniformOutput',false);
            gaussD = cellfun(@(coords,gauss,sigma) ...
                -gauss.*coords./(sigma.^2),coords,gauss,sigma, ...
                'UniformOutput',false);
            gaussD2 = cellfun(@(coords,gauss,sigma) ...
                gauss.*(coords.^2./(sigma.^2)-1)./(sigma.^2), ...
                coords,gauss,sigma,'UniformOutput',false);
            H = cell(3,3);
            for dim1 = 1:3
                for dim2 = 1:dim1
                    if(dim1==dim2)
                        H{dim1,dim2} = convn(raw,gaussD2{dim1},'same');
                    else
                        H{dim1,dim2} = convn(raw,gaussD{dim1},'same');
                        H{dim1,dim2} = convn(H{dim1,dim2},gaussD{dim2}, ...
                                             'same');
                    end
                    for dim = setdiff(1:3,[dim1,dim2])
                        H{dim1,dim2} = convn(H{dim1,dim2},gauss{dim}, ...
                                             'same');
                    end
                end
            end
            [nx,ny,nz] = size(raw);
            ev = SynEM.Aux.eig3S([H{1,1}(:)'; H{2,1}(:)'; H{3,1}(:)'; ...
                                  H{2,2}(:)'; H{3,2}(:)'; H{3,3}(:)']);
            
            % sort by absolute value
            [~,sortIds] = sort(abs(ev),1);
            linearIds = bsxfun(@plus, sortIds, ...
                (0:(size(sortIds,2) - 1)).*3);
            ev = ev(linearIds);
            fm = cell(3,1);
            fm{1} = reshape(ev(1,:),nx,ny,nz);
            fm{2} = reshape(ev(2,:),nx,ny,nz);
            fm{3} = reshape(ev(3,:),nx,ny,nz);
        end
    end

end
