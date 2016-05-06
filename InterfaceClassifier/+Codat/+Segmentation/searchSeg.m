function [ idx, name, params ] = searchSeg( result, mode, varargin )
%SEARCHSEG Search a segmentation from a result structure. Search will always be
% done for the training segmentation.
% INPUT result: See Codat.Segmentation.parameterSearchSeg output.
%       mode: String specifying the segmentation mode
%           'merge': Get segmentation with average merge free path length
%               closest to the specified one.
%               varargin:
%               	1: Double specifying the desired merge rate in um.
%                   2: Integer specifing the index of the sm curve
%                   (corresponding to the index of the desired node
%                   threshold.) (Default: 1)
%           'split': Get segmentation with average split free path length
%               closest to the specified one.
%               varargin:
%                   1: Double specifying the desired split rate in um.
%                   2: Integer specifing the index of the sm curve
%                   (corresponding to the index of the desired node
%                   threshold.) (Default: 1)
%       varargin: See mode.
% OUTPUT idx: Integer index of the corresponding smCurve. For the modes
%           'merge' and 'split idx are the indices of all segmentation
%           sorted by increasing distance to the specified target rates.
%        name: String containing the name of the segmentation file.
%        params: Cell array of the corresponding parameters.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

switch mode
    case 'merge'
        if length(varargin) < 2
            varargin{2} = 1;
        end
        ms = result.smCurveTrain{varargin{2}}{:,:};
        [~,idx] = sort(abs(ms(:,1) - varargin{1}),'ascend');
        name = result.paramsTrain{idx,3};
        params = result.paramsTrain{idx,2};
    case 'split'
        if length(varargin) < 2
            varargin{2} = 1;
        end
        ms = result.smCurveTrain{varargin{2}}{:,:};
        [~,idx] = sort(abs(ms(:,2) - varargin{1}),'ascend');
        name = result.paramsTrain{idx,3};
        params = result.paramsTrain{idx,2};
    otherwise
        error('Mode %s not implemented.',mode)
end


end

