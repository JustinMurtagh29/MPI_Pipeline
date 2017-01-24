function names = createNameStrings(obj)
%CREATENAMESSTRINGS Create a name string cell array encoding the exact
%nature of each feature.
% OUTPUT names: [Nx1] cell array of length obj.numFeatures. Each cell
%           contains the name of the respective feature (column) in the
%           feature output matrix. The form of the name for texture
%           features is
%               FeatureName_id_Vol_SumStat_Channel,
%           where Vol \in \{v0,...,v_n\}
%                 SumStats \in \{p_0,...p_n\}
%                 Channel \in \{c1_,...,c_n\}.
%           For shape feature the form is
%               FeatureName_Id_Channel
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>


%prepare texture name strings
pStr = arrayfun(@(x)['_p', int2str(x)], ...
    1:length(obj.quantiles) + sum(obj.moments),'UniformOutput',false)';
vStr = arrayfun(@(x)strcat(['_v' num2str(x)],{'s1','s2'}), ...
	1:obj.numSubvolumes,'UniformOutput',false);
vStr = [{'_v0'} [vStr{:}]];
cStr = permute(arrayfun(@(x)['_c', int2str(x)], ...
    1:max(cellfun(@(x)x.numChannels, obj.featTexture)), ...
    'UniformOutput',false),[1 3 2]);
tBase = strcat(repmat(vStr,length(pStr),1,length(cStr)), ...
    repmat(pStr,1,length(vStr),length(cStr)), ...
    repmat(cStr,length(pStr),length(vStr),1));

%create texture feature names
nT = length(obj.featTexture);
namesT = cellfun(@(x,y)strcat([x.name, '_', int2str(y)],...
    tBase(:,:,1:x.numChannels)), ...
    obj.featTexture,num2cell(1:nT)','UniformOutput',false);
namesT = cellfun(@(x)x(:), namesT, 'UniformOutput',false);
namesT = vertcat(namesT{:});

%shape featue names
namesS = cellfun(@(x,y)strcat([x.name, '_' int2str(y)], ...
    arrayfun(@(x)strcat('_c',int2str(x)),1:x.numFeatures, ...
    'UniformOutput', false)), ...
    obj.featShape, num2cell(1:length(obj.featShape))', ...
    'UniformOutput',false);
namesS = cellfun(@(x)x(:), namesS, 'UniformOutput',false);
namesS = vertcat(namesS{:});
obj.names = [namesT; namesS];
names = obj.names;
end