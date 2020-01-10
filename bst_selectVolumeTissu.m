function varargout = bst_selectVolumeTissu(varargin)
% varargout = bst_selectVolumeTissu(varargin)
% cfg = bst_selectVolumeTissu(cfg)
% cfg = bst_selectVolumeTissu(cfg,layerToKeep)
% [node, elem] = bst_selectVolumeTissu(node,elem,layerToKeep)
% Function that keep the desired tissus index specified on the booleqn vector layerToKeep or cfg.layerToKeep.
% Important ... the tissu have to be labled in ascendant manner, from the
% most inner layer to the most outer layer.
% example : layerToKeep = [0 1 0 1]; will keep the 2nd and the fourth layer.

% Takfarinas MEDANI

%input : vector with boolean index for each tissu

if isstruct(varargin{1})
    cfg = varargin{1};
    if length(varargin) == 2
        layerToKeep =  varargin{2};
    elseif isfield(cfg,'layerToKeep')
        layerToKeep = cfg.layerToKeep;
    elseif ~isfield(cfg,'layerToKeep')
        error('The layer(s) to keep are not specified on cfg.layerToKeep');
    end
    func = 1;
else
    if length(varargin) < 3
        error('Please specify the three input argunents  : [node, newelem] = bst_selectVolumeTissu(node,elem,layerToKeep)' )
    end
    node = varargin{1};
    elem = varargin{2};
    layerToKeep = varargin{3};
    func = 2;
end

if func ==1
    tissuLabel = unique(cfg.elem(:,end));
    indexToKeep = (zeros(length(cfg.elem),1));
    numberOfLayer = length(tissuLabel);
    if sum(layerToKeep) == numberOfLayer
        varargout{1} = cfg;
        disp('All the tissus are selected ... don''t need to run this function')
        return;
    end
    
    if length(layerToKeep) ~= numberOfLayer
        error('Error on the head model, the number of layer available is diffrent from the specified  in the input');
    end
    
    % Find the index of the layers to keep
    labelToKeep = tissuLabel((boolean(layerToKeep(:))));
    for ind = 1 : length(labelToKeep)
        indexToKeep = ((cfg.elem(:,end) == labelToKeep(ind)) | (indexToKeep));
    end
    
    % Extract the element to keep
    newelem = cfg.elem(indexToKeep,:);
    [newelem, ~ ] = meshreorient(cfg.node,newelem(:,1:end-1));
    
    % Updates with the new list of eleme
    cfg.old_elem = cfg.elem; %% NOT GOOD IDEA TO STORE BOTH DATA... MEMORY CONSUMPTION
    cfg.elem = [newelem  cfg.elem(indexToKeep,5)];
%     figure; plotmesh(cfg.node, cfg.elem,'x>0')
    varargout{1} = cfg;
else % fun == 2
    tissuLabel = unique(elem(:,end));
    indexToKeep = (zeros(length(elem),1));
    numberOfLayer = length(tissuLabel);
    if sum(layerToKeep) == numberOfLayer;
        varargout{1} = node;
        varargout{2} = elem;
        disp('All the tissus are selected ... don''t need to run this function')
        return;
    end
    if length(layerToKeep) ~= numberOfLayer
        error('Error on the head model, the number of layer available is diffrent from the specified  in the input');
    end
    % Find the index of the layers to keep
    labelToKeep = tissuLabel((boolean(layerToKeep(:))));
    for ind = 1 : length(labelToKeep)
        indexToKeep = ((elem(:,end) == labelToKeep(ind)) | (indexToKeep));
    end
    % Extract the element to keep
    newelem = elem(indexToKeep,:);
    [newelem, ~ ] = meshreorient(node,newelem(:,1:end-1));
    % Updates with the new list of eleme
    %     old_elem = elem; %% NOT GOOD IDEA TO STORE BOTH DATA... MEMORY CONSUMPTION
    elem = [newelem  elem(indexToKeep,5)];    clear newelem
%     figure; plotmesh(node,elem,'x>0');
    varargout{1} = node;
    varargout{2} = elem;
end
%%%% DO NOT REMOVE THE ISOLATED NODE ... THIS WILL AFFECT THE NODE INDEX,
%%%% THEN WILL DISTTURB THE MESH ELEMENT AND ALSO THE CONDUCTIVITY LABEL
end
