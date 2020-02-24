function cfg = bst_prepare_head_model(cfg)
% cfg = bst_prepare_head_model(cfg)
% write the head model file to the current path
% The format can be either ".msh" or ".geo" according
% to the simulation.
% cfg.useTensor = 0;  In this case the conductivity is
%                                 be a scalar, then ".msh" file is used.
% cfg.useTensor = 1;  In this case the conductivity is
%                                 be a tensor, then ".geo" file is used.
%
% File created by Takfarinas MEDANI, November 2019;


% TODO : Optimisation ... use binary input/output
%               Use a standard format for all modalities
%               Discuss with duneuro team to add the I/O files

%% Set the minimal parameters
if ~isfield(cfg,'filename'); cfg.filename = 'head_model'; end  % Use the default name.
if ~isfield(cfg,'useTensor'); cfg.useTensor = 0; end  % Use the tensor model or not
if ~isfield(cfg,'isotrop'); cfg.isotrop = 1; end  %
if ~isfield(cfg,'modality'); cfg.moadlity = 'eeg'; end  %
if ~isfield(cfg,'layerToKeep'); cfg.layerToKeep = ones(length(unique(cfg.elem(:,5))),1);  end  %
if ~isfield(cfg,'savefile'); cfg.savefile = 1; end
if ~isfield(cfg,'saveBstFormat'); cfg.saveBstFormat = 0; end
%% Check if the wole head model will be used in the case of the MEG

if strcmp(cfg.modality,'meg')
    % Apply the reduced volume for the MEG, this function will remove the
    % unselected tissu from the head model, this could be also done for
    % sEEG later
    % 
    cfg = bst_selectVolumeTissu(cfg,cfg.layerToKeep);
end

if strcmp(cfg.modality,'meeg')
    if sum(cfg.layerToKeep) ~= length(unique(cfg.elem(:,5)))
        warning(['THE REDUCED MEG HEAD MODEL WILL BE ALSO USED BY THE EEG'...
                   ' >> NOT CORRECT : SEPARATE HEAD MODEL SHOULD BE USED >> NOT IMPLEMENTED YET'])
    end
%% @@@ WARNING @@@
% For the combined MEEG .... if the user select specific layers for the MEG
% the new model will be also used by the EEG .... and this is not correct
% FOR the MEEG and cfg.layerToKeep contient des zeros .... == > add
% condition that generate two head model and run duneuro for each
% modalities. ==> Open an issue ... this could be also a proble for
% openmeeg 
end

%% Check the conductivity
if cfg.isotrop == 1    
    if cfg.useTensor == 0
        if strcmp(cfg.dnMeshElementType,'hexahedron')
            cfg.head_filename = [cfg.filename '.dgf'];
            cfg.saveDfgFormat =1;
        else
            % Duneuro uses the msh file
            cfg.head_filename = [cfg.filename '.msh'];
            cfg.saveMshFormat = 1;
        end
    else % case cfg.useTensor =1        
        if strcmp(cfg.dnMeshElementType,'hexahedron')
            error('Using the tensor Model with hexahedron mesh is not supported for now')
        end
        % Duneuro uses the Cauchy files
        cfg.head_filename = [cfg.filename '.geo'];
        cfg.saveCauFormat = 1;
    end
    
else % otherwise it anisotrop, and for  sure the geo is used
    if strcmp(cfg.dnMeshElementType,'hexahedron')
        error('Using the anisotropy model with hexahedron mesh is not supported for now');
        % TODO : Check if the geo file could be used for the hexa
    end
    % Duneuro uses the Cauchy files
    cfg.head_filename = [cfg.filename '.geo'];
    cfg.saveCauFormat = 1;
    % Just to make sure for nestx steps
    cfg.useTensor = 1;
end

%% Write the mesh model
bst_write_mesh_model(cfg);
end