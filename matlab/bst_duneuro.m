function [Gain, errMsg] = bst_duneuro(OPTIONS)
% BST_DUNEURO: Call Duneuro to compute a FEM solution for Brainstorm.
%
% USAGE:  [Gain, errMsg] = bst_duneuro(OPTIONS)
%  Todo?         DuneuroDir = bst_duneuro('update')
%   Todo?         DuneuroDir = bst_duneuro('download')
%
% INPUT: TODO
%     - OPTIONS: Todo :structure with the following fields
% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c)2000-2019 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Takfarinas MEDANI, December 2019;

% ===== DOWNLOAD DUNEURO  =====
if exist('bst_duneuro', 'file')
    disp([10, 'bst-duneuro path: ', bst_fileparts(which('bst_duneuro')), 10]);
else
    error('bst-duneuro toolnbox is not installed')
end
% ===== Set number of cores  =====
% TODO : check whith Duneuro team if it's possibel

% ===== Start duneuro process  =====
bst_progress('start', 'Head modeler', 'Initialization...');
% Intialize variables
Dims = 3;
nv = length(OPTIONS.GridLoc);
Gain = NaN * zeros(length(OPTIONS.Channel), Dims * nv);errMsg = '';
% Save current folder
curdir = pwd;
solverMethod = 'duneuro';
isEeg  = strcmpi(OPTIONS.EEGMethod, solverMethod )  && ~isempty(OPTIONS.iEeg);
isMeg  = strcmpi(OPTIONS.MEGMethod, solverMethod)  && ~isempty(OPTIONS.iMeg);
isEcog = strcmpi(OPTIONS.ECOGMethod, solverMethod) && ~isempty(OPTIONS.iEcog);
isSeeg = strcmpi(OPTIONS.SEEGMethod, solverMethod) && ~isempty(OPTIONS.iSeeg);
% Get temp folder
TmpDir = bst_get('BrainstormTmpDir');
% Open log file
logFile = bst_fullfile(TmpDir, 'duneuro_log.txt');
fid_log = fopen(logFile, 'w');

%% ===== PREPARE THE INTERFACE TO DUNEURO =====
%% 0 - Initialisation for duneuro computation (configuration used by the duneuro interface)
%  find the bst_duneuro_toolbox path
str = which('bst_unique_readme.txt','-all'); % <== later use : bst_fullfile(bst_get('BrainstormUserDir'), 'bst_duneuro');
filepath = fileparts(str{1});                      % < == //                       //                   //                    //
bst_progress('text', 'Duneuro: prepare the input ...');
bst_progress('setimage', fullfile(filepath,'matlab','external', 'logo-duneuro.png'));
% Get current subject
sStudy = bst_get('Study');
[sSubject, iSubject] = bst_get('Subject', sStudy.BrainStormSubject);
if isempty(sSubject) || isempty(sSubject.iCortex)
    errMessage = 'No cortex surface available for this subject.';
end
if isempty(sSubject) || isempty(sSubject.iFEM)
    errMessage = 'No FEM head model available for this subject.';
end
% Path to the head fem model
OPTIONS.FemFiles = file_fullpath(sSubject.Surface(sSubject.iFEM).FileName);
% Get the information from the head model
fields = whos('-file',OPTIONS.FemFiles);
% Get number of  layers
ivar = find(strcmpi({fields.name}, 'TissueLabels'));
numberOfLayer =  max(fields(ivar).size);
% Get mesh element type
ivar = find(strcmpi({fields.name}, 'Elements'));
numberOfEdges =  (fields(ivar).size(2));
if numberOfEdges == 4
    dnMeshElementType = 'tetrahedron';
elseif numberOfEdges == 8
    dnMeshElementType = 'hexahedron';
else
    error('Mesh element type is not defined')
end
% Get the mesh method generation
load(OPTIONS.FemFiles, 'Comment');
if contains(Comment,'iso2mesh')
    meshMethod = 'iso2mesh';
elseif contains(Comment,'SimNibs')
    meshMethod = 'SimNibs';
else
    meshMethod = 'others';
end
% Get the default conductivity values
OPTIONS.Conductivity  = get_standard_conductivity((numberOfLayer));
% Get the names of the tissues
load(OPTIONS.FemFiles, 'TissueLabels');
OPTIONS.FemNames = TissueLabels;

%% ===== DUNEURO =====

cfg = []; % the most important variable that will pass all the variable to prepare the interface to duneuro
cfg.runFromBst  = 1; % update the bst_progress text
cfg.currentPath = curdir;
cfg.pathOfTempOutPut = TmpDir;
cfg.pathOfDuneuroToolbox = filepath; % USE : bst_fullfile(bst_get('BrainstormUserDir'), 'bst_duneuro');
cfg.displayComment  = 0; % update the bst_progress text
cfg.dnMeshElementType = dnMeshElementType;
cfg.deleteOutputFolder = 0; % brainstorm will manage the rest
cfg.meshMethod = meshMethod;
%% Add THE BIG PANEL HERE and Move the previous panel to HERE
% DuneuroMeegOptions = gui_show_dialog('Duneuro Options', @panel_duneuro, 1, [], OPTIONS);
% for now we are going to use this following part
isDuneuro = 1;
if isDuneuro
    %% FEM solver parametrs advanced:
    [res, isCancel] =  java_dialog('question', '<HTML><B> DUNEuro Parameters <B>', ...
        'Set Duneuro Parameters', [],{'Use Default Options (recommended)', ...
        'Set Advanced Options (Expert Mode)' }, ...
        'Use Default Options (recommended)');
    if isCancel ==1
        return;
    end
    if strcmpi(res,'Use Default Options (recommended)')
        cfg.advancedMode = 0;
    else
        cfg.advancedMode = 1;
    end
    %% Advanced users
    if cfg.advancedMode ==1
        %% Ask for the source FEM Method
        dnFemMethodType = {'fitted', 'unfitted'};
        [res, isCancel] = java_dialog('question', '<HTML><B> DUNEuro : Select FEM Method Type (advanced) <B>', ...
            'FEM Source Model', [], dnFemMethodType, 'fitted');
        if isCancel
            return
        end
        cfg.dnFemMethodType = res;
        %% Ask for the source FEM soler Type
        dnFemSolverType = {'CG : Continious Galerkin','DG : Discontinious Galerkin'}; % Continious Galerkin and Discountinious Galerkin : Default cg,
        [res, isCancel] =  java_dialog('question', '<HTML><B> DUNEuro : Select FEM Solver Type (advanced) <B>', ...
            'FEM Source Model', [], dnFemSolverType, 'fitted');
        if isCancel
            return
        end
        cfg.dnFemSolverType = lower(res(1:2));
        % Get the meshing type from the inputs
        if strcmp(cfg.dnMeshElementType,'hexahedron')
            cfg.dnGeometryAdapted = 'true'; % 'true, or false' ==> should be tuned ffom the panel
            dnGeometryAdapted = {'true', 'false'};
            res = java_dialog('question', '<HTML><B> DUNEuro : Use dnGeometryAdapted <B>', ...
                'FEM Geometry', [], dnGeometryAdapted, 'true');
            cfg.dnGeometryAdapted = res;
        end
        %% Ask for the source FEM Model
        femSourceModel =  {'Venant','Subtraction','Partial_Integration'};
        [res, isCancel] =  java_dialog('question', '<HTML><B> DUNEuro : Select FEM Source Model <B>', ...
            'FEM Source Model', [],femSourceModel, 'Venant');
        if isCancel
            return
        end
        cfg.femSourceModel = res;
        % In the case of the Venant source model
        if strcmpi((cfg.femSourceModel),'venant')
            % Ask Venant options
            [res, isCancel] =  java_dialog('input', ...
                {'Number Of Moments (1-5):', ...
                'Reference Length (1-100):', ...
                'Weighting Exponent (1-5):', ...
                'Relaxation Factor (1-5):', ...
                'Mixed Moments (true-false):', ...
                'Restricted (true-false):'}, ...
                'Venant Options', [], {'3', '20','1','1e-6','true','true'});
            if isCancel
                return
            end
            cfg.femSourceModelNumberOfMoments =  str2num(res{1});
            cfg.femSourceModelReferenceLength =  str2num(res{2});
            cfg.femSourceModelWeightingExponent =  str2num(res{3});
            cfg.femSourceModelRelaxationFactor =  str2num(res{4});
            cfg.femSourceModelMixedMoments = res{5};
            cfg.femSourceModelRestrict = res{6};
        elseif strcmpi((cfg.femSourceModel),'subtraction')
            % todo
        elseif strcmpi((cfg.femSourceModel),'partial_integration')
            % todo
        elseif strcmpi((cfg.femSourceModel),'whitney')
            % todo
        elseif strcmpi((cfg.femSourceModel),'other')
            % todo
        else
            error('FEM source model is not defined');
        end
    else % regular user
        disp('All the DUNEuro options will be set to the default values')
    end
    %% Ask for the conductivity of each layer
    for ind = 1 : length(OPTIONS.Conductivity)
        cond{ind} = num2str(OPTIONS.Conductivity(ind));
    end
    [res, isCancel] = java_dialog('input', TissueLabels, 'Conductivity Values', [], cond);
    if isCancel
        return
    end
    for ind = 1 : length(OPTIONS.Conductivity)
        OPTIONS.FemConductivity(ind) = str2num(res{ind});
    end
    % TODO : The conductivities values in the case of the
    % combined model should be the same for eeg and meg
    %% Ask for the layers to keep for the MEG computation
    % We can not use different volme conductor in the case of combined
    % eeg/meg ... otherwise we need to call separately the binaries ....
    % TODO : discuss this point with th the team.
    if strcmpi(OPTIONS.EEGMethod, 'duneuro') && strcmpi(OPTIONS.MEGMethod, 'duneuro')
        [res, isCancel] = java_dialog('checkbox', ...
            '<HTML>Select the layers to consider for the combined MEG/EEG head modeling <BR>', 'Select Volume', [], ...
            TissueLabels, [ones(1, length(TissueLabels))]);
        if isCancel
            return
        end
        layerToKeep =  res;
    elseif  strcmpi(OPTIONS.EEGMethod, 'duneuro')
        [res, isCancel] = java_dialog('checkbox', ...
            '<HTML>Select the layers to consider for the EEG head modeling <BR>', 'Select Volume', [], ...
            TissueLabels, [ones(1, length(TissueLabels))]);
        layerToKeep =  res;
        if isCancel
            return
        end
    elseif strcmpi(OPTIONS.MEGMethod, 'duneuro')
        [res, isCancel] = java_dialog('checkbox', ...
            '<HTML>Select the layers to consider for the MEG head modeling <BR>', 'Select Volume', [], ...
            TissueLabels, [1 zeros(1, length(TissueLabels)-1)]);
        if isCancel
            return
        end
        layerToKeep =  res;
    else
        % TODO : add similar option in the case of iEEG and sEEG
    end
end

% Get the modality from the OPTIONS structure
if isEeg; cfg.modality = 'eeg'; goodChannel = OPTIONS.iEeg; end
if isMeg; cfg.modality = 'meg'; goodChannel = OPTIONS.iMeg; end
if isEeg && isMeg; cfg.modality = 'meeg'; goodChannel = [OPTIONS.iMeg OPTIONS.iEeg]; end
if isEcog; cfg.modality = 'ecog'; goodChannel = OPTIONS.iEcog; end
if isSeeg; cfg.modality = 'seeg'; goodChannel = OPTIONS.iSeeg; end

%% ===== DUNEURO PREPARE GEOMETRY =====
%% 1 - Head Model
femhead = load(OPTIONS.FemFiles);
% femhead : contains the bst FEM mesh format
cfg.node = femhead.Vertices;
cfg.elem = [femhead.Elements femhead.Tissue];
cfg.tissuLabel = femhead.TissueLabels;
% Get here the information about the tissu to keep (similar to the openmeeg panel that check the scalp, inner and outer skull )
% vector of one and zeros, the one mean keep the layer acording to the rank
% of the index on the vector.
% [1 1 0 0 ] ==> there is four layer and here the model will keep the label
% 1 and 2 ===> ALWAYS the  label are from inner to outer  ....
% if isEeg; cfg.eegLayerToKeep = eegLayerToKeep; end % < == situation of 1 1 1 1 1 1 <== keep all, not needed if all 1;
% if isMeg; cfg.megLayerToKeep = megLayerToKeep; end
% if isMeg; cfg.megLayerToKeep = megLayerToKeep; end
% if isEcog; cfg.eCogLayerToKeep =  eCogLayerToKeep; end
% if isSeeg; cfg.seegLayerToKeep =  seegLayerToKeep; end
% temporary option for testing the rest of the code
cfg.layerToKeep = layerToKeep; % or may be this one should be sufficient

%% 2- Source space
% only in the case of SimNibs or equivalent : TODO : check the input
% get this parameters from the user panel
if ~strcmp(cfg.meshMethod, 'iso2mesh')
    cfg.shrinkSourceSpace = 1;
    cfg.sourceSpaceDepth =1.5;
    if cfg.advancedMode ==1
        res = java_dialog('input', 'Source Space depth from Pial Surface (mm)', 'FEM Source Space', [], num2str(cfg.sourceSpaceDepth));
        cfg.sourceSpaceDepth = str2num(res);
        if isCancel ==1
            return;
        end
    end
else
    cfg.shrinkSourceSpace = 0;
end

% Source space type
switch (OPTIONS.HeadModelType)
    case 'volume'
        % TODO or keep it as it's now....
    case 'surface'
        % Read cortex file
        sCortex = bst_memory('LoadSurface', OPTIONS.CortexFile);
        if cfg.shrinkSourceSpace
            % Shrink the cortex surface by XX mm
            nrm_sur_nodes  =  sCortex.VertNormals;
            % I- Trouver les composante spherique de chaque normale à la surface:
            [azimuth, elevation]  =  cart2sph(nrm_sur_nodes(:,1),nrm_sur_nodes(:,2),nrm_sur_nodes(:,3));
            %% II - Choisir la profondeur de la postion de l'espace des sources dans la couche du cortex :
            profondeur_sources = cfg.sourceSpaceDepth/1000; % en mm
            % III- Trouver les composantes dans les trois directions à partir du point centroide de chaque facette:
            profondeur_x  =  profondeur_sources .* cos(elevation) .* cos(azimuth);
            profondeur_y  =  profondeur_sources .* cos(elevation) .* sin(azimuth);
            profondeur_z  =  profondeur_sources .* sin(elevation);
            % verif = sqrt(profondeur_x.^2+profondeur_y.^2+profondeur_z.^2);
            % IV- Appliquer cette profondeur dans chque direction à partir du centroide de chaque feacette:
            pos_source_x = sCortex.Vertices(:,1) - profondeur_x;
            pos_source_y = sCortex.Vertices(:,2) - profondeur_y;
            pos_source_z = sCortex.Vertices(:,3) - profondeur_z;
            %             figure;
            %             plotmesh([pos_source_x pos_source_y pos_source_z],sCortex.Faces,'edgecolor','none','facecolor','r');hold on
            %             plotmesh(sCortex.Vertices,sCortex.Faces,'edgecolor','none','facecolor','k');hold on
            %             % check if the nodes are on the GM volume
            %             iVertOut = find(~inpolyhd([pos_source_x pos_source_y pos_source_z], ...
            %                                                                 sCortex.Vertices, sCortex.Faces));
            % Surface: Use the cortex surface
            OPTIONS.GridLoc    = [pos_source_x pos_source_y pos_source_z];
        else
            OPTIONS.GridLoc =  sCortex.Vertices;
        end
    case 'mixed'
        %TODO : not used ?
end
cfg.sourceSpace = OPTIONS.GridLoc;

%% 3- Electrode Position / Coils positon / Orientation
% === EEG ===
if isEeg
    %     EegChannel = [];
    %     for iChan = 1: length(OPTIONS.iEeg)
    %         sChan = OPTIONS.Channel(OPTIONS.iEeg(iChan));
    %         if ~isempty(OPTIONS.Channel(iChan).Loc) %
    %             EegChannel(iChan,:) = OPTIONS.Channel(iChan).Loc;
    %         else
    %             EegChannel(iChan,:) = [];
    %         end
    %     end
    % Channel location :
    eegloc = cat(2, OPTIONS.Channel(OPTIONS.iEeg).Loc)';
    cfg.channelLoc = eegloc;
end
% === MEG ===
% get here the infos about the sensors
if isMeg
    MegChannel = [];
    for iChan = 1 : length(OPTIONS.iMeg)
        sChan = OPTIONS.Channel(OPTIONS.iMeg(iChan));
        for iInteg = 1:size(sChan.Loc, 2)
            MegChannel =  [MegChannel; iChan, sChan.Loc(:,iInteg)', sChan.Orient(:,iInteg)', sChan.Weight(iInteg)];
        end
    end
    cfg.MegChannel = MegChannel;
end
%% 4- Conductivity/tensor
% TODO : Adapt this for anisotrpy // maybe could be done from outside and
% put the tensor on OPTIONS.Conductivity or create new field OPTIONS.ConductivityTensor
cfg.conductivity = OPTIONS.FemConductivity; % TODO : we can do it from here instead to modify the bst file

%% 5- Duneuro-Brainstorm interface
%% ===== DUNEURO LEADFIELD COMPUTATION   =====
bst_progress('text', 'Duneuro: computing ...');
cfg = bst_duneuro_interface(cfg);
% displayOutput = 1; % 1 yes, 0 no
% if displayOutput == 1; tic; cfg = bst_duneuro_interface(cfg); t_direct = toc; end
% if displayOutput == 0; tic; evalc('cfg = bst_duneuro_interface(cfg);'); t_evalc = toc; end % could be longer ... not really ==> todo check for high resolution

%% 6- Read the lead field
% fil the bad channel with nan (not possible within the duneuro)
if isMeg;     Gain(OPTIONS.iMeg, :) = cfg.fem_meg_lf; end % the final value ==>  cfg.fem_meg_lf = B = Bp + aBs ? check if the case for bst
if isEeg;     Gain(OPTIONS.iEeg, :) = cfg.fem_eeg_lf; end

%% ===== CLEANUP =====
% TODO : adapt it to bst
bst_progress('text', 'Duneuro: Emptying temporary folder...');
% TODO : Delete intermediary files
% /!\ Warning : the file "transferOut.dat" : should be stored somewhere
% it contains the transfer matrix, this matrix could be re-used if the users change,
% ONLY and ONLY, the source space and/or the source model.
%  If the head geometry (mesh), the electrode position and or the conductivity value are changed,
% the forwar model need to be recomputed again (de not reuse the "transferOut.dat")

% Close log file
if ~isempty(fid_log) && (fid_log >= 0) && ~isempty(fopen(fid_log))
    fclose(fid_log);
end
% Go back to initial folder
cd(curdir);

% Remove Duneuro image
bst_progress('removeimage');
bst_progress('stop');
end
