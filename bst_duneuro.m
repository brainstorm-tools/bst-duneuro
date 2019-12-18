function [Gain, errMsg] = bst_duneuro(OPTIONS)
% BST_DUNEURO: Call Duneuro to compute a FEM solution for Brainstorm.
%
% USAGE:  [Gain, errMsg] = bst_duneuro(OPTIONS)
%   todo?         DuneuroDir = bst_duneuro('update')
%   todo?         DuneuroDir = bst_duneuro('download')
%
% INPUT:
%     - OPTIONS: structure with the following fields
%        |- MEGMethod    : 'duneuro', else ignored
%        |- EEGMethod    : 'duneuro', else ignored
%        |- ECOGMethod   : 'duneuro', else ignored
%        |- SEEGMethod   : 'duneuro', else ignored
%        |- Channel      : Brainstorm channel structure
%        |- iMeg         : Indices of MEG sensors in the Channel structure
%        |- iEeg         : Indices of EEG sensors in the Channel structure
%        |- GridLoc      : Dipoles locations
%        |- FemHeadFile         : structure with the name of the fem head model
%        |- GridLoc      : Dipoles locations
%        |- conductivity      : The values of the conductivity %% < === conductivity and no Conductivity <== TODO : shoud to fixe it  
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

% tic;

% ===== DOWNLOAD DUNEURO  =====
% TODO :  TO DECIDE & integrate on bst disribution or download ?

% ===== Set number of cores  =====
% TODO : check whith Duneuro team if it's possibel

% ===== Set number of cores  =====
% TODO : check whith Duneuro team if it's possibel
bst_progress('setimage', 'logo-duneuro.png');
bst_progress('start', 'Head modeler', 'Initialization...');
% Intialize variables
Dims = 3;
nv = length(OPTIONS.GridLoc);
Gain = NaN * zeros(length(OPTIONS.Channel), Dims * nv);errMsg = '';
% Save current folder
curdir = pwd;

solverMethod = 'duneuro'; % TODO : to change to 'duneuro' later
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
str = which('bst_duneuro','-all'); % <== May be not needed if we include the main function
[filepath,~,~] = fileparts(str{1});                      % < == //                       //                   //                    //

% TODO : some of these parameters and other should be tuned from outside
bst_progress('text', 'Duneuro: prepare the input ...');
cfg = [];
cfg.pathOfTempOutPut = TmpDir;
cfg.pathOfDuneuroToolbox = filepath;
cfg.currentPath = curdir;
cfg.runFromBst  = 1; % update the bst_progress text
cfg.displayComment  = 0; % update the bst_progress text
cfg.femSourceModel = OPTIONS.FemSourceModel;% TODO : from outside :  partial_integration, venant, subtraction ==> check set_minifile
cfg.deleteOutputFolder = 0; % brainstorm will manage the rest
%%%
if isEeg; cfg.modality = 'eeg'; goodChannel = OPTIONS.iEeg; end
if isMeg; cfg.modality = 'meg'; goodChannel = OPTIONS.iMeg; end
if isEeg && isMeg; cfg.modality = 'meeg'; goodChannel = [OPTIONS.iMeg OPTIONS.iEeg]; end

if isEcog; cfg.modality = 'ecog'; goodChannel = OPTIONS.iEcog; end
if isSeeg; cfg.modality = 'seeg'; goodChannel = OPTIONS.iSeeg; end

%% ===== DUNEURO PREPARE GEOMETRY =====
%% 1 - Head Model
ProtocolInfo = bst_get('ProtocolInfo');
FemHeadFile = fullfile(ProtocolInfo.SUBJECTS,OPTIONS.FemHeadFile.FileName);
femhead = load(FemHeadFile);
% femhead : contains the bst FEM mesh format
cfg.node = femhead.Vertices;cfg.elem = [femhead.Elements femhead.Tissue];
cfg.tissuLabel = femhead.TissueLabels;
% clear femhead

% TODO : Get here the information about the tissu to keep (similar to the openmeeg panel that check the scalp, inner and outer skull )
% vector of one and zeros, the one mean keep the layer acording to the rank
% of the index on the vector. 
% [1 1 0 0 ] ==> there is four layer and here the model will keep the label
% 1 and 2 ===> ALWAYS the  label are from inner to outer  .... 
if isEeg; cfg.layerToKeep = ones(length(unique(femhead.Tissue))); end % < == situation of 1 1 1 1 1 1 <== keep all, not needed if all 1;
if isMeg; cfg.layerToKeep = OPTIONS.layerToKeep; end 

%% 2- Source space
cfg.sourceSpace = OPTIONS.GridLoc;

% %%% TODO adapt it for duneuro fem
% %Center all the points on the center of the envelope
%     bfs_center = bst_bfs(TessMat.Vertices);
%     vDipoles = bst_bsxfun(@minus, OPTIONS.GridLoc, bfs_center(:)');
%     vInner = bst_bsxfun(@minus, TessMat.Vertices, bfs_center(:)');
% %Check if any dipole is outside of the innermost layer%     % Check if any dipole is outside of the innermost layer
%     iDipInside = find(~inpolyhd(vDipoles, vInner, TessMat.Faces));
%     if ~isempty(iDipInside)
%         errMsg = sprintf(['Some dipoles are outside the BEM layers (%d dipoles).\n' ...
%                           'The leadfield for these dipoles is probably incorrect.\n\n'], length(iDipInside));
%         if strcmpi(OPTIONS.HeadModelType, 'surface')
%             errMsg = [errMsg, 'To fix the cortex surface:', 10, 'Right-click on the surface file > Force inside skull.'];
%         end
%         disp([10 'WARNING: ' errMsg 10]);
%     end
%  %%%%   

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
    clear channelLoc  eegloc;
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
    clear MegChannel;
end
%% 4- Conductivity/tensor
% TODO : Adapt this for anisotrpy // maybe could be done from outside and
% put the tensor on OPTIONS.Conductivity or create new field OPTIONS.ConductivityTensor
cfg.conductivity = OPTIONS.conductivity;

%% 5- Duneuro-Brainstorm interface
%% ===== DUNEURO LEADFIELD COMPUTATION   =====
bst_progress('text', 'Duneuro: computing ...');
cfg = bst_duneuro_interface(cfg); 
% displayOutput = 1; % 1 yes, 0 no
 % if displayOutput == 1; tic; cfg = bst_duneuro_interface(cfg); t_direct = toc; end
% if displayOutput == 0; tic; evalc('cfg = bst_duneuro_interface(cfg);'); t_evalc = toc; end % could be longer ... not really ==> todo check for high resolution

%% 6- Read the lead field
% fil the bad channel with nan (not possible within the duneuro)
if isMeg;     Gain(OPTIONS.iMeg, :) = cfg.fem_meg_lf; end % the final value ==>  cfg.fem_meg_lf = B = Bp + Bs ? check if the case for bst
if isEeg;     Gain(OPTIONS.iEeg, :) = cfg.fem_eeg_lf; end 
% if isEeg; Gain(goodChannel,:) = cfg.fem_eeg_lf; end

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

% t_fem = toc
end
