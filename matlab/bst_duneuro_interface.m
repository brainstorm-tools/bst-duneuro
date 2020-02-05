function cfg = bst_duneuro_interface(cfg)
% BST_DUNEURO_INTERFACE : Writes the arguments from bst to duneuro and run the FEM
%
% USAGE:      cfg = bst_duneuro_interface(cfg)
%
% INPUT:
%     - cfg: structure with the fields:
%              Run the example\demo  in order to have the full liste 
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


%% 0 - Initialisation 
% Set the default parameters used in this function
if ~isfield(cfg,'modality'); cfg.modality  = 'eeg'; end
if ~isfield(cfg,'runFromBst'); cfg.runFromBst = 0; end                    % Works only if called from brainstorm 
if ~isfield(cfg,'currentPath'); cfg.currentPath = pwd; end                 % This function will cd to a temporary file and then return here (pwd) 
if ~isfield(cfg,'useTransferMatrix'); cfg.useTransferMatrix = 1; end % use the transfer matrix is recommended (choice 0 only for duneuroVersion 1)
if ~isfield(cfg,'BstDuneuroVersion'); cfg.BstDuneuroVersion = 2; end                 % 1 previous with separate files, 2the new version combined eeg and meg and binary + txt output,   
if ~isfield(cfg,'isotrop'); cfg.isotrop =1; end                                       % important to specify in order to write the correct file format (1 will use MSH, 0 will use Cauchy)
if ~isfield(cfg,'lfAvrgRef'); cfg.lfAvrgRef = 0; end                              %  compute average reference 1, otherwise the electrode 1 is the reference and set to 0
                                                                                                          %   in that case the for duneuro the electrod 1 is the reference and set to 0                              
if ~isfield(cfg,'displayComment'); cfg.displayComment  = 0; end
if ~isfield(cfg,'pathOfTempOutPut'); cfg.pathOfTempOutPut =  cfg.currentPath; end  

if cfg.runFromBst == 1
    cfg.lfAvrgRef = 0;
    cfg.brainstormModality = cfg.modality; %
    if ~isfield(cfg,'brainstormOutputFolder') % ==> the output folder could be different from the temporary folder in the case where we want to save the transfer for example
        cfg.brainstormOutputFolder = cfg.pathOfTempOutPut;
    end % should be done from the bst
    %%%%% UPDATES these values from here :
    % subpart  [brainstorm]
    % if ~isfield(cfg,'brainstormEegSaveTransfer'); cfg.brainstormEegSaveTransfer = 'false'; end %
    % if ~isfield(cfg,'brainstormMegSaveTransfer'); cfg.brainstormMegSaveTransfer = 'false'; end %
end
                                                                            
                                                                   
cfg.displayComment = 0;
cfg.BstDuneuroVersion = 3;


%% ------------- DUNEURO INTERFACE ------------- %%

% cd(fullfile(cfg.pathOfTempOutPut));
%% 1- Prepare the head model : 
% Write the head file according to the configuration cfg
if cfg.runFromBst == 1; bst_progress('text', 'Duneuro (1/7): 1- Prepare the head model ... '); end
if cfg.displayComment ==1;disp(['duneruo >>1 - Writing the head file  to ' ((cfg.pathOfTempOutPut))]);end
if ~isfield(cfg,'filename');    cfg.filename = fullfile(cfg.pathOfTempOutPut,'head_model'); end  % Use the default name.
cfg = bst_prepare_head_model(cfg);

%% 2- Prepare the Source Model : Same for EEG and MEG 
% Write the source/dipole file
if cfg.runFromBst == 1; bst_progress('text', 'Duneuro (2/7): 2- Prepare the Source Model ... '); end
if cfg.displayComment ==1;disp(['duneruo >>2 - Writing the dipole file  to ' ((cfg.pathOfTempOutPut))]);end
cfg.dipole_filename  = fullfile(cfg.pathOfTempOutPut,'dipole_model.txt');
write_duneuro_dipole_file(cfg.sourceSpace,cfg.dipole_filename);

%% 3- Prepare the Sensor Model : 
if cfg.runFromBst == 1; bst_progress('text','Duneuro (2/7): 3- Prepare the Sensor Model ... '); end
if cfg.displayComment ==1;disp(['duneruo >>3 - Writing the sensor file  to ' ((cfg.pathOfTempOutPut))]);end
cfg = bst_prepare_sensor_model(cfg);

%% 4- Prepare the Conductivity/tensor Model
if cfg.runFromBst == 1; bst_progress('text', 'Duneuro (4/7): 4- Prepare the Conductivity/tensor Model ... '); end
if cfg.displayComment ==1;disp(['duneruo >>4 - Writing the conductivity/tensor file  to ' ((cfg.pathOfTempOutPut))]);end
cfg = bst_prepare_conductivity_model(cfg);

%%  5- Prepare the Duneuro Configuration file / the minifile
if cfg.runFromBst == 1; bst_progress('text',  'Duneuro (5/7): 5- Prepare the Duneuro Configuration file  ... '); end
if cfg.displayComment ==1;disp(['duneruo >>5 - Writing the duneuro configuration file  to ' ((cfg.pathOfTempOutPut))]);end
cfg = bst_prepare_minifile(cfg);

%% 6- Prepapre & Run the Duneuro Application
% cd(fullfile(cfg.pathOfDuneuroToolbox,'bin'));

if cfg.runFromBst == 1; bst_progress('text', 'Duneuro (6/7):  6- Prepapre & Run the Duneuro Application ... '); end
% define the command line
cfg = bst_set_duneuro_cmd(cfg);
%%%% @@ Run Duneuro @@ %%%%%%
tic; cfg = bst_run_duneuro_cmd(cfg); cfg.time_fem = toc;
if cfg.displayComment ==1;disp(['duneruo >>6 - FEM computation :  ' num2str(cfg.time_fem) ' s']);end

%% 7- Read the lead field matrix (EEG or/and MEG)
if cfg.runFromBst == 1; bst_progress('text', 'Duneuro (7/7):  7- Post-process the LeadField ... '); end

if cfg.displayComment ==1;disp(['duneuro >>7 - Read the leadfield from ' (fullfile(cfg.pathOfTempOutPut))]);end
cfg = bst_read_duneuro_leadfield(cfg);

%% 8- Post Process the leadfield
if cfg.displayComment ==1;disp(['duneuro >>8 - Post-process the leadfield ' ]);end
cfg = bst_postprocess_lf(cfg);

%% 9- Remove the temporary folder
if cfg.displayComment ==1;disp(['duneruo >>9 - Clean the folder  : ' (fullfile(cfg.pathOfTempOutPut))]);end
if ~isfield(cfg,'deleteOutputFolder'); cfg.deleteOutputFolder = 0; end
if cfg.deleteOutputFolder == 1
    % TODO or leave bst to do it
    disp(['remove the '  (fullfile(cfg.pathOfTempOutPut)) ' from the hard disc']);
end

%% 10- Go back to the work space
if cfg.displayComment ==1;disp(['duneruo >>10 - Going back to the current folder ' cfg.currentPath ]);end
cd(cfg.currentPath)
if ~isfield(cfg,'writeLogFile'); cfg.writeLogFile = 0; end
if cfg.writeLogFile == 1; diary off; end
end