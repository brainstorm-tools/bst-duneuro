function cfg = bst_prepare_sensor_model(cfg)
% cfg = bst_prepare_sensor_model(cfg)
% Prepare the sensors for the duneuro computation

% Takfarinas MEDANI, December 2019,


%% Set minimal configuration
if ~isfield(cfg,'runFromBst'); cfg.runFromBst = 0; end  

%% Write the files

% Write the EEG electrode file in the case of EEG and MEEG
if strcmp(cfg.modality,'eeg') || strcmp(cfg.modality,'meeg')
    cfg.electrode_filename = 'electrode_model.txt';
    write_duneuro_electrode_file(cfg.channelLoc, cfg.electrode_filename);
end

% Write the MEG sensors file in the case of MEG and MEEG
if strcmp(cfg.modality,'meg')  || strcmp(cfg.modality,'meeg')
    if cfg.runFromBst == 1 % otherwise, read the files autrement
        %MegChannel =  [iChan, sChan.Loc(:,iInteg)', sChan.Orient(:,iInteg)', sChan.Weight(iInteg)];
        cfg.coilsLoc =    cfg.MegChannel(:,2:4);
        cfg.coil_filename = 'coil_model.txt';
        write_duneuro_coil_file(cfg.coilsLoc, cfg.coil_filename);
        % TODO : need to fixe ... it seems that duneuro uses 3 orientations per
        % coil and brainstrom only one orientation .... maybe it needs local
        % projection (scalar product to the correct orientation)
        cfg.coilsProjection =    cfg.MegChannel(:,5:7);
        cfg.projection_filename = 'projection_model.txt';
        write_duneuro_projection_file(cfg.coilsProjection, cfg.projection_filename);
        % elseif strcmp(cfg.modality,'seeg')
        % TODO
    end    
end

end