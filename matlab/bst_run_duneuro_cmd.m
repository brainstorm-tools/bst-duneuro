function cfg = bst_run_duneuro_cmd(cfg)
% cfg = bst_run_duneuro_cmd(cfg)

% set the arguments and run the duneuro computation.
% Author , Takfarinas MEDANI, December 2019,


%% Se minimal configuration 
if ~isfield(cfg,'BstDuneuroVersion');  cfg.BstDuneuroVersion = 2; end 

% the configuration file.
if cfg.runFromBst
    arg1 = fullfile(cfg.pathOfTempOutPut, cfg.mini_filename); 
else
    arg1 = cfg.mini_filename;
end
arg2 = [' '  '--' cfg.modality];

% run the command
% To avoid the display on the terminal : use : evalc('[status,cmdout] = system([cfg.cmd '.exe' ' '  arg1]);;')
% version 1
if cfg.BstDuneuroVersion == 1
    if ispc
        [status,cmdout] = system([cfg.cmd '.exe' ' '  arg1]);
    elseif isunix
        [status,cmdout] = system(['./' cfg.cmd ' '  arg1]);
    elseif ismac
        [status,cmdout] = system(['./' cfg.cmd ' '  arg1]);
    end
end
% version 2
if cfg.BstDuneuroVersion == 2
    if ispc
        [status,cmdout] = system([cfg.cmd '.exe' ' '  arg1 arg2])  ;
    elseif isunix
        [status,cmdout] = system(['./' cfg.cmd ' '  arg1 arg2]);
    elseif ismac
        [status,cmdout] = system(['./' cfg.cmd ' '  arg1 arg2]);
    end
end

%% Check status
if status ~= 0
    duneuro_logfile = 'duneuro_log.txt';
    fid = fopen(duneuro_logfile , 'w');
    fprintf(fid, '%st', cmdout);
    fprintf(fid, '%st', ' Load the file duneuro_configuration for manual debuging');
    fclose(fid);
    save('duneuro_configuration', 'cfg');
    error('Something was wrong during duneuro computation.\nPlease check %s', duneuro_logfile);
else
    disp('DuNeuro FEM computation completed without error ');
end

end
