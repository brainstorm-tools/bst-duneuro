function cfg = bst_run_duneuro_cmd(cfg)
% cfg = bst_run_duneuro_cmd(cfg)

% set the arguments and run the duneuro computation.
% Author , Takfarinas MEDANI, December 2019,
% Updates :   - Remove the previous version and keep only the version3                                  
%                   - Add the option to call the binaries even from weird  path 
%                      https://github.com/brainstorm-tools/bst-duneuro/issues/1#issuecomment-583076206
                 
%% 0- Set minimal configuration
if ~isfield(cfg,'BstDuneuroVersion');  cfg.BstDuneuroVersion = 3; end


%% 1- Get the Duneuro configuration file.
DuneuroConfiguraionFile = cfg.mini_filename;

%% 2- Construct & run the Duneuro comande line.
% version 3 the modality is included on the ini file
% This version take account of the weired files names and paths
[filepath,name,ext] = fileparts(cfg.cmd);
bstfpath = filepath;bstfname = name; bstext = ext;
[filepath,name,ext] = fileparts(DuneuroConfiguraionFile);
inifpath = filepath;inifname = [name ext]; arghelp = '--h';
% callStr5 = [ '"' bstfpath  filesep bstfname bstext '"' ' ' '"' inifpath...
% filesep inifname '"' ' ' '"' arghelp '"'];   % Call duneuro Help 
callStr6 = [ '"' bstfpath  filesep bstfname  '"' ' ' '"' inifpath filesep inifname '"'];
 
[status,cmdout] = system(callStr6);

%% 3- Check status & outputs
if status ~= 0
    duneuro_logfile = fullfile(cfg.pathOfTempOutPut,'binaries_duneuro_log.txt');
    fid = fopen(duneuro_logfile , 'w');
    fprintf(fid, '%st', cmdout);
    fprintf(fid, '%st', ' Load the file duneuro_configuration for manual debuging/data checking');
    fclose(fid);
    save(fullfile(cfg.pathOfTempOutPut,'duneuro_configuration'), 'cfg');
    fclose('all');
    if cfg.runFromBst == 1; bst_progress('stop'); end
    error('Something was wrong during duneuro computation.\nPlease check %s', duneuro_logfile);
else
    if cfg.displayComment ==1;     disp('DuNeuro FEM computation completed without error '); end
end

end
