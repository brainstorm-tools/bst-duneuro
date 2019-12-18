function cfg = bst_prepare_minifile(cfg)
% cfg = bst_prepare_minifile(cfg)
% set the configration values used by duneuro computation.
% This function will also write the mini file.

% Takfarinas MEDANI, October 2019,

% Set mini file to default parameter /configuration
% Put all the paramater in the cfg structure this step could  be modified
% from the bst gui via cfg structure. 
cfg = bst_set_minifile(cfg);

% Write the file
cfg.mini_filename = [ cfg.filename '_minifile.mini' ];
write_duneuro_minifile2(cfg);

end