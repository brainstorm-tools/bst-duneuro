function cfg = bst_postprocess_eeg_lf(cfg)

%% substract the mean or not from the electrode
%     cfg.fem_eeg_lf
%     cfg.fem_meg_lf
%TODO : check the minifile parameters and adapt this code
if cfg.lfAvrgRef == 1
    if cfg.useTransferMatrix == 1
        if sum(cfg.lf_fem(1,:)) == 0
            disp('Transforming from elec1 reference to average reference');
            cfg.fem_eeg_lf  = cfg.fem_eeg_lf  - (mean(cfg.fem_eeg_lf,1));
        else
            disp('The average reference is the output of duneuro, please check the mini file');
        end    
    end
end

cfg.fem_eeg_lf = cfg.fem_eeg_lf;
end