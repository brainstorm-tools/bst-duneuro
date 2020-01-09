function cfg = bst_postprocess_lf(cfg)
% call the post-processing of the LF.
% meg : apply the wight oc sensor (gradiometre, magnetometre)
% eeg : apply the desired reference/ reference == > TODO

if strcmp(cfg.modality,'eeg')
    cfg = bst_postprocess_eeg_lf(cfg);
end

if strcmp(cfg.modality,'meg')
    cfg = bst_postprocess_meg_lf(cfg);
end

if strcmp(cfg.modality,'meeg')
    cfg = bst_postprocess_meg_lf(cfg);
    cfg = bst_postprocess_eeg_lf(cfg);    
end

end


