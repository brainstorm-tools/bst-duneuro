function cfg = bst_postprocess_meg_lf(cfg)
%% Compute the total magnetic field B = Bp + Bs
% a- Compute the MEG Primary Magnetic field  % apply formula of Sarvas
dipoles_pos_orie = [kron(cfg.sourceSpace,ones(3,1)), kron(ones(length(cfg.sourceSpace),1), eye(3))];
[Bp] = compute_B_primary(cfg.coilsLoc, dipoles_pos_orie, cfg.coilsProjection);
% b- Get the MEG Primary Magnetic field
% [Bs] = cfg.fem_meg_lf;
% c- The total magnetic field B = Bp + Bs;
Bs =  cfg.fem_meg_lf;
B = Bp + Bs;

% d- Apply the weight :
if cfg.runFromBst == 1
    [channelIndex] = unique(cfg.MegChannel(:,1));
    nbChannel = length(channelIndex);
    wighted_B =zeros(nbChannel,size(B,2));
    for iCh = 1 : nbChannel
        communChannel = find(iCh==cfg.MegChannel);
        BcommunChannel = B(communChannel(:),:);
        WcommunChannel =  cfg.MegChannel(communChannel(:), 8: end);%(communChannel',:);
        wighted_B(iCh,:) = sum (BcommunChannel.*WcommunChannel,1);
    end
    
    cfg.fem_meg_lf = wighted_B;
end
% cfg = rmfield(cfg,'lf_fem');
end


% channelIndex    = cfg.MegChannel(:,1);
% channelWeight = cfg.MegChannel(:, 8: end);
% % % testing for the comming data
% channelIndex = [ones(1,4), 2*ones(1,3), 3*ones(1,1), 4*ones(1,3), 5*ones(1,5),]'; %
% channelIndex = [ones(1,1), 2*ones(1,3), 3*ones(1,1)]'; %
% channelWeight = repmat(1,size(channelIndex));
% [channelIndex channelWeight];
% Bt = rand(length(channelIndex),2);
% [~, groupChanIndex] = unique(cfg.MegChannel(:,1));

% [channelIndex channelWeight B]
%     if iCh < nbChannel
%         communChannel =     ([(groupChanIndex(iCh):groupChanIndex(iCh+1)-1)]);
%     else
%         communChannel =    (groupChanIndex(iCh):length(cfg.MegChannel(:,1)));
%     end