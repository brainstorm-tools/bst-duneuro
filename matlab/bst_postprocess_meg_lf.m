function cfg = bst_postprocess_meg_lf(cfg)
%% Compute the total magnetic field 
dipoles_pos_orie = [kron(cfg.sourceSpace,ones(3,1)), kron(ones(length(cfg.sourceSpace),1), eye(3))];

% a- Compute the MEG Primary Magnetic field  % apply formula of Sarvas
%primary B-field
[Bp] = compute_B_primary(cfg.coilsLoc, dipoles_pos_orie, cfg.coilsProjection);

% b- The total magnetic field B = Bp + Bs;
%  full B-field
Bs =  cfg.fem_meg_lf;

mu = 4*pi*1e-4; % check the value of the units maybe it needs to be mu = 4*pi*1e-7
Bfull = (mu/(4*pi)) * (Bp - Bs);

% c- Apply the weight :
if cfg.runFromBst == 1
    [channelIndex] = unique(cfg.MegChannel(:,1));
    nbChannel = length(channelIndex);
    wighted_B =zeros(nbChannel,size(Bfull,2));
    for iCh = 1 : nbChannel
        communChannel = find(iCh==cfg.MegChannel);
        BcommunChannel = Bfull(communChannel(:),:);
        WcommunChannel =  cfg.MegChannel(communChannel(:), 8: end);
        wighted_B(iCh,:) = sum (BcommunChannel.*WcommunChannel,1);
    end    
    cfg.fem_meg_lf = wighted_B;
end
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