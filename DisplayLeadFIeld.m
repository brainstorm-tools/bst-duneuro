function DisplayLeadFIeld(bstNodes)


% Display lead field as vectors ... Jonh Mosher Lead field representation
% discuss with him for the cas of the MEG and the notion of the reference
% and target electrode ... does it make sens for MEG ?

% Ask Francois how to activate the righ click on the figure in order to
% display the EEG sessor

ProtocolInfo = bst_get('ProtocolInfo');

% number of selected files
nbSelectedFiles = length(bstNodes);
indexColor = [0 0 1; 0 1 0;0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1]; % assuming 7 values of colors 

if nbSelectedFiles > length(indexColor)
    error('Too much files ... ')
end

for indFile = 1 :  nbSelectedFiles
    %% Get study description
    iStudy = bstNodes(indFile).getStudyIndex();
    sStudy = bst_get('Study', iStudy);
    
    %% Get Head model data
    iHeadModel = bstNodes(indFile).getItemIndex();
    %     HeadFileNames{indFile} = fullfile(bstDataPath,char(bstNodes(indFile).getFileName()));
    HeadModelFileNames{indFile} = fullfile(ProtocolInfo.STUDIES, sStudy.HeadModel(iHeadModel).FileName);        
    % Get the size of the Lead field matrix
    data = whos('-file', HeadModelFileNames{indFile}); % headModelData{indFile} = load(HeadModelFileNames{indFile})
    index = find(strcmp({data.name}, 'Gain')==1); % to check       data(index).name
    sizeOfGainMatrix{indFile} =  data(index).size; % 
    index = find(strcmp({data.name}, 'SurfaceFile')==1); % to check       data(index).name
    sourceSpaceFile = load(HeadModelFileNames{indFile},'SurfaceFile');
    % Get the modalities used in the study
    isEEG(indFile) = ~isempty(sStudy.HeadModel(iHeadModel).EEGMethod);
    isMEG(indFile) = ~isempty(sStudy.HeadModel(iHeadModel).MEGMethod);
    isSEEG(indFile) = ~isempty(sStudy.HeadModel(iHeadModel).SEEGMethod);
    isECOG(indFile) = ~isempty(sStudy.HeadModel(iHeadModel).ECOGMethod);      
    
    %% Get channel data
    ChannelFile{indFile} = bst_fullfile(ProtocolInfo.STUDIES, sStudy.Channel.FileName);
    data = whos('-file', ChannelFile{indFile});
    index = find(strcmp({data.name}, 'Channel')==1); % to check       data(index).name
    sizeOfSensorSpace{indFile} =  data(index).size(2); % 

     %% Get source space data
     %     % Source Space file name : Cortex File
    CortexModelFile{indFile}  = bst_fullfile(ProtocolInfo.SUBJECTS,sourceSpaceFile.SurfaceFile);
    data = whos('-file', CortexModelFile{indFile});
    index = find(strcmp({data.name}, 'Vertices')==1); % to check       data(index).name
    sizeOfSourceSpace{indFile} =  3*data(index).size(1); % multiply by 3 for unconstrained source
    
    %% Check the size of the different data     
    if (sizeOfGainMatrix{indFile}(1) ~= sizeOfSensorSpace{indFile}) && ...
            (sizeOfGainMatrix{indFile}(2) ~= sizeOfSourceSpace{indFile})
        error('The size of the LeadField Matrix does not match the number of source and sensors');
    end
    
    if indFile >1
        % TODO : does it make sens to compare different LF comming from different source space and sensor??
        if sum(sizeOfGainMatrix{indFile} == sizeOfGainMatrix{indFile-1}) ~= 2
            error('The two head models are differents');
        end
    end
    
    %% Load the data
    headModelData{indFile}  = load(HeadModelFileNames{indFile});
%     channelModelData{indFile}  = load(ChannelFile{indFile});
%     cortexModelData{indFile}  = load(CortexModelFile{indFile} );    
% LF(indFile) = {headModelData{indFile}.Gain(iEeg{indFile},:)}; LFnames(indFile) = {headModelData{indFile}.Comment }; LFcolor(indFile) = {indxColor(indFile,:)};
LF(indFile) = {headModelData{indFile}.Gain}; LFnames(indFile) = {[headModelData{indFile}.Comment ' LF'] }; LFcolor(indFile) = {indexColor(indFile,:)};
end

%% Ask the user for the modality to display
allModalities = {'EEG', 'MEG','sEEG','ECOG'}; % 1: EEG, 2 : MEG, 3 : sEEG, 4 : ECOG
notAvailable = [];
if isEEG == 0; notAvailable = [notAvailable, 1];end
if isMEG == 0; notAvailable = [notAvailable, 2] ;end
if isSEEG == 0; notAvailable = [notAvailable, 3];end
if isECOG == 0; notAvailable = [notAvailable, 4] ;end
allModalities(notAvailable) = '';
[selectedModality, isCancel] =  java_dialog('question', '<HTML><B> Select the modality <B>', ...
    'Display the Lead Field', [],allModalities, allModalities{1});
if isCancel ==1
    return;
end

%% Get Channel index and location per modalities:
% assuming that all the data has the same size (GAIN), the n we can fixe
% indFile = 1, 
hotindFile =1;
%% Load only one file for source and sensor
channelModelData{indFile}  = load(ChannelFile{indFile});
cortexModelData{indFile}  = load(CortexModelFile{indFile});


if isMEG(indFile) && strcmpi(selectedModality,'MEG')
    iMeg{indFile}   = good_channel(channelModelData{indFile}.Channel,[],'MEG');
    % Get the 3D location of the channels
    megloc{indFile}  = cat(2, channelModelData{indFile}.Channel(iMeg{indFile} ).Loc)';
    goodChannel = iMeg;
    sensorLocalisation = megloc;
end
if isEEG(indFile) && strcmpi(selectedModality,'EEG')
    iEeg{indFile}   = good_channel(channelModelData{indFile}.Channel,[],'EEG');
    % Get the 3D location of the channels
    eegloc{indFile}  = cat(2, channelModelData{indFile}.Channel(iEeg{indFile} ).Loc)';
    goodChannel = iEeg;
    sensorLocalisation = eegloc;
end
if isECOG(indFile) && strcmpi(selectedModality,'ECOG')
    iEcog{indFile}  = good_channel(channelModelData{indFile}.Channel,[],'ECOG');
    % Get the 3D location of the channels
    ecogloc{indFile}  = cat(2, channelModelData{indFile}.Channel(iEcog{indFile} ).Loc)';
    goodChannel = iEcog;
    sensorLocalisation = ecogloc;
end
if isSEEG(indFile) && strcmpi(selectedModality,'SEEG')
    iSeeg{indFile}  = good_channel(channelModelData{indFile}.Channel,[],'SEEG');
    % Get the 3D location of the channels
    seegloc{indFile}  = cat(2, channelModelData{indFile}.Channel(iSeeg{indFile} ).Loc)';
    goodChannel = iSeeg;
    sensorLocalisation = seegloc;
end


%% Update the LF according to the selected channels
for  indLF = 1 : length(LF)
LF_finale{indLF} = LF{indLF}(goodChannel{indFile},:);
end
LF = LF_finale;
%%Ask the user for the reference for the reference electrode to use 
if 1%~strcmp(selectedModality,'MEG')    
    referenceMode =  {'avgref','ref'}; % ask the user
    modality = { ' Use average reference  '};
    [res, isCancel] = java_dialog('checkbox', ...
        '<HTML><B> Do you want to use the average refence for the EEG (sEEG/ECOG) ?<BR>Otherwise you will choose the reference electrode <B>', 'Save Transfer Matrix', [], ...
        modality, [ones(1, 1)]);
    if isCancel ==1
        return;
    end
    if res == 1
        ref_mode = referenceMode{1};
    else
        ref_mode = referenceMode{2};
         % Ask for the reference electrode
            [res, isCancel] =  java_dialog('input', ...
                ['Set the Reference electrode [1 , '  num2str(length(sensorLocalisation{1})) ' ] : '], ...                
                'eference Electrode', [], '1');   
            if isCancel ==1
                return;
            end
            iref = str2num(res); %opts.referenceElectrode ; % 3;
    end   
end

%% Ask the user for the Source and Sink sensor
% Ask for the reference electrode
[res, isCancel] =  java_dialog('input', ...
    ['Set the target electrode [1 , '  num2str(length(sensorLocalisation{1})) ' ] : '], ...
    'Reference Electrode', [], '1');
if isCancel ==1
    return;
end
ielec = str2num(res); %opts.referenceElectrode ; % 3;% ielec = 15; %opts.targetElectrode; %2;
if (ielec>length(sensorLocalisation{1})) || (ielec<1)
    error('The selected value is out of range');
end
 
%% Start from here
    % Source Space
    GridLoc = cortexModelData{indFile}.Vertices';
    % sensor locations
    SensorLoc = sensorLocalisation{indFile}';
    %     % head model surfaces
    %     HeadFV =  struct('faces',cfg.bemFace{3}, 'vertices',cfg.bemNode{3});
    %     OuterFV = struct('faces',cfg.bemFace{2}, 'vertices',cfg.bemNode{2});
    %     InnerFV =  struct('faces',cfg.bemFace{1}, 'vertices',cfg.bemNode{1});
%     CortexFV =  struct('faces', cortexModelData{indFile}.Faces, 'vertices',cortexModelData{indFile}.Vertices);
    
    %     LF(1) = {headModelData{indFile}.Gain(iEeg{indFile},:)}; LFnames(1) = {'test'}; LFcolor(1) = {[0 0 1]};
    %     LF(1) = {opts.computed_solution1}; LFnames(1) = {opts.computed_solution1_name}; LFcolor(1) = {[0 0 1]};
    %     LF(2) = {opts.computed_solution2}; LFnames(2) = {opts.computed_solution2_name};LFcolor(2) = {[0 1 0]};
    %     if isfield(opts,'computed_solution3')
    %         LF(3) = {opts.computed_solution3}; LFnames(3) ={ opts.computed_solution3_name};LFcolor(3) = {[1 0 0]};
    %     end
    %     if isfield(opts,'computed_solution4')
    %         LF(4) = {opts.computed_solution4}; LFnames(4) ={ opts.computed_solution4_name};LFcolor(4) = {[1 0 1]};
    %     end
% end
% view the head
% figure
% clf reset
% h_h = patch(HeadFV,'facecolor','none','edgealpha',0.1);
% h_o = patch(OuterFV,'facecolor','none','edgealpha',0.15);
% h_i = patch(InnerFV,'facecolor','none','edgealpha',0.2);
% h_i = patch(CortexFV,'facecolor','none','edgealpha',0.2);

% [hFig, iDS, iFig] = 
% view_surface(CortexModelFile{indFile})
SurfaceFile = CortexModelFile{indFile};
SurfAlpha = 0.5;
SurfColor = [0.5 0.5 0.5] ;
hFig = 'NewFigure';
view_surface(SurfaceFile, SurfAlpha, SurfColor, hFig)
% axis vis3d, axis equal, rotate3d on, shg

% convenient for plotting command
X = GridLoc(1,:);Y = GridLoc(2,:); Z=GridLoc(3,:);

h = zeros(length(LF),1);
LeadField = cell(length(LF),1);

%% Plotting
hold on
for imodel = 1 : length(LF)
    switch ref_mode % {'avgref','ref'}
        case 'ref'
            LeadField{imodel} = LF{imodel}(ielec,:)- LF{imodel}(iref,:); %  leadfield row
            
            %legendName = [{['TargetElectrode' num2str(ielec) ] ,'ReferenceElectrode'} LFnames ];
            titleName = [' Reference Electrode index = ' num2str(iref) ];
        case 'avgref'
            AvgRef = mean(LF{imodel},1);            
            LeadField{imodel} = LF{imodel}(ielec,:) - AvgRef; %  leadfield row
            titleName = [' Reference Electrode average '];
            %legendName = [{['TargetElectrode' num2str(ielec) ]} LFnames ];
    end
    LeadField{imodel} = reshape(LeadField{imodel},3,[]); % each column is a vector
    U = LeadField{imodel}(1,:);
    V = LeadField{imodel}(2,:);
    W =LeadField{imodel}(3,:); % convenient
    h(imodel) = quiver3(X,Y,Z,U,V,W, 5);
    set(h(imodel),'linewidth',1,'color',LFcolor{imodel},'linewidth',1); % mark the electrode)
end

%hl =legend(h,LFnames);
% add the sensor locations
% 
% hold on
% h(imodel+1) = plot3(SensorLoc(1,:),SensorLoc(2,:),SensorLoc(3,:),'wo','markersize',10);


hold on
h(imodel+1) = plot3(SensorLoc(1,ielec),SensorLoc(2,ielec),SensorLoc(3,ielec),'r*','markersize',15,'linewidth',2); % mark the electrode
hold on; plot3(SensorLoc(1,ielec),SensorLoc(2,ielec),SensorLoc(3,ielec),'ro','markersize',10,'linewidth',2); % mark the electrode
% hl =legend(h,[LFnames,{'Sensor'},{'Target'}], 'Location','northwestoutside','TextColor','w');
hl =legend(h,[LFnames,{'Target'}], 'TextColor','w');
legend('boxoff')

if strcmp(ref_mode,'ref')
    hold on
    h(imodel+2) = plot3(SensorLoc(1,iref),SensorLoc(2,iref),SensorLoc(3,iref),'yo','markersize',10,'linewidth',2); % mark the electrode
    h(imodel+2) = plot3(SensorLoc(1,iref),SensorLoc(2,iref),SensorLoc(3,iref),'y+','markersize',15,'linewidth',2); % mark the electrode
%     hl =legend(h,[LFnames,{'Sensor'},{'Target'},{'Reference'}], 'Location','northwestoutside','TextColor','w');
    hl =legend(h,[LFnames,{'Target'},{'Reference'}], 'TextColor','w');
    legend('boxoff')
end
set(hl,'fontsize',14)
% rotate3d on;
title([titleName ' LF'])
%% check this one
% bst_display_leadField_vector(cfg,opts)
end
