function bst_view_tensor(iSubject)

%% SECTION 1 : Get the data
disp('Generate fem tensors')

%% Get the input data
% get the mesh
% get the EIG-DTI
% get the isotropic conductivity
% if anisotrop as for the method to use
% Get Protocol information
ProtocolInfo     = bst_get('ProtocolInfo');
ProtocolSubjects = bst_get('ProtocolSubjects');
% Default subject
if (iSubject == 0)
    sSubject = ProtocolSubjects.DefaultSubject;
    % Normal subject
else
    sSubject = ProtocolSubjects.Subject(iSubject);
end

%% Get the mesh file
% Get the conductivity values
FemFiles = file_fullpath(sSubject.Surface(sSubject.iFEM).FileName);
% Get name and the number of  layers
% Load the mesh
femHead=  load(FemFiles);
numberOfLayer = length(femHead.TissueLabels);

%% Display the tensore within a mesh slice : defined by a plan
% define the cutting plan
z0 = mean(femHead.Vertices(:,3)); % range(femHead.Vertices(:,3))
z0 = 0.038;
plane = [min(femHead.Vertices(:,1)) min(femHead.Vertices(:,2)) (z0 + z0/2)
    min(femHead.Vertices(:,1)) max(femHead.Vertices(:,2)) (z0 + z0/2)
    max(femHead.Vertices(:,1)) min(femHead.Vertices(:,2)) (z0 + z0/2)];
% 
% y0 = 0.022;
% plane=[min(femHead.Vertices(:,1)) y0 min(femHead.Vertices(:,3))
%     min(femHead.Vertices(:,1)) y0 max(femHead.Vertices(:,3))
%     max(femHead.Vertices(:,1)) y0 min(femHead.Vertices(:,3))];
% 
% 
% x0 = 0.0;
% plane=[x0  min(femHead.Vertices(:,2))  min(femHead.Vertices(:,3))
%     x0  max(femHead.Vertices(:,2)) max(femHead.Vertices(:,3))
%     x0  min(femHead.Vertices(:,2)) min(femHead.Vertices(:,3))];

% run qmeshcut to get the cross-section information at z=mean(node(:,1))
% use the x-coordinates as the nodal values

[cutpos, cutvalue, facedata,elemid] = ...
                    qmeshcut(femHead.Elements,femHead.Vertices,zeros(length(femHead.Vertices),1),plane);

figure;
plotmesh(femHead.Vertices,[femHead.Elements femHead.Tissue], 'facealpha', 0.3,'edgecolor','none');
hold on
plotmesh(femHead.Vertices,[femHead.Elements(elemid,:) femHead.Tissue(elemid,:)]);
view([0 90 0])
title('Slice of mesh where the tensor will be displayed')

%% Display intererface
% 1- ask user which tissue to include
% ask for the layer to consider as anisotrop
[res, isCancel] = java_dialog('radio', ...
    '<HTML> Select the layers to consider for displaying the tensor <BR>', 'Select Tissues', [], ...
    [femHead.TissueLabels {'all'}], 1);
if isCancel;         return;    end
includedTissueIndex = (res);
if includedTissueIndex > length(unique(femHead.TissueLabels))
    includedTissueIndex = unique(femHead.Tissue);
end
elemToDisplay = num2str(length(find(femHead.Tissue(elemid) == includedTissueIndex')));

% linespace : equally spaced points between all the tensors 
[res, isCancel] = java_dialog(  'input' , ['step : equally space the '  num2str(elemToDisplay) ' tensor by a step =  ' ],...
    'Reduce the number of tensor to display' ,[],'10');
if isCancel;         return;    end
regularStep = str2double( res);


% display the head model
[res, isCancel] = java_dialog('radio', '<HTML><B> Display the head model (overlay with tensors) <B>', ...
        'Select Head Model', [],{'Yes','No'}, 1);
if isCancel;         return;    end
if res == 1
    displayHeadModel = 1;
else
    displayHeadModel = 0;
end

% display the tensors
[res, isCancel] = java_dialog('radio', '<HTML><B> Display the tenosr as <B>', ...
        'Select tensor view', [],{'Ellipse','Arrow (main eigen vector, useful only on anistropic case) '}, 1);
if isCancel;         return;    end
if res == 1
    displayAsEllipse = 1;
    displayAsArrow = 0;
    
else
    displayAsEllipse = 0;
    displayAsArrow = 1;
end

cfg = [];
cfg.node = femHead.Vertices; % list of nodes 
cfg.elem = [femHead.Elements femHead.Tissue] ; % list of element with lable 
cfg.elemid = elemid;
cfg.elemid = elemid(find(sum((femHead.Tissue(elemid) == includedTissueIndex'),2))); % reduced to the element to display

cfg.elemid = cfg.elemid(1:regularStep:end); % reduced to the element to display    find(cfg.elem(cfg.elemid,5) ==3)
cfg.eigen = femHead.tensors; % thei eigen data  
cfg.elem_centroide = femHead.tensors.position; % the centoride of the element
cfg.ellipse = displayAsEllipse; % display as an ellipse
cfg.arrow = displayAsArrow; % display as an arrow
cfg.plotMesh = displayHeadModel; % plot the head model
bst_display_fem_tensors(cfg) % the displaying function


% 
% view([0 -90 0])
% axis([ -0.063 0.095 ...
%     min(cfg.node(:,2))  max(cfg.node(:,2)) ...
%     -0.02  0.105])
% hold on
% plotmesh(femHead.Vertices,femHead.Elements(elemid,:),'facealpha',0.2,'edgecolor','k');
% view([0 90 0])
% 
% hold on
% plotmesh(femHead.Vertices,femHead.Elements(:,:),'facealpha',0.2);
% 
% view([0 0 90])
% 
% 
% figure;
% plotmesh(femHead.Vertices,femHead.Elements,'facealpha',0.2,'edgecolor','none');
% view([0 0 90])

% 
% figure;
% plotmesh(femHead.Vertices,femHead.Elements(elemid,:),'facealpha',0.2);
% view([0 0 90])

end