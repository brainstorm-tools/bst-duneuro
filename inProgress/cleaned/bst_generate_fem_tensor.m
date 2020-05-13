function bst_generate_fem_tensor(iSubject)

disp('Generate fem tensors')


% TODO : need to extend this function for multiple tissue anisotrop (wm and gm or skull for example)

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

% Get the EIG-DTI
% Load the DTI tensor and convert to conductivity tensors
temp1 = []; temp2 = [];
for ind = 1 : length(sSubject.Anatomy)
    temp1 = [temp1 ~isempty(strfind(sSubject.Anatomy(ind).Comment, 'DTI-EIG')) ];
    temp2 = [temp2 ~isempty(strfind(sSubject.Anatomy(ind).Comment, 'tissue')) ];
end
iEigenTissue = find(temp1);
iMaskTissue = find(temp2);

eigenFile = file_fullpath(sSubject.Anatomy(iEigenTissue).FileName);
maskFile = file_fullpath(sSubject.Anatomy(iMaskTissue).FileName);


% Get the scalar conductivity values
default_iso_conductivity = get_standard_conductivity(numberOfLayer) ;


%% Ask which tensor tensors to use
[res, isCancel] =  java_dialog('radio', '<HTML><B> Tensor model <B>', ...
    'Set Tensor model', [],{'Isotropic [expeimental]', 'Anisotropic' }, 1);
if isCancel ==1;     return; end
if res == 2
    isAnisotropic = 1;
else
    isAnisotropic = 0;
end

%% Ask for the conductivity values
for ind = 1 : length(default_iso_conductivity)
    cond{ind} = num2str(default_iso_conductivity(ind));
end

[res, isCancel] = java_dialog('input', femHead.TissueLabels, 'Isotropic conductivities values', [], cond);
if isCancel;     return; end

if length(default_iso_conductivity) == 1
    default_iso_conductivity(1) =  str2num(res(1));
else
    for ind = 1 : length(default_iso_conductivity)
        default_iso_conductivity(ind) = str2num(res{ind});
    end
end


%% Ask for anisotropic options
if isAnisotropic == 1
%     anisotropyOption = [];
    % ask for the layer to consider as anisotrop ==> TODO maybe extended
    % for more than one layer
    [res, isCancel] = java_dialog('radio', ...
        '<HTML>Select the layers to consider for the anisotropy <BR>', 'Select Volume', [], ...
        femHead.TissueLabels, 1);
    if isCancel;         return;    end
    anisoTissueIndex = res;
    anisotropyOption.anisoTissueIndex = anisoTissueIndex;
    % We can mix between the methods ? ask francois if we can buif a panel that can propose this option
    % for example the WM will be from the DWI and the skull can be
    % artifical
    
    % ask for the method
    AnisotropyModel =  {'Effective Medium Approach [EMA]','EMA with Volume Constraint','Artificial/simulated anisotropy'};
    [res, isCancel] =  java_dialog('radio', '<HTML><B> Tensor Estimation : Select Method <B>', ...
        'FEM Source Model', [],AnisotropyModel, 1);
    if isCancel ==1;         return;    end
    anisotropyMethod = res;
    anisotropyOption.anisotropeMethod = anisotropyMethod;
    % if artificial ask for the ratio between longitudinal and transverse
    % conductivity
    isArtificial = 0;
    if anisotropyMethod == 3
        isArtificial = 1;
        artificialRatio = 10;
        % ask for the ratio
        [res, isCancel] = java_dialog('input', 'ratio: long / tran' , 'Artificial Anisotropy ratio', [], num2str(artificialRatio));
        artificialRatio = str2num(res);
        if isCancel ==1;         return;    end
        anisotropyOption.artificialRatio = artificialRatio;
        
        % ask for constrained method
        [res, isCancel] = java_dialog('radio', ...
            '<HTML> Constraint method <BR>', 'Select Volume', [], {'Wang''s constraint [sig_r*sig_t = sig^2]',...
            'Wolter''s consraint [(4/3)*pi*(sig_r*sig_t^2) = (4/3)*pi*(sig^3)]'}, 1);
        if isCancel ==1;         return;    end
        if res == 1
            constraintMethod = 1; % wang
        else
            constraintMethod = 2; % wolters
        end
        anisotropyOption.constraintMethod = constraintMethod;
        
        % ask for orientation to use
        [res, isCancel] = java_dialog('radio', ...
            '<HTML> Tensors Orientation <BR>', 'Select Vectors', [], {'Artificial Eigen Vectors',...
            'DTI Eigen Vectors [Need the DTI-EIG in the database]'}, 1);
        if isCancel ==1;         return;    end
        if res == 1
            isArtificalOrientation = 1;
        else
            isArtificalOrientation = 0;
        end
    end
end

% get the index of the max conductivity for normalisation
[tmp,  ref_tissu_index] = max(default_iso_conductivity);
anisotropyOption.ref_tissu_index = ref_tissu_index;
%%%%%%%%%%%%%%%%%%%%% %%
%% Build the tensor for the isotropic tissues
cfg = [];
cfg.elem = [femHead.Elements femHead.Tissue];
cfg.node = femHead.Vertices;
cfg.conductivity = default_iso_conductivity;
cfg = bst_generate_tensor_on_elem(cfg);

% output results
tensors = cfg.eigen;
tensors.position = cfg.elem_centroide;
tensors.param = [];

% The function ends here in the case of the isotropic model

%% Get the anisotropy information
if isAnisotropic == 1
    temp = [];
    if ~ isArtificial % not artificial values
        
        % Load the DTI and the Masks
        %         eigenFile = file_fullpath(sSubject.Anatomy(iEigenTissue).FileName);
        %         maskFile = file_fullpath(sSubject.Anatomy(iMaskTissue).FileName);
        sEigDti =  in_mri_bst(eigenFile);
        sMask =  in_mri_bst(maskFile);
        
        DTI = {};
        DTI{1}.hdr.dim.pixdim = sMask.Header.dim.pixdim;
        DTI{1}.anatomy = sEigDti.Cube(:,:,:,1:3);
        DTI{2}.anatomy = sEigDti.Cube(:,:,:,4:6);
        DTI{3}.anatomy = sEigDti.Cube(:,:,:,7:9);
        DTI{4}.anatomy = sEigDti.Cube(:,:,:,10);
        DTI{5}.anatomy = sEigDti.Cube(:,:,:,11);
        DTI{6}.anatomy = sEigDti.Cube(:,:,:,12);
        
        %% a- Convert the coordinate to SCS and project the tensors the centroide
        % of the element
        [V1rot,V2rot,V3rot, L1a, L2a, L3a ] = bst_tensor_cs_convert(sMask, femHead, DTI,anisotropyOption.anisoTissueIndex);
        
        %% Convert diffusion to conductivity tensors
        
        % Call the main function : bst_compute_anisotropy_tensors
        % METHOD 1 : direct transformation approach [Güllmar et al   NeuroImage 2010]
        % METHOD 2 : direct transformation approcah with volume constraint  [Güllmar et al NeuroImage 2010]
        % METHOD 3 : Artificial anisotropy with volume constraint [Güllmar et al NeuroImage 2010]
        % METHOD 4 : The volume normalized approcah  [Won Hee Lee et all  IEEE 2015]
        % METHOD 5 : Wang-constraint model and volume-constraint model  [Wolters CH 2006]
        % METHOD 6 : Wang-constraint model and volume-constraint model  [Wolters CH 2006]
        % METHOD 7 : As SimBio toolbox  [Rullmann et al 2008 / Vorwerk et al 2014 ]
        usedMethod = [2, 7];
        options.aniso_method = usedMethod(anisotropyOption.anisotropeMethod);%  [ 1, 7]
        options.aniso_tissu_index = anisotropyOption.anisoTissueIndex;
        options.ref_tissu_index = anisotropyOption.ref_tissu_index;
        options.transversal_factor = 1;
        options.longitudinal_factor = 1;
    
    
        [aniso_conductivity, anisotropicTensor,param ] = ...
            bst_compute_anisotropy_tensors(femHead,  default_iso_conductivity, tensors, L1a,L2a,L3a,V1rot,V2rot,V3rot,options);
        
        
        %% The final output
        tensors = anisotropicTensor;
        tensors.param = param;
        tensors.position = cfg.elem_centroide;

    else % it's artificial conductivity
        %        IMPORTANT :  Add two option : artificial orientation, or/and artificial eigen valued
        iso_conductivity = default_iso_conductivity(anisotropyOption.anisoTissueIndex);
        ratio = artificialRatio;
        % compute the new eigen value
        
        if  anisotropyOption.constraintMethod == 1 % Wang's
            lm1 = (iso_conductivity^2 * ratio)^(1/2); % sigma longitidunal
            lm2 = (iso_conductivity^2 / ratio)^(1/2);     % sigma transversal
            lm3 = lm2;
        end
        
        if anisotropyOption.constraintMethod == 2 % Wolters's
            lm1 = (iso_conductivity^3 * ratio^2)^(1/3); % sigma longitidunal
            lm2 = (iso_conductivity^3 / ratio)^(1/3);     % sigma transversal
            lm3 = lm2;
            % Check
            if round((lm1*lm2*lm3),8) == round(iso_conductivity^3,8)
                disp('Volume is conserved');
            end
        end
        
        % get the desired orientation
        if isArtificalOrientation
            % keep the orientation within the tensors.eigen_vector
            
        else % realistric orientation
            sEigDti =  in_mri_bst(eigenFile);
            sMask =  in_mri_bst(maskFile);
            
            DTI = {};
            DTI{1}.hdr.dim.pixdim = sMask.Header.dim.pixdim;
            DTI{1}.anatomy = sEigDti.Cube(:,:,:,1:3);
            DTI{2}.anatomy = sEigDti.Cube(:,:,:,4:6);
            DTI{3}.anatomy = sEigDti.Cube(:,:,:,7:9);
            DTI{4}.anatomy = sEigDti.Cube(:,:,:,10);
            DTI{5}.anatomy = sEigDti.Cube(:,:,:,11);
            DTI{6}.anatomy = sEigDti.Cube(:,:,:,12);
            
            %% a- Convert the coordinate to SCS and project the tensors the centroide
            % of the element
            [V1rot,V2rot,V3rot, L1a, L2a, L3a] = bst_tensor_cs_convert(sMask, femHead, DTI,anisotropyOption.anisoTissueIndex);
        end
        
        
        % The output
        allOfIndexAnisoTissu = find(femHead.Tissue == anisotropyOption.anisoTissueIndex);
        for ind = 1 : length(allOfIndexAnisoTissu)
            if isArtificalOrientation
                V = tensors.eigen_vector{allOfIndexAnisoTissu(ind)};
            else % use realistric orientation
                V = [V1rot(ind,:)' V2rot(ind,:)' V3rot(ind,:)'];
                %                 V = [V1rot(ind,:); V2rot(ind,:); V3rot(ind,:)];
            end
            %             V = eigenVectors{allOfIndexAnisoTissu(ind)};
            %         aniso_conductivity(:,:,allOfIndexAnisoTissu(ind)) =  V* diag([lm1, lm2,lm3])*V';
            tensors.eigen_vector{allOfIndexAnisoTissu(ind)} = V;
            tensors.eigen_value{allOfIndexAnisoTissu(ind)} = diag([lm1 lm2 lm3]);
            meanConductivity(ind) = mean([lm1 lm2 lm3]);
        end
    end
    % summerize the outputs
    tensors.param.anisoTissueIndex = anisotropyOption.anisoTissueIndex;
    tensors.param.anisotropeMethod = anisotropyMethod;
end



%% Save to the database / update the headmodel with the new field tensor
Comment = femHead.Comment;
Vertices = femHead.Vertices;
Elements = femHead.Elements;
Tissue = femHead.Tissue;
TissueLabels = femHead.TissueLabels;
History = femHead.History;

save(FemFiles,'Comment',...
    'Vertices',...
    'Elements',...
    'Tissue',...
    'TissueLabels',...
    'History',...
    'tensors') ;

disp('Generate fem tensors ... done')

end