function [aniso_conductivity, anisotropicTensor ] = bst_compute_anisotropy_tensors(femHead,  default_iso_conductivity, IsotropicTensor, L1a,L2a,L3a,V1rot,V2rot,V3rot,options)
% Compute or convert the DTI tensors to anisotropic conductivity tensors.
% Input :
% femHead : FEM mesh within brainstorm format
% default_iso_conductivity : conductivity values of the tissues, shoud have
% the same length as the number of tissues of the femHead
% IsotropicTensor : The conuctivity rensor compited on the femHead model (outpuf from bst_generate_tensor_on_elem)
% based on the scal value of the conductivity
% L1a to L3a : eigenvalue computed from brainsuite and then interpolated to the fem element
% V1a to V3a : eigenv vector interpolated to the centroide of the element.
% options : options that contains the method to use and the index of the
% anisotropic tissue.
%           options.aniso_method : integer that specify the index of the
%           method
%           options.aniso_tissu_index: integer that refer to the anisotrop
%           tissue
%           options.ratio_long_trans =  ratio between the
%           longitidunal and transversal conductivities (method : 4, 5 and 6 )
% List of implemented method
% output :
% aniso_conductivity : 3D matrix containing the tensors
% anisotropicTensor : cell containing the eigen vector and value computed
% on each element (useful foe the display) ==> same format as IsotropicTensor

% Takfarinas MEDANI

%% Get the inputs arguments
if nargin < 10 % the options are not specified
    aniso_method = 7; % select the method to use in order to compute the anisotropy
    aniso_tissu_index = 1; % mainly the white matter but it can be also gray matter and skull
    ref_tissu_index = 3; % mainlt the CSF hat have the biggest values
    transversal_factor = 1;
    longitudinal_factor = 10;
    ratio = longitudinal_factor/transversal_factor;
else
    aniso_method = options.aniso_method; % select the method to use in order to compute the anisotropy
    aniso_tissu_index = options.aniso_tissu_index; % mainly the white matter but it can be also gray matter and skull
    ref_tissu_index = options.ref_tissu_index;
    transversal_factor = options.transversal_factor;
    longitudinal_factor = options.transversal_factor;
    ratio = longitudinal_factor/transversal_factor;
end

iso_conductivity = default_iso_conductivity(aniso_tissu_index);  % get this value from the input
max_conductivity = default_iso_conductivity(ref_tissu_index);  %should be the CSF

%% 1- Ensure that all the eigen values are positives 
L1a = abs(L1a); L2a = abs(L2a); L3a = abs(L3a);
% The tensors will be symetric by construction
% ensures that all conductivity tensors are positive definite
%%

%% Initialization
aniso_conductivity= [];
dti_factor = 1; % maybe it needs to convert the units ... need check
L1a = L1a * dti_factor;
L2a = L2a * dti_factor;
L3a = L3a * dti_factor;
isIsotrope = [];

% index for the anisotrpic tissu
ind = 1;
meanConductivity =[];
maxL1 = max(L1a); minL1 =  min(L1a);  meanL1 =  mean(L1a);
maxL2 = max(L2a); minL2 =  min(L2a);  meanL2 =  mean(L2a);
maxL3 = max(L3a); minL3 =  min(L3a);  meanL3 =  mean(L3a);


for globalIndex = 1 : length(femHead.Tissue)
    % Check if the tissue ID if it's belonging to the anisotropy tissue
    if femHead.Tissue(globalIndex) == aniso_tissu_index % 1 is the label of the wm
        L1 = L1a(ind);
        L2 = L2a(ind);
        L3 = L3a(ind);              
        L = diag([L1, L2, L3]);
        %% Ensure maximal ratio of 10 between the  largest and smallest conductivity eigenvalues
        largest_ratio = 10;
        if  (sum(sum(L)) ~= 0)
            if ((max([L1 L2 L3]) / min([L1 L2 L3])) >= largest_ratio)
                [vMax, iMax] = max([L1 L2 L3]);
                [vMin, iMin] = min([L1 L2 L3]);
                L(iMin,iMin) = L(iMax,iMax)/largest_ratio;
                % update the values
                L1 = L(1,1);L2 = L(2,2);L1 = L(3,3);
            end
        end
        % Construct the DTI tensors 
        V = [V1rot(ind,:)', V2rot(ind,:)', V3rot(ind,:)'];
        T = V*L*V'; % the dti tensor      
        
        %% Assigne the tensors for each element 
        if sum(sum(L)) == 0 % useful in the case where BDP fails or it's not part of the wm mask
            % in this case the isotropic model is used will be used
            %             disp('case 1 : sum(L) = 0')
            isIsotrope = [isIsotrope ind];
            aniso_conductivity(:,:,ind) = diag([iso_conductivity,iso_conductivity,iso_conductivity]);
            eigen.eigen_vector{ind} = IsotropicTensor.eigen_vector{globalIndex};
            eigen.eigen_value{ind} =  IsotropicTensor.eigen_value{globalIndex};
            
        elseif L1 < L2 || L1<L3 % this may not happen ... BDP fails ==> replace by isotropic case
            %             disp(' case 2 : L1 < L2 || L1<L3')
            isIsotrope = [isIsotrope ind];
            aniso_conductivity(:,:,ind) = diag([iso_conductivity,iso_conductivity,iso_conductivity]);
            eigen.eigen_vector{ind} = IsotropicTensor.eigen_vector{globalIndex};
            eigen.eigen_value{ind} =  IsotropicTensor.eigen_value{globalIndex};
        elseif L1 == 0 % this may not happen ... BDP fails ==> replace by isotropic case
            %             disp('case 3 : L1 == 0 ')
            isIsotrope = [isIsotrope ind];
            aniso_conductivity(:,:,ind) = diag([iso_conductivity,iso_conductivity,iso_conductivity]);
            eigen.eigen_vector{ind} = IsotropicTensor.eigen_vector{globalIndex};
            eigen.eigen_value{ind} =  IsotropicTensor.eigen_value{globalIndex};            
        else %% APPLY the transformation from DTI tensors to conductivity tensors
            %disp('case 4 : anisotrop ')
            if L2 == 0
                L2 = L1/largest_ratio;
            end
            if L3 == 0
                L3 = L2;
            end
            
            %% METHOD 1 : direct transformation approach [Güllmar et al   NeuroImage 2010]
            if aniso_method == 1  % Use the direct approach
                % Tuch parameters use - unstead of +
                k = 0.844 ;% +/- 0.0545;
                de = 0.124;
                lm1 = k*(L1 - de);
                lm2 = k*(L2- de);
                lm3 = k*(L3 - de);
                % The output
                aniso_conductivity(:,:,ind) = V * diag([lm1, lm2,lm3]) * V';
                eigen.eigen_vector{ind} = V;
                eigen.eigen_value{ind} = diag([lm1 lm2 lm3]);
                
                 meanConductivity = [meanConductivity mean([lm1 lm2 lm3]) ];
            end
            %% METHOD 2 : direct transformation approcah with volume constraint  [Güllmar et al NeuroImage 2010]
            if aniso_method == 2
                % apply the volume approcah
                % Tuch parameters
                k = 0.844 ;% +/- 0.0545;
                de = 0.124;
                lm1 = k*(L1 - de);
                lm2 = k*(L2 - de);
                lm3 = k*(L3 - de);
                % Apply the normalized volume
                lm1n = iso_conductivity * lm1/((lm1*lm2*lm3)^(1/3));
                lm2n = iso_conductivity * lm2/((lm1*lm2*lm3)^(1/3));
                lm3n = iso_conductivity * lm3/((lm1*lm2*lm3)^(1/3));
                % Check
                if ( lm1*lm2*lm3) == (lm1n*lm2n*lm3n)
                    disp('Volume is conserved');
                end
                % The output
                aniso_conductivity(:,:,ind) = V* diag([lm1n, lm2n,lm3n])*V';
                eigen.eigen_vector{ind} = V;
                eigen.eigen_value{ind} = diag([lm1n lm2n lm3n]);                
                meanConductivity = [meanConductivity mean([lm1 lm2 lm3])];
            end            
            %% METHOD 3 : Artificial anisotropy with volume constraint [Güllmar et al NeuroImage 2010]
            if aniso_method == 4
                % compute the new eigen value
                lm1 = (iso_conductivity^3 * ratio^2)^(1/3); % sigma longitidunal
                lm2 = (iso_conductivity^3 / ratio)^(1/3);      % sigma transversal
                if lm1 < lm2
                    temp = lm1;
                    lm1 = lm2;
                    lm2 = temp;
                end
                lm3 = lm2;
                % Check
                if (lm1*lm2*lm3) == (iso_conductivity^3)
                    disp('Volume is conserved');
                end
                % The output
                aniso_conductivity(:,:,ind) =  V* diag([lm1, lm2,lm3])*V';
                eigen.eigen_vector{ind} = V;
                eigen.eigen_value{ind} = diag([lm1 lm2 lm3]);                
                meanConductivity = [meanConductivity mean([lm1 lm2 lm3]) ];
            end
            %% METHOD 4 : The volume normalized approcah  [Won Hee Lee et all  IEEE 2015]
            if aniso_method == 4
                % apply the volume approcah
                lm1 = iso_conductivity * [(L1*L2*L3)^(-1/3)]*L1;
                lm2 = iso_conductivity * [(L1*L2*L3)^(-1/3)]*L2;
                lm3 = iso_conductivity * [(L1*L2*L3)^(-1/3)]*L3;
                if lm1 < lm2
                    temp = lm1;
                    lm1 = lm2;
                    lm2 = temp;
                end
                aniso_conductivity(:,:,ind) =  V * diag([lm1, lm2,lm3]) * V';
                eigen.eigen_vector{ind} = V;
                eigen.eigen_value{ind} = diag([lm1 lm2 lm3]);                
                meanConductivity = [meanConductivity mean([lm1 lm2 lm3])];
            end
            
            %% METHOD 5 : Wang-constraint model and volume-constraint model  [Wolters CH 2006]
            if aniso_method == 5
                % Wang constrained
                % sigma_rad * sigma_tang = sigma_tissu^2
                % sigma_rad : sigma_tang = ratio
                % sigma_tang = sqrt(sigma_tissu^2/r)
                % sigma_rad = sqrt(r*sigma_tissu^2)
                %     ratio = longitudinal_factor/transversal_factor;
                sigma_tissu = iso_conductivity;
                sigma_rad = sqrt((sigma_tissu^2)*ratio);
                sigma_tang = sqrt((sigma_tissu^2)/ratio);
                if sigma_rad > sigma_tang
                    lm1 = sigma_rad ;
                    lm2 =sigma_tang ;
                else
                    lm1 =sigma_tang;
                    lm2 = sigma_rad;
                end
                lm3 = lm2;
                aniso_conductivity(:,:,ind) =  V* diag([lm1, lm2 , lm3])*V';
                eigen.eigen_vector{ind} = V;
                eigen.eigen_value{ind} = diag([lm1 lm2 lm3]);                
                meanConductivity = [meanConductivity mean([lm1 lm2 lm3]) ];
            end
            
            %% METHOD 6 : Wang-constraint model and volume-constraint model  [Wolters CH 2006]
            % Noe : seems to be similar to the method 5
            if aniso_method == 6  % Wolters constrained
                % 4/3*pi*sigma_rad * (sigma_tang^2) = 4/3*pi*sigma_tissu^3
                % sigma_rad : sigma_tang = ratio
                % sigma_tang = ()
                % sigma_rad =
                % ratio = longitudinal_factor/transversal_factor;
                sigma_tissu = iso_conductivity;
                sigma_rad = ((sigma_tissu^3)*ratio^2)^(1/3);
                sigma_tang = ((sigma_tissu^3)/ratio)^(1/3);
                if sigma_rad > sigma_tang
                    lm1 = sigma_rad ;
                    lm2 =sigma_tang ;
                else
                    lm1 =sigma_tang;
                    lm2 = sigma_rad;
                end
                lm3 = lm2;
                aniso_conductivity(:,:,ind) =  V* diag([lm1, lm2 , lm3])*V';
                eigen.eigen_vector{ind} = V;
                eigen.eigen_value{ind} = diag([lm1 lm2 lm3]);                
                meanConductivity = [meanConductivity mean([lm1 lm2 lm3]) ];
            end
            
            %% METHOD 7 : As SimBio toolbox  [Rullmann et al 2008 / Vorwerk et al 2014 ]
            %% PART 1 : compute the scaling factor 's'
            if aniso_method == 7 % SimBio implementation
                % initialization
                if ind == 1
                    sumProductEigenValue = 0;
                    counter = 0;
                    g_index = [];
                    l_index = [];
                end
                % compute dcomp
                productEigenValue = L1 * L2 * L3;
                sumProductEigenValue = sumProductEigenValue + productEigenValue;
                counter = counter + 1;
                g_index = [g_index,  globalIndex]; % indexes for the whole head,global index
                l_index = [l_index,  ind]; % index for only the isotropic tissue or local index
                
                aniso_conductivity(:,:,ind) = nan(3,3);
                eigen.eigen_vector{ind} = nan;
                eigen.eigen_value{ind} = nan;
            end
            
        end
        anisotropicTensor.eigen_vector{globalIndex} = eigen.eigen_vector{ind} ;
        anisotropicTensor.eigen_value{globalIndex} = eigen.eigen_value{ind};
        ind = ind +1;
    else
        anisotropicTensor.eigen_vector{globalIndex} = IsotropicTensor.eigen_vector{globalIndex}  ;
        anisotropicTensor.eigen_value{globalIndex} = IsotropicTensor.eigen_value{globalIndex} ;
    end
end

%% METHOD 7 : Assign the the tensors
%%% PART 2 : Assigne the tensors
if aniso_method == 7
    meanDiffusity =  (sumProductEigenValue / counter)^(1/3);
    scalingFactor = iso_conductivity / meanDiffusity;        
    %% 1 - Apply the Tuch process Sigm = sD
    isNan = [];
    meanConductivity =[];
    for gind = l_index % the last value of l_index is not belonging to anistropic tissue
        if femHead.Tissue(gind) == aniso_tissu_index % 1 is the label of the wm
            try
                L = diag([L1a(gind), L2a(gind), L3a(gind)]);
                V = [V1rot(gind,:)', V2rot(gind,:)', V3rot(gind,:)'];
            catch ME
                disp(ME)
            end
        else
            continue
        end
        
        lm1 = scalingFactor * L(1,1);
        lm2 = scalingFactor * L(2,2);
        lm3 = scalingFactor * L(3,3);
        % rearange just in case
        if lm1 < lm2
            temp = lm1;
            lm1 = lm2;
            lm2 = temp;
        end
        
        %% 2 - Apply the normalized volume
        lm1n = iso_conductivity * (lm1/((lm1*lm2*lm3)^(1/3)));
        lm2n = iso_conductivity * (lm2/((lm1*lm2*lm3)^(1/3)));
        lm3n = iso_conductivity * (lm3/((lm1*lm2*lm3)^(1/3)));        
        % Check
        if (iso_conductivity*iso_conductivity*iso_conductivity) == (lm1n*lm2n*lm3n)
            disp('Volume is conserved');
        end
        useNormalizedVolume = 1;
        if useNormalizedVolume
            lm1 = lm1n;
            lm2 = lm2n;
            lm3 = lm3n;
        end        
        % check : the max value of the conductivity should not be larger
        % thant the SCF or the reference value
        % a maximal conductivity of CSF
        if lm1>max_conductivity; lm1 = max_conductivity; end
        if lm2>max_conductivity; lm2 = max_conductivity; end
        if lm3>max_conductivity; lm3 = max_conductivity; end
        
        meanConductivity = [meanConductivity mean([lm1 lm2 lm3]) ];
        isNan = [ isNan isnan(aniso_conductivity(1,1,gind)) ];
        
        % final output
        aniso_conductivity(:,:,gind) =  V* diag([lm1, lm2 , lm3])*V';
        anisotropicTensor.eigen_vector{gind} = V;
        anisotropicTensor.eigen_value{gind} = diag([lm1 lm2 lm3]);
    end
    % check the values for debugging
    [iso_conductivity max(meanConductivity)            mean(meanConductivity)    min(meanConductivity)]    ;      
    %% Scale the isotropic WM(others) tensors within the anisotropic tissue by the mean of the computed value 
    %%from the Wolter approach
    for ind = isIsotrope     
            aniso_conductivity(:,:,ind) = diag([mean(meanConductivity),mean(meanConductivity),mean(meanConductivity)]);
            eigen.eigen_value{ind} =  aniso_conductivity(:,:,ind) ;
    end    
end


end