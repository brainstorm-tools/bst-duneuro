function [conductivity, tissuLabel]  = get_standard_conductivity(numberOfLayer,reference)

% use :
% [conductivity, tissuLabel] = bst_standard_conductivity(numberOfLayer, reference)
%                                             return the default values of
%                                             the conductivity and the name
%                                             of the layers
% define the std conductivity
% refers to standard_cond
% from inner to outer ... always....


if nargin == 1
    reference =1; % by default use simbio value
end
% reference :
% simNibsValue = 0; % use value published by SimNibs team, Axel
% simBioValue = 1; % use value published by SimBio team, Carsten and Johannes : A guidline for head volume ....


if reference == 0 %simNibsValue = 0
    switch numberOfLayer
        case 6
            tissu = {'WM','GM','CSF','Skull Spongia','Skull Compacta','Scalp'};% simnibs paper & soft
            conductivity = [ 0.126 0.275 1.654 0.010 0.465];% simnibs paper & soft
            
        case 5
            tissu = {'WM','GM','CSF','Skull','Scalp'};% simnibs paper & soft
            conductivity = [ 0.126 0.275 1.654 0.010 0.465];% simnibs paper & soft
            
        case 4
            conductivity = [0.275 1.654 0.010 0.465];% simnibs paper & soft
            tissu = {'GM','CSF','Skull','Scalp'};% simnibs paper & soft
        case 3
            conductivity = [0.275 0.010 0.465];% simnibs paper & soft
            tissu = {'GM','Skull','Scalp'};% simnibs paper & soft
        case 2
            conductivity = [0.275 0.465];% simnibs paper & soft
            tissu = {'GM','Scalp'};% simnibs paper & soft
        case 1
            conductivity = [ 0.465];% simnibs paper & soft
            tissu = {'Scalp'};% simnibs paper & soft
        otherwise
            error('error on your model')
    end
end

if reference == 1 %simBioValue = 0
    switch numberOfLayer
        case 6
            tissu = {'WM','GM','CSF','Skull Spongia','Skull Compacta','Scalp'}; %SimBio team
            conductivity = [ 0.14 0.33 1.79 0.025 0.008 0.43];
            
        case 5
            tissu = {'WM','GM','CSF','Skull','Scalp'};
            conductivity = [ 0.14 0.33 1.79  0.01 0.43];
            
        case 4
            conductivity = [ 0.33 1.79 0.01 0.43];
            tissu = {'GM','CSF','Skull','Scalp'};
        case 3
            conductivity = [0.33 0.01 0.43];
            tissu = {'GM','Skull','Scalp'};
        case 2
            conductivity = [0.33 0.43];
            tissu = {'GM','Scalp'};
        case 1
            conductivity = [0.43];
            tissu = {'Scalp'};
        otherwise
            error('error on your model')
    end
end


tissuLabel = tissu;

end
