function sample_eig2nifti
% EIG2NIFTI convertes .eig.nii.gz files (as saved by BDP - BrainSuite
% Diffusion Pipeline) to standard NIfTI-1 format.
%
% 'src' directory must be added to MATLAB path before usage. Also, if
% applicable, remove other NIfTI readers from MATLAB path.
%
% Usage:
%      eig2nifti(input_eig_file, output_base)
%
% See http://neuroimage.usc.edu/neuro/Resources/eig2nifti for more details.
%

addpath('src'); % replace by correct path of 'src' directory

eig_filename = 'sub.eig.nii.gz'; % filename can be specified with full path
output_base = 'sub.tensor'; % file prefix string

eig2nifti(eig_filename, output_base);

% Above example will save following files: 
%  sub.tensor.V1.nii.gz - Eigenvector corresponding to largest eigenvalue
%  sub.tensor.V2.nii.gz - Eigenvector corresponding to second eigenvalue
%  sub.tensor.V3.nii.gz - Eigenvector corresponding to third eigenvalue
%  sub.tensor.L1.nii.gz - largest eigenvalue
%  sub.tensor.L2.nii.gz - second eigenvalue
%  sub.tensor.L3.nii.gz - third eigenvalue
