
function process_dwi2conductivity(iSubject)

% Generate the diffusion tensor
[diffusion_tensor, errMsg] = bst_generate_diffusion_tensor(iSubject);
 
% Generate the conductivity tensor
[aniso_conductivity_tensor, eigenData, errMsg ] = bst_generate_conductivity_tensor(iSubject, DTI);


% Display the tensors
bst_display_tensors(iSubject,eigenData)

end