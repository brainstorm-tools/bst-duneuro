function ftMri = bstMri2ftMri(bstMri, vox2ras)
% basic convertion from the brainstorm MRI file to Filedtrip MRI file (compatibel with ft_read_mri)
% ftMri = bstMri2ftMri(bstMri) 
% ftMri = bstMri2ftMri(bstMri, vox2ras) 
% input bstMri : MRI struvture in the brainstorm format
%         vox2ras : transfomation 4x4 matrix 
% output  ftMri : MRI struvture in the fieltrip format

% Takfarinas MEDANI, 03/10/2020

if nargin == 1
    vox2ras = [];
end
% Re-write to the fieldtrip format
ftMri.dim = size(bstMri.Cube);
ftMri.anatomy = bstMri.Cube;
ftMri.hdr = bstMri.Header;
ftMri.transform = vox2ras;
ftMri.unit = 'mm'; % we assume that milimiters are used
end