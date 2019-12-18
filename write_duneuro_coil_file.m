function write_duneuro_coil_file(coil_loc, coil_filename)
% write_duneuro_coil_file(coil_loc, coil_filename)
% Write the electrode file for Duneuro application
% coil_loc : 3D Cartesien position of the coils (or integration point), Ncoil x 3

% Authors: Takfarinas MEDANI, December 2019;     

[filepath,name,ext] = fileparts(coil_filename);
if isempty(ext) || ~strcmp(ext,'.txt')
    ext = '.txt';
end
coil_filename = [filepath,name,ext];

fid = fopen(coil_filename, 'wt+');
fprintf(fid, '%d %d %d  \n', coil_loc');
fclose(fid); 
end