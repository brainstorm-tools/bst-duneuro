function write_duneuro_projection_file(coil_orientation,coilOrientation_filename)
%  write_duneuro_projection_file(coil_orientation,coilOrientation_filename)

% Authors: Takfarinas MEDANI, Decembre 2019;     

[filepath,name,ext] = fileparts(coilOrientation_filename);
if isempty(ext) || ~strcmp(ext,'.txt')
    ext = '.txt';
end
coilOrientation_filename = fullfile(filepath,[name,ext]);
% Nb_coils = size(coil_orientation, 1); 
% generate triedre orientation for each dipole
% coil_pos_orie = [kron(coil_pos,ones(3,1)), kron(ones(Nb_dipole,1), eye(3))];
% coil_orie = [ kron(ones(Nb_coils,1), eye(3))];
% coil_orie = repmat(([1 0 0 0 1 0 0 0 1]),Nb_coils,1);
fid = fopen(coilOrientation_filename, 'wt+');
% fprintf(fid, '%d %d %d %d %d %d \n', coil_pos_orie');
fprintf(fid, '%d %d %d  \n', coil_orientation');
fclose(fid); 
end