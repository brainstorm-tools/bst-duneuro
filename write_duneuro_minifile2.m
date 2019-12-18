function write_duneuro_minifile2(cfg, minifile_name)
%  write_duneuro_minifile2(cfg, minifile_name)
%  
% This funcition write in a text file with all the confuguration (initialisation) parameters
% that will be used by the Duneuro application to compute the FEM solution
% The input is a structure that contains the parameters and the output is a
% file with *.mini extention.
% The name of the saved file is  output_filename if specified,
% if not, this function will use the name included in the
% minifile_data.name.
% ex : output_filename = 'myMiniFile.mini'
% This function is limited for the CG at this time
% Dependencies : nothing
%
% August 28, 2019, file created by Takfarinas MEDANI
%             updates : add the MEG configurations, december, 8th, 2019,
%             updates : add the MEEG option, for combined options , december, 13th,  2019,

% Add MEG mainly (to validate before merge to write_duneuro_minifile)
% TODO : add the unfitted method and its options
% refer to this : https://docs.google.com/spreadsheets/d/1MqURQsszn8Qj3-XRX_Z8qFFnz6Yl2-uYALkV-8pJVaM/edit#gid=0

%% Check the input arguments
if ~isstruct(cfg)
    error(['The format of %s is not correct, it should be a structure... check the help' ], inputname(1))
end

if nargin == 1
    minifile_name = cfg.mini_filename;
end

% Check the solver type :
if strcmp(cfg.minifile.solver_type,'udg')
    error (['solver_type : %s is not suppored for now ... we work on it', minifile.solver_type ])
end

% check if the extention is included 
[~,~,ext] = fileparts(minifile_name);
if isempty(ext)
    minifile_name = [minifile_name '.mini'];
end

% Open the mini file
fid = fopen(minifile_name, 'wt+');

%% 0 - Subpart general setting
fprintf(fid, '__name = %s\n\n',cfg.minifile.name);
if strcmp(cfg.minifile.solver_type,'cg')
    fprintf(fid, 'type = %s\n',cfg.minifile.type);
end
fprintf(fid, 'element_type  =%s\n',cfg.minifile.element_type );
fprintf(fid, 'solver_type = %s\n',cfg.minifile.solver_type);
fprintf(fid, 'geometry_adapted  = %s\n',cfg.minifile.geometry_adapted);
fprintf(fid, 'tolerance = %d\n',cfg.minifile.tolerance);



%% 1 - Subpart sensors
if strcmp(cfg.modality,'eeg') || strcmp(cfg.modality,'meeg') 
% subpart electrode : [electrodes]
fprintf(fid, '[electrodes]\n');
fprintf(fid, 'filename  = %s\n',cfg.minifile.electrode.filename);
fprintf(fid, 'type = %s\n',cfg.minifile.electrode.type);
end

if strcmp(cfg.modality,'meg') || strcmp(cfg.modality,'meeg') 
% subpart electrode : [meg]
fprintf(fid, '[meg]\n');
fprintf(fid, 'intorderadd  = %d\n', cfg.minifile.meg.intorderadd); % =0
fprintf(fid, 'type  = %s\n', cfg.minifile.meg.type); % = physical
% subpart coils : [coils]
fprintf(fid, '[coils]\n');
fprintf(fid, 'filename  = %s\n',cfg.coil_filename);
% subpart [projections]
fprintf(fid, '[projections]\n');
fprintf(fid, 'filename  = %s\n',cfg.minifile.projection_filename); % = text file
end

%% 2 -  subpart dipoles : [dipoles]
fprintf(fid, '[dipoles]\n');
fprintf(fid, 'filename  = %s\n',cfg.minifile.dipole.filename);

%% 3 - Subpart [volume_conductor.grid]
fprintf(fid, '[volume_conductor.grid]\n');
fprintf(fid, 'filename  = %s\n',cfg.minifile.volume_conductor_grid.filename);
% subpart  [volume_conductor.tensors]
fprintf(fid, '[volume_conductor.tensors]\n');
fprintf(fid, 'filename  = %s\n',cfg.minifile.volume_conductor_tensors.filename);

%% 4 - Subpart  [solver]
fprintf(fid, '[solver]\n');
fprintf(fid, 'solver_type  = %s\n',cfg.minifile.solver.solver_type); % cg pour conjugate gradient et pas continious galerkin
fprintf(fid, 'preconditioner_type  = %s\n',cfg.minifile.solver.preconditioner_type);
if strcmp(cfg.minifile.solver_type,'cg')
fprintf(fid, 'cg_smoother_type  = %s\n',cfg.minifile.solver.cg_smoother_type);
end
fprintf(fid, 'intorderadd  = %d\n',cfg.minifile.solver.intorderadd);
% case of the dg discontinious galerkin
if strcmp(cfg.minifile.solver_type,'dg')
fprintf(fid, 'dg_smoother_type  = %s\n',cfg.minifile.solver.dg_smoother_type);
fprintf(fid, 'scheme  = %s\n',cfg.minifile.solver.scheme);%/sipg
fprintf(fid, 'penalty  = %d\n',cfg.minifile.solver.penalty);%20
fprintf(fid, 'edge_norm_type  = %s\n',cfg.minifile.solver.penalty);%houston
fprintf(fid, 'weights  = %s\n',cfg.minifile.solver.penalty);%true
fprintf(fid, 'reduction  = %s\n',cfg.minifile.solver.reduction);%true
end

%% 5 - Subpart  [solution]
fprintf(fid, '[solution]\n');
fprintf(fid, 'post_process  = %s\n',cfg.minifile.solution.post_process); % true/false
fprintf(fid, 'subtract_mean  = %s\n',cfg.minifile.solution.subtract_mean); % boolean
% subpart  [solution.solver]
fprintf(fid, '[solution.solver]\n');
fprintf(fid, 'reduction  = %d\n',cfg.minifile.solution.solver.reduction);

%% subpart  [solution.source_model]
fprintf(fid, '[solution.source_model]\n');
fprintf(fid, 'type  = %s\n',cfg.minifile.solution.source_model.type );
fprintf(fid, 'intorderadd  = %d\n',cfg.minifile.solution.source_model.intorderadd);
fprintf(fid, 'intorderadd_lb  = %d\n',cfg.minifile.solution.source_model.intorderadd_lb);
fprintf(fid, 'numberOfMoments  = %d\n',cfg.minifile.solution.source_model.numberOfMoments);
fprintf(fid, 'referenceLength  = %d\n',cfg.minifile.solution.source_model.referenceLength);
fprintf(fid, 'weightingExponent  = %d\n',cfg.minifile.solution.source_model.weightingExponent );
fprintf(fid, 'relaxationFactor  = %d\n',cfg.minifile.solution.source_model.relaxationFactor);
fprintf(fid, 'mixedMoments  = %s\n',cfg.minifile.solution.source_model.mixedMoments);
fprintf(fid, 'restrict  = %s\n',cfg.minifile.solution.source_model.restrict );
fprintf(fid, 'initialization  = %s\n',cfg.minifile.solution.source_model.initialization);

% The reste is not needed... just in case for further use
% subpart [analytic_solution]
fprintf(fid, '[analytic_solution]\n');
Nb_Layer = length(cfg.minifile.solution.analytic_solution.radii);
fprintf(fid, 'radii =\t');
for ind = 1 : Nb_Layer
    fprintf(fid, '%d \t', cfg.minifile.solution.analytic_solution.radii(ind));
end
fprintf(fid, '\n');
fprintf(fid, 'center =\t');
for ind = 1 : 3
    fprintf(fid, '%d \t', cfg.minifile.solution.analytic_solution.center(ind));
end
fprintf(fid, '\n');
fprintf(fid, 'conductivities =\t');
for ind = 1 : Nb_Layer
    fprintf(fid, '%d \t', cfg.minifile.solution.analytic_solution.conductivities(ind));
end
fprintf(fid, '\n');
% subpart  [output]
fprintf(fid, '[output]\n');
fprintf(fid, 'filename  = %s\n',cfg.minifile.output.filename);
fprintf(fid, 'extension  = %s\n',cfg.minifile.output.extension);
% subpart [wrapper.outputtreecompare]
fprintf(fid, '[wrapper.outputtreecompare]\n');
fprintf(fid, 'name  = %s\n',cfg.minifile.wrapper.outputtreecompare.name);
fprintf(fid, 'extension  = %s\n',cfg.minifile.wrapper.outputtreecompare.extension);
fprintf(fid, 'reference  = %s\n',cfg.minifile.wrapper.outputtreecompare.reference);
fprintf(fid, 'type  = %s\n',cfg.minifile.wrapper.outputtreecompare.type);
fprintf(fid, 'absolute  = %s\n',cfg.minifile.wrapper.outputtreecompare.absolute);
fclose(fid);

end