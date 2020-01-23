function cfg = bst_set_minifile(cfg)
%% Script to fill the vol.cfg.minifile with the duneuro parameters
% Some of the parameters will be classified as :
% 1 : Fixed parameters that the user can't see and change
% 2 : Advanced parameters that advanced used can change
% 3 : Basic parameters that user can tune and that we will set to default.
% This function should be tuned from outside via the GUI
% 
%Author :Takfarinas MEDANI, Creation  August 28, 2019,
%             updates : add the MEG configurations, december 209,

% Load the default value for the mnifile // duneuro configuration file 
% TODO : add the unfitted method and its options
% refer to this : https://docs.google.com/spreadsheets/d/1MqURQsszn8Qj3-XRX_Z8qFFnz6Yl2-uYALkV-8pJVaM/edit#gid=0

% Set the default parameters
cfg = bst_default_duneuroConfiguration(cfg);

%% UPDATES TO THE REAL VALUES
%% 0 - Subpart general setting
cfg.minifile.name = cfg.duneuro_configuration_filename;
cfg.minifile.type  = cfg.dnFemMethodType; % 'fitted' or 'unfitted'
cfg.minifile.solver_type  = cfg.dnFemSolverType ; % cg, dg, udg (ucg or cut to check), maybe the mini file should change according to solver_type.
cfg.minifile. element_type  = cfg.dnMeshElementType; %  'tetrahedron' or  'hexahedron'
cfg.minifile.geometry_adapted  =  cfg.dnGeometryAdapted ; % what does it mean ? 
cfg.minifile.tolerance  = cfg.dnTolerance;

%% 1 - Subpart sensors
% subpart electrode : [electrodes]
cfg.minifile.electrode.filename  = cfg.electrode_filename;
cfg.minifile.electrode.type   = cfg.dnElectrodType; % 'normal' what are the other options ?

%% subpart [meg]
cfg.minifile.meg.intorderadd = cfg.dnMegIntorderadd;
cfg.minifile.meg.type = cfg.dnMegType;
% subpart coils : [coils]
cfg.minifile.coil_filename = cfg.coil_filename;
% subpart [projection]
cfg.minifile.projection_filename = cfg.projection_filename;

%% 2 -  subpart dipoles : [dipoles]
cfg.minifile.dipole.filename = cfg.dipole_filename ;

%% 3 - Subpart [volume_conductor.grid]
cfg.minifile.volume_conductor_grid.filename = cfg.head_filename ;% subpart  [volume_conductor.tensors]
cfg.minifile.volume_conductor_tensors.filename = cfg.cond_filename;

%% 4 - Subpart  [solver]
cfg.minifile.solver.solver_type = cfg.dnSolverSolverType; %'cg';  what are the others 
cfg.minifile.solver.preconditioner_type = cfg.dnSolverPreconditionerType; % 'amg' 
cfg.minifile.solver.cg_smoother_type = cfg.dnSolverCgSmootherType; %  'ssor';
cfg.minifile.solver.intorderadd = cfg.dnSolverIntorderadd; %0;

% case of the dg discontinious galerkin
cfg.minifile.solver.dg_smoother_type = cfg.dnGgSmootherType;
cfg.minifile.solver.scheme = cfg.dnScheme;%/sipg
cfg.minifile.solver.penalty = cfg.dnPenalty;%20
cfg.minifile.solver.penalty = cfg.dnEdgeNormType;%houston
cfg.minifile.solver.penalty = cfg.dnWeights ;%true
cfg.minifile.solver.reduction = cfg.dnWeights;%true

%% 5 - Subpart  [solution]
cfg.minifile.solution.post_process = cfg.dnSolutionPostProcess; %  'true'; cfg.dnSolutionPostProcess =  'true'
cfg.minifile.solution.subtract_mean = cfg.dnSolutionSubstractMean; % 'true'; cfg.dnSolutionSubstractMean =  'true'
%subpart  [solution.solver]
cfg.minifile.solution.solver.reduction = cfg.dnSolutionSolverReduction; %  1e-10 ; cfg.dnSolutionSolverReduction = 1e-10;

%% 6 - subpart  [solution.source_model]
cfg.minifile.solution.source_model.type = cfg.femSourceModel;% partial_integration, venant, subtraction | expand smtype
cfg.minifile.solution.source_model.intorderadd =  cfg.femSourceModelIntorderadd ;%0;
cfg.minifile.solution.source_model.intorderadd_lb =  cfg.femSourceModelIntorderadd_lb;%2;
cfg.minifile.solution.source_model.numberOfMoments = cfg.femSourceModelNumberOfMoments;%3;
cfg.minifile.solution.source_model.referenceLength = cfg.femSourceModelReferenceLength;%20;
cfg.minifile.solution.source_model.weightingExponent = cfg.femSourceModelWeightingExponent;%1;
cfg.minifile.solution.source_model.relaxationFactor = cfg.femSourceModelRelaxationFactor; %1e-6;
cfg.minifile.solution.source_model.mixedMoments = cfg.femSourceModelMixedMoments;% 'false';
cfg.minifile.solution.source_model.restrict = cfg.femSourceModelRestrict;%'true';
cfg.minifile.solution.source_model.initialization = cfg.femSourceModelInitialization;%'closest_vertex';

% The reste is not needed... just in case
% subpart [analytic_solution]
cfg.minifile.solution.analytic_solution.radii = [1 2 3 4 ];
cfg.minifile.solution.analytic_solution.center = [0 0 0];
cfg.minifile.solution.analytic_solution.conductivities = [1 0.0125 1 1];
% subpart  [output]
cfg.minifile.output.filename = 'ns';
cfg.minifile.output.extension = 'ns';
% subpart [wrapper.outputtreecompare]
cfg.minifile.wrapper.outputtreecompare.name = 'ns';
cfg.minifile.wrapper.outputtreecompare.extension = 'ns';
cfg.minifile.wrapper.outputtreecompare.reference = 'ns';
cfg.minifile.wrapper.outputtreecompare.type = 'ns';
cfg.minifile.wrapper.outputtreecompare.absolute = 'ns'; %1e-2;
end