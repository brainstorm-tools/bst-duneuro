function cfg = bst_default_duneuroConfiguration(varargin)
% cfg = bst_load_default_duneuroConfiguration(cfg);
% This function will load or update the default duneuro values for the configuration file.
% If there is an argument, this function will complete the cfg otherwise it
% will create new structure with default parameters.

% Takfarinas MEDANI
%             updates : add the MEG configurations, december 209,

% TODO : add the unfitted method and its options
% refer to this : https://docs.google.com/spreadsheets/d/1MqURQsszn8Qj3-XRX_Z8qFFnz6Yl2-uYALkV-8pJVaM/edit#gid=0


if nargin == 0
    cfg = [];
else
    cfg = varargin{1};
end    
 
%% 0 - Subpart general setting
if ~isfield(cfg,'dnFemMethodType'); cfg.duneuro_configuration_filename = 'duneuro_minifile.mini'; end % 'fitted' or 'unfitted'
if ~isfield(cfg,'dnFemMethodType'); cfg.dnFemMethodType = 'fitted'; end % 'fitted' or 'unfitted'
if ~isfield(cfg,'dnFemSolverType'); cfg.dnFemSolverType = 'cg'; end % 'fitted' or 'unfitted'
if ~isfield(cfg,'dnMeshElementType'); cfg.dnMeshElementType = 'tetrahedron'; end %  'tetrahedron' or  'hexahedron'
if ~isfield(cfg,'dnGeometryAdapted'); cfg.dnGeometryAdapted = 'false'; end %  'true' or  'false'
if ~isfield(cfg,'dnTolerance'); cfg.dnTolerance = 1e-8; end %  

%% 1 - Subpart sensors
% subpart electrode : [electrodes]
if ~isfield(cfg,'electrode_filename'); cfg.electrode_filename = 'electrode_model.txt'; end %  
if ~isfield(cfg,'dnElectrodType'); cfg.dnElectrodType = 'normal'; end %  

%% subpart [meg]
if ~isfield(cfg,'dnMegIntorderadd');cfg.dnMegIntorderadd  = 0; end
if ~isfield(cfg,'dnMegType');cfg.dnMegType  = 'physical'; end
% subpart coils : [coils]
if ~isfield(cfg,'coil_filename'); cfg.coil_filename ='coil_model.txt'; end %  
% subpart [projection]
if ~isfield(cfg,'projection_filename');cfg.projection_filename  = 'projection_model.txt'; end

%% 2 -  subpart dipoles : [dipoles]
if ~isfield(cfg,'dipole_filename'); cfg.dipole_filename ='dipole_model.txt'; end %  what are the others 

%% 3 - Subpart [volume_conductor.grid]
if ~isfield(cfg,'head_filename'); cfg.head_filename ='head_model.msh'; end %  what are the others 
if ~isfield(cfg,'cond_filename'); cfg.dnSolverSolverType ='conductivity_model.con'; end %  what are the others 

%% 4 - Subpart  [solver] ==> refers to the linear system solver ?
if ~isfield(cfg,'dnSolverSolverType'); cfg.dnSolverSolverType ='cg'; end %  what are the others 
if ~isfield(cfg,'dnSolverPreconditionerType'); cfg.dnSolverPreconditionerType ='amg'; end %  what are the others 
if ~isfield(cfg,'dnSolverCgSmootherType'); cfg.dnSolverCgSmootherType ='ssor'; end %  what are the others 
if ~isfield(cfg,'dnSolverIntorderadd'); cfg.dnSolverIntorderadd =0; end %  what are the others 

% case of the dg discontinious galerkin
if ~isfield(cfg,'dnGgSmootherType'); cfg.dnGgSmootherType = 'ssor'; end %  what are the others 
if ~isfield(cfg,'dnScheme'); cfg.dnScheme = 'sipg'; end %  what are the others 
if ~isfield(cfg,'dnPenalty'); cfg.dnPenalty = 20; end %  what are the others 
if ~isfield(cfg,'dnEdgeNormType'); cfg.dnEdgeNormType = 'houston'; end %  what are the others 
if ~isfield(cfg,'dnWeights'); cfg.dnWeights = 'true'; end %  what are the others 
if ~isfield(cfg,'dnReduction'); cfg.dnWeights = 'true'; end %  what are the others 

%% 5 - Subpart  [solution]
if ~isfield(cfg,'dnSolutionPostProcess'); cfg.dnSolutionPostProcess ='true'; end %  what are the others 
if ~isfield(cfg,'dnSolutionSubstractMean'); cfg.dnSolutionSubstractMean ='false'; end %  what are the others 
% subpart  [solution.solver]
if ~isfield(cfg,'dnSolutionSolverReduction'); cfg.dnSolutionSolverReduction = 1e-10; end %  what are the others 

%% 6 - subpart  [solution.source_model]
if ~isfield(cfg,'femSourceModel'); cfg.femSourceModel = 'venant'; end % partial_integration, venant, subtraction | expand smtype
if ~isfield(cfg,'femSourceModelIntorderadd'); cfg.femSourceModelIntorderadd = 0; end % 
if ~isfield(cfg,'femSourceModelIntorderadd_lb'); cfg.femSourceModelIntorderadd_lb = 2; end % 
if ~isfield(cfg,'femSourceModelNumberOfMoments'); cfg.femSourceModelNumberOfMoments = 3; end % 
if ~isfield(cfg,'femSourceModelReferenceLength'); cfg.femSourceModelReferenceLength = 20; end % 
if ~isfield(cfg,'femSourceModelWeightingExponent'); cfg.femSourceModelWeightingExponent = 1; end % 
if ~isfield(cfg,'femSourceModelRelaxationFactor'); cfg.femSourceModelRelaxationFactor = 1e-6; end % 
if ~isfield(cfg,'femSourceModelMixedMoments'); cfg.femSourceModelMixedMoments = 'true'; end % 
if ~isfield(cfg,'femSourceModelRestrict'); cfg.femSourceModelRestrict = 'true'; end % 
if ~isfield(cfg,'femSourceModelInitialization'); cfg.femSourceModelInitialization = 'closest_vertex'; end % 

end