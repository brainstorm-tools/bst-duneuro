function [iNewSurfaces, OutputFiles] = fem_hexa2tetra(iSubject, FemFiles, FileFormat, isInteractive)
% IMPORT_FEMLAYERS: Extracts surfaces from FEM 3D mesh and saves them in the database
%
% USAGE: iNewSurfaces = import_surfaces(iSubject, FemFiles, FileFormat)
%        iNewSurfaces = import_surfaces(iSubject)   : Ask user the files to import
%
% INPUT:
%    - iSubject     : Indice of the subject where to import the surfaces
%                     If iSubject=0 : import surfaces in default subject
%    - FemFiles     : Cell array of full filenames of the surfaces to import (format is autodetected)
%                     => if not specified : files to import are asked to the user
%    - FileFormat   : String representing the file format to import.
%                     Please see in_tess.m to get the list of supported file formats
%    - isInteractive: {0,1} If 0, do not ask any question to the user and use default values
% OUTPUT:
%    - iNewSurfaces : Indices of the surfaces added in database

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c)2000-2020 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Francois Tadel, 2020


%% ===== PARSE INPUTS =====
% Check command line
if ~isnumeric(iSubject) || (iSubject < 0)
    error('Invalid subject indice.');
end
if (nargin < 4) || isempty(isInteractive)
    isInteractive = 0;
end
if (nargin < 3) || isempty(FemFiles)
    FemFiles = {};
    FileFormat = [];
else
    if ischar(FemFiles)
        FemFiles = {FemFiles};
    end
    if (nargin == 2) || ((nargin >= 3) && isempty(FileFormat))
        error('When you pass a FemFiles argument, FileFormat must be defined too.');
    end
end
iNewSurfaces = [];
OutputFiles = {};
nVertices = [];

% Get Protocol information
ProtocolInfo = bst_get('ProtocolInfo');
% Get subject directory
sSubject = bst_get('Subject', iSubject);
subjectSubDir = bst_fileparts(sSubject.FileName);


%% ===== INSTALL ISO2MESH =====
% Install iso2mesh if needed
if ~exist('iso2meshver', 'file') || ~isdir(bst_fullfile(bst_fileparts(which('iso2meshver')), 'doc'))
    errMsg = process_generate_fem('InstallIso2mesh', isInteractive);
    if ~isempty(errMsg) || ~exist('iso2meshver', 'file') || ~isdir(bst_fullfile(bst_fileparts(which('iso2meshver')), 'doc'))
        warning('Could not find Iso2mesh on your computer... the extracted surface may have some isolated faces.')
    end
end


%% ===== SELECT INPUT FILES =====
% If surface files to load are not defined : open a dialog box to select it
if isempty(FemFiles)
    % Get last used directories and formats
    LastUsedDirs = bst_get('LastUsedDirs');
    % Get Surface files
    [FemFiles, FileFormat, FileFilter] = java_getfile( 'open', ...
        'Import surfaces...', ...     % Window title
        LastUsedDirs.ImportAnat, ...   % Default directory
        'multiple', 'files', ...      % Selection mode
        {{'_fem'}, 'Brainstorm (*.mat)', 'BSTFEM'}, ...
        'BSTFEM');
    % If no file was selected: exit
    if isempty(FemFiles)
        return
    end
    % Save default import directory
    LastUsedDirs.ImportAnat = bst_fileparts(FemFiles{1});
    bst_set('LastUsedDirs', LastUsedDirs);
end


%% ===== LOAD EACH SURFACE =====
% Process all the selected surfaces
for iFile = 1:length(FemFiles)
    % Load file
    FemFile = FemFiles{iFile};
    bst_progress('start', 'Extract surfaces', ['Loading file "' FemFile '"...']);
    FemMat = load(FemFile);   
    
    if size(FemMat.Elements,2) == 4
        warning('This mesh is already on TETRA')
        bst_progress('stop');
        return;
        %         disp('Convert from Tetra 2 Hexa')
        %         % Subdeviding the hexahedral element
        %         [Es,Vs] = tet2hex(FemMat.Elements,FemMat.Vertices);
        %         hexaLabel = repmat(FemMat.Tissue,1,4); hexaLabel = hexaLabel';
        %         hexaLabel= hexaLabel(:);
        %         FemMat.Vertices = Vs;
        %         FemMat.Tissue =    hexaLabel;
        %         FemMat.Elements = Es;
        %         % Create output structure
        %         FemMat.Comment = sprintf('FEM %dV (%s, %d layers)',length(FemMat.Vertices),'Hexa mesh : Tetra2Hexa', length(unique( FemMat.Tissue)) );
        %         % Add history
        %         FemMat = bst_history('add', FemMat, 'fem_tetra2hexa');
    elseif size(FemMat.Elements,2) == 8
         disp('Convert from Hexa 2 Tetra ')
        % convert the mesh to tetra for diplay purpose
        [tetraElem,tetraNode,tetraLabel]=hex2tet(FemMat.Elements,FemMat.Vertices ,FemMat.Tissue,4);
        % updates FemMat for display purpose
        FemMat.Vertices =  tetraNode;
        FemMat.Elements = tetraElem;
        FemMat.Tissue = tetraLabel;
        % Create output structure
        FemMat.Comment = sprintf('FEM %dV (%s, %d layers)',length(FemMat.Vertices),'Tetra mesh : Hexa2Tetra ', length(unique( FemMat.Tissue)) );
        % Add history
        FemMat = bst_history('add', FemMat, 'fem_hexa2tetra');
    end    
    
    % ===== SAVE FEM MESH =====
    bst_progress('text', 'Saving FEM mesh...');    
    % Save to database
    % Produce a default surface filename
    BstFemFile = bst_fullfile(ProtocolInfo.SUBJECTS, subjectSubDir, 'tess_fem.mat');
    % Make this filename unique
    BstFemFile = file_unique(BstFemFile);
    % Save new surface in Brainstorm format
    bst_save(BstFemFile, FemMat, 'v7');
    db_add_surface(iSubject, BstFemFile, FemMat.Comment)
    % db_add_surface(iSubject, FemFile, FemMat.Comment);    
    % Return success
    isOk = 1;
    bst_progress('stop');
end
end

