function varargout = panel_duneuro(varargin)
% PANEL_OPENMEEG: Options for OpenMEEG BEM (GUI).
% 
% USAGE:  bstPanelNew = panel_openmeeg('CreatePanel', OPTIONS)           : Call from the interactive interface
%         bstPanelNew = panel_openmeeg('CreatePanel', sProcess, sFiles)  : Call from the process editor
%                   s = panel_openmeeg('GetPanelContents')

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
% Authors: Francois Tadel, 2011-2019

eval(macro_method);
end

%% ===== CREATE PANEL =====
function [bstPanelNew, panelName] = CreatePanel(sProcess, sFiles)  %#ok<DEFNU>  
    panelName = 'DuneuroOptions';
    % Java initializations
    import java.awt.*;
    import javax.swing.*;
    % CALL: From GUI
    if (nargin == 1)
        OPTIONS = sProcess;
    % CALL: From Process
    else
%         OPTIONS = sProcess.options.openmeeg.Value;
    end
    % Default options
    if isempty(OPTIONS)
        OPTIONS = struct();
    end
    defOPTIONS = bst_get('DuneuroMEEGOptions');
    OPTIONS = struct_copy_fields(OPTIONS, defOPTIONS, 0);
    OPTIONS.FemCond = OPTIONS.Conductivity;
    % Create main panel
    jPanelNew = gui_river();
    
    % ===== FEM LAYERS =====
    jPanelLayers = gui_river([4,4], [3,15,10,10], 'FEM Layers & Conductivities');
        nLayers = length(OPTIONS.FemNames);
        jCheckLayer = javaArray('javax.swing.JCheckBox', nLayers);
        jLabelLayer = javaArray('javax.swing.JLabel', nLayers);
        jTextCond   = javaArray('javax.swing.JTextField', nLayers);
        % Loop on each layer
        for i = 1:nLayers
            % Add components
            jCheckLayer(i) = gui_component('checkbox', jPanelLayers, 'br',  OPTIONS.FemNames{i}, [], [], @UpdatePanel, []);
%             jLabelLayer(i) = gui_component('label', jPanelLayers, 'tab', strVert, [], [], [], []);
            jTextCond(i) = gui_component('texttime', jPanelLayers, 'tab', num2str(OPTIONS.FemCond(i), '%g'), [], [], [], []);
            % EEG: Select all layers; MEG: Select only the innermost layer
            jCheckLayer(i).setSelected(OPTIONS.FemSelect(i));
        end
    jPanelNew.add('br hfill', jPanelLayers);

    % ===== DUNEURO OPTIONS ======
    jPanelDuneuro = gui_river([3,3], [3,15,10,10], 'DuneuroMEEG options');
        % Venant
        jCheckVenant = gui_component('radio', jPanelDuneuro, 'br', '<HTML> Venant <FONT COLOR="#808080"></FONT>', [], [], @UpdatePanel, []);
        jCheckVenant.setSelected(OPTIONS.isVenant);
        % Subtraction
        jCheckSubtraction = gui_component('radio', jPanelDuneuro, 'br', '<HTML> Subtraction <FONT COLOR="#808080"></FONT>', [], [], @UpdatePanel, []);
        jCheckSubtraction.setSelected(OPTIONS.isSubtraction);
        % PartialIntegration
        jCheckPartialIntegration = gui_component('radio', jPanelDuneuro, 'br',  '<HTML> Partial Integration  <FONT COLOR="#808080"></FONT>', [], [], @UpdatePanel, []);
        jCheckPartialIntegration.setSelected(OPTIONS.isPartialIntegration);
    jPanelNew.add('br hfill', jPanelDuneuro);
        
    % ===== VALIDATION BUTTONS =====
    gui_component('button', jPanelNew, 'br right', 'Cancel', [], [], @ButtonCancel_Callback, []);
    gui_component('button', jPanelNew, [], 'OK', [], [], @ButtonOk_Callback, []);

    % ===== PANEL CREATION =====
    % Return a mutex to wait for panel close
    bst_mutex('create', panelName);
    % Controls list
    ctrl = struct('jCheckLayer',      jCheckLayer, ...
                  'jTextCond',        jTextCond, ...
                  'jCheckVenant',    jCheckVenant, ...
                  'jCheckSubtraction', jCheckSubtraction, ...
                  'jCheckPartialIntegration',      jCheckPartialIntegration);
    ctrl.FemFiles = OPTIONS.FemFiles;
    ctrl.FemCond  = OPTIONS.FemCond;
    ctrl.FemNames = OPTIONS.FemNames;
    % Create the BstPanel object that is returned by the function
    % => constructor BstPanel(jHandle, panelName, sControls)
    bstPanelNew = BstPanel(panelName, jPanelNew, ctrl);
    % Update panel
    UpdatePanel();
    

%% =================================================================================
%  === INTERNAL CALLBACKS ==========================================================
%  =================================================================================
%% ===== CANCEL BUTTON =====
    function ButtonCancel_Callback(hObject, event)
        % Close panel without saving (release mutex automatically)
        gui_hide(panelName);
        bst_mutex('release', panelName);
    end

%% ===== OK BUTTON =====
    function ButtonOk_Callback(varargin)       
        % Release mutex and keep the panel opened
        bst_mutex('release', panelName);
    end

%% ===== UPDATE PANEL =====
 function UpdatePanel(varargin)
        % Venant
        isVenant = jCheckVenant. isSelected();
        if isVenant
             jCheckPartialIntegration.setEnabled(1);
             jCheckSubtraction.setEnabled(1);
        end
        % PartialIntegration
        isPartialIntegration = jCheckPartialIntegration.isSelected();
        if isPartialIntegration
            jCheckVenant.setEnabled(1);
            jCheckSubtraction.setEnabled(1);
        end        
        % Subtraction
        isSubtraction = jCheckSubtraction.isSelected();
        if isSubtraction
            jCheckVenant.setEnabled(1);
            jCheckPartialIntegration.setEnabled(1);
        end

%         jCheckVenant.setSelected(OPTIONS.isVenant);
%         jCheckVenant.setEnabled();
end

%% =================================================================================
%  === EXTERNAL CALLBACKS ==========================================================
%  =================================================================================   
%% ===== GET PANEL CONTENTS =====
function s = GetPanelContents() %#ok<DEFNU>
    % Get panel controls
    ctrl = bst_get('PanelControls', 'DuneuroMEEGOptions');
    % Fem layers
    for i = 1:length(ctrl.jCheckLayer)
        s.FemSelect(i) = ctrl.jCheckLayer(i).isSelected();
        s.FemCond(i) = str2num(char(ctrl.jTextCond(i).getText()));
    end
    s.FemNames = ctrl.FemNames;
    s.FemFiles = ctrl.FemFiles;
    % DuneuroMEEG Options 
    s.isVenant    = ctrl.jCheckVenant.isSelected();
    s.isSubtraction = ctrl.jCheckSubtraction.isSelected();
    s.isPartialIntegration      = ctrl.jCheckPartialIntegration.isSelected();
end

end