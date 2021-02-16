function varargout = CellReg(varargin)
% This GUI is an implementation of a probabilistic approach for the
% identification of the same neurons (cell registration) across multiple sessions
% in Ca2+ imaging data, developed by Sheintuch et al., 2017.

% Input: The inputs for the cell registration method are the spatial footprints of
% cellular activity of the cells that were detected in the different
% sessions. Each spatial footprint is a matrix the size of the frame and
% each pixel's value represents its contribution to the
% cell's fluorescence.

% Output: The main output for the cell registration method is the obtained mapping of
% cell identity across all registered sessions. It is a matrix the size of
% the final number of registered cells by the number of registered
% sessions. Each entry holds the index for the cell in a given session.
% Other outputs include:
% 1) register scores - providing with the registration quality of each cell register
% 2) log file - with all the relevant information regarding the data, registration
% configurations, and a summary of the registration results and quality.
% 3) figures - important figures that are saved automatically.

% The GUI includes the following stages:
% 1) Loading the spatial footprints of cellular activity from the different sessions.

% 2) Transforming all the sessions to a reference coordinate system using
% rigid-body transformation.

% 3) Computing a probabilistic model of the spatial footprints similarities
% of neighboring cell-pairs from different sessions using the centroid
% distances and spatial correlations.

% 4) Obtaining an initial cell registration according to an optimized registration threshold.

% 5) Obtaining the final cell registration based on a correlation clustering algorithm.

% CellReg MATLAB code for CellReg.fig
%      CellReg, by itself, creates a new CellReg or raises the existing
%      singleton*.
%
%      H = CellReg returns the handle to a new CellReg or the handle to
%      the existing singleton*.
%
%      CellReg('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CellReg.M with the given input arguments.
%
%      CellReg('Property','Value',...) creates a new CellReg or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CellReg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CellReg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CellReg

% Last Modified by GUIDE v2.5 19-Mar-2018 16:44:19

% reset figure properties to default:
if verLessThan('matlab','8.4')
    % MATLAB R2014a and earlier
    % if there are problems with GUI and figures change properties back to default
else
    % MATLAB R2014b and later 
    reset(0); 
end

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @CellReg_OpeningFcn, ...
    'gui_OutputFcn',  @CellReg_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT


% --- Executes just before CellReg is made visible.
function CellReg_OpeningFcn(hObject,~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)
% varargin   command line arguments to CellReg (see VARARGIN)

% Choose default command line output for CellReg
handles.output = hObject;

% Update handles struct
guidata(hObject, handles);

% UIWAIT makes CellReg wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% defining data struct:
data_struct=struct;
data_struct.sessions_list=[];

% Reseting figures and GUI parameters:
cla(handles.axes1,'reset')
axes(handles.axes1);
plot([],[])
set(gca,'xtick',[])
set(gca,'ytick',[])
cla(handles.axes2,'reset')
axes(handles.axes2);
plot([],[])
set(gca,'xtick',[])
set(gca,'ytick',[])
cla(handles.axes3,'reset')
axes(handles.axes3);
plot([],[])
set(gca,'xtick',[])
set(gca,'ytick',[])
cla(handles.axes4,'reset')
axes(handles.axes4);
plot([],[])
set(gca,'xtick',[])
set(gca,'ytick',[])
cla(handles.axes5,'reset')
axes(handles.axes5);
plot([],[])
set(gca,'xtick',[])
set(gca,'ytick',[])
cla(handles.axes6,'reset')
axes(handles.axes6);
plot([],[])
set(gca,'xtick',[])
set(gca,'ytick',[])
cla(handles.axes7,'reset')
axes(handles.axes7);
logo=imread('CellReg_Logo.png');
imagesc(logo);
set(gca,'xtick',[])
set(gca,'ytick',[])

data_struct.sessions_list=[];
set(handles.list_of_sessions,'value',1)
set(handles.list_of_sessions,'string',[]);
set(handles.green_session,'string','2')
set(handles.blue_session,'string','3')
set(handles.decision_thresh,'string','0.5')
set(handles.initial_p_same_slider,'value',0.5);
set(handles.initial_p_same_threshold,'string','0.5');
set(handles.final_p_same_slider,'value',0.5);
set(handles.model_maximal_distance,'string','14')
set(handles.distance_threshold,'string','5')
set(handles.correlation_threshold,'string','0.65')
set(handles.simple_distance_threshold,'string','5')
set(handles.simple_correlation_threshold,'string','0.65')
set(handles.figures_visibility_on,'Value',1);
set(handles.translations_rotations,'Value',1);
set(handles.spatial_correlations_2,'Value',1);
set(handles.spatial_correlations,'Value',1);
set(handles.use_model,'Value',1);
set(handles.microns_per_pixel,'string',[])
set(handles.microns_per_pixel,'value',0)
set(handles.microns_per_pixel,'backgroundColor',[1 1 1]);
set(handles.reference_session_index,'string','1')
set(handles.maximal_rotation','string','30')
set(handles.maximal_rotation','enable','on')
set(handles.transformation_smoothness,'string','2')
set(handles.transformation_smoothness,'enable','off')
set(handles.distance_threshold,'enable','off')
set(handles.correlation_threshold,'enable','on')
set(handles.decision_thresh,'enable','on')
set(handles.simple_distance_threshold,'enable','off')
set(handles.simple_correlation_threshold,'enable','off')
set(handles.comments,'string',[])
handles.data_struct=data_struct;
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = CellReg_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% Get default command line output from handles struct
varargout{1} = handles.output;


% --------------------------------------------------------------------
function load_new_data_Callback(hObject,~, handles)
% hObject    handle to load_new_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% Stage 1: Loading the spatial footprints of cellular activity from the different
% sessions.

% This callback loads a new data set which includes several sessions with
% their spatial footprints, centroid locations (optional), and events (optional). A
% single folder should be selected with all the mat files with number:
% example: "finalFiltersMat_1", "finalFiltersMat_2", "finalEventsMat_1", "finalEventsMat_2"

data_struct=handles.data_struct;

% choosing the files to load:
msgbox('Please choose the files containing the spatial footprints from all the sessions: ')
pause(3)
[file_names,files_path]=uigetfile('*.mat','MultiSelect','on');

number_of_sessions=size(file_names,2);
sessions_list=cell(1,number_of_sessions);
temp_file_names=cell(1,number_of_sessions);
for n=1:number_of_sessions
    sessions_list{1,n}=['Session ' num2str(n) ' - ' files_path file_names{1,n}];
    temp_file_names{1,n}=[files_path file_names{1,n}];
end
file_names=temp_file_names;

% defining the microns per pixel ratio:
microns_per_pixel=str2num(get(handles.microns_per_pixel,'string'));
if isempty(microns_per_pixel)
    msgbox('Please insert the pixel size in microns and press enter ')
    set(handles.microns_per_pixel,'backgroundColor',[1 0.5 0.5]);
    waitfor(handles.microns_per_pixel,'value',1);
    microns_per_pixel=str2num(get(handles.microns_per_pixel,'string'));
end
set(handles.microns_per_pixel,'string',num2str(round(100*microns_per_pixel)/100));
set(handles.microns_per_pixel,'backgroundColor',[1 1 1]);

% defining the results directory:
msgbox('Please select the folder in which the results will be saved')
pause(3)
results_directory=uigetdir(files_path); % the directory which the final results will be saved
figures_directory=fullfile(results_directory,'Figures');
if exist(figures_directory,'dir')~=7
    mkdir(figures_directory);
end

if get(handles.figures_visibility_on,'Value');
    figures_visibility='On';
else
    figures_visibility='Off';
end

% loading the spatial footprints:
disp('Stage 1 - Loading sessions')
[spatial_footprints,number_of_sessions]=load_multiple_sessions(file_names);
[footprints_projections]=compute_footprints_projections(spatial_footprints);
plot_all_sessions_projections(footprints_projections,figures_directory,figures_visibility)

% saving the loaded data into the data struct for the GUI
data_struct.results_directory=results_directory;
data_struct.figures_directory=figures_directory;
data_struct.microns_per_pixel=microns_per_pixel;
data_struct.spatial_footprints=spatial_footprints;
data_struct.footprints_projections=footprints_projections;
data_struct.number_of_sessions=number_of_sessions;
data_struct.sessions_list=sessions_list;
data_struct.file_names=file_names;
set(handles.list_of_sessions,'string',data_struct.sessions_list);
handles.data_struct=data_struct;
guidata(hObject, handles)
disp('Done')
msgbox('Finished loading sessions')


% --- Executes on button press in add_session.
function add_session_Callback(hObject,~, handles)
% hObject    handle to add_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% This callback adds another session to the list of sessions to be
% registered. The folder containg the filters, centroid_locations, and events
% (optional) should be selected

data_struct=handles.data_struct;

if get(handles.figures_visibility_on,'Value');
    figures_visibility='On';
else
    figures_visibility='Off';
end

if isfield(data_struct,'spatial_footprints') % some sessions were already loaded
    spatial_footprints=data_struct.spatial_footprints;
    number_of_sessions=data_struct.number_of_sessions;
    file_names=data_struct.file_names;
    sessions_list=data_struct.sessions_list;
    results_directory=data_struct.results_directory;
    
    % loading the session:
    msgbox('Please choose the file with the spatial footprints for this session: ')
    pause(1)
    [file_name,file_path]=uigetfile(results_directory,'*.mat','MultiSelect','off');
    number_of_sessions=number_of_sessions+1;
    file_names{number_of_sessions}=[file_path file_name];
    sessions_list{number_of_sessions}=['Session ' num2str(number_of_sessions) ' - ' file_path file_name];
    disp('Stage 1 - Loading sessions')
    [added_spatial_footprints]=load_single_session(file_names{number_of_sessions});
    spatial_footprints{number_of_sessions}=added_spatial_footprints;
    [added_footprints_projection]=compute_footprints_projections({added_spatial_footprints});
    footprints_projections{number_of_sessions}=added_footprints_projection;
    plot_single_session_projections(added_footprints_projection,num2str(number_of_sessions),figures_visibility)
else % first loaded session
    data_struct=handles.data_struct;
    
    % choosing the file to load:
    msgbox('Please choose the file with the spatial footprints for this session: ')
    pause(3)
    [file_name,file_path]=uigetfile('*.mat','MultiSelect','off');
    number_of_sessions=1;
    sessions_list={['Session 1 - ' file_path file_name]};
    file_names={[file_path file_name]};
    
    % defining the microns per pixel ratio:
    microns_per_pixel=str2num(get(handles.microns_per_pixel,'string'));
    if isempty(microns_per_pixel)
        msgbox('Please insert the pixel size in microns and press enter ')
        set(handles.microns_per_pixel,'backgroundColor',[1 0.5 0.5]);
        waitfor(handles.microns_per_pixel,'value',1);
        microns_per_pixel=str2num(get(handles.microns_per_pixel,'string'));
    end
    set(handles.microns_per_pixel,'string',num2str(round(100*microns_per_pixel)/100));
    set(handles.microns_per_pixel,'backgroundColor',[1 1 1]);
    
    % defining the results directory:
    msgbox('Please select the folder in which the results will be saved')
    pause(3)
    results_directory=uigetdir(file_path); % the directory which the final results will be saved
    figures_directory=fullfile(results_directory,'Figures');
    if exist(figures_directory,'dir')~=7
        mkdir(figures_directory);
    end
        
    % loading the spatial footprints:
    disp('Stage 1 - Loading sessions')
    [spatial_footprints]={load_single_session(file_names{1})};
    [footprints_projections]=compute_footprints_projections(spatial_footprints);
    plot_single_session_projections(footprints_projections,1,figures_visibility)
    
    % saving the loaded data into the data struct for the GUI
    data_struct.results_directory=results_directory;
    data_struct.figures_directory=figures_directory;
    data_struct.microns_per_pixel=microns_per_pixel;
end

% saving the loaded data into the data struct for the GUI
data_struct.spatial_footprints=spatial_footprints;
data_struct.footprints_projections=footprints_projections;
data_struct.number_of_sessions=number_of_sessions;
data_struct.sessions_list=sessions_list;
data_struct.file_names=file_names;
set(handles.list_of_sessions,'string',data_struct.sessions_list);
handles.data_struct=data_struct;
guidata(hObject, handles)
disp('Done')
msgbox('Finished loading session')

% --- Executes on button press in remove_session.
function remove_session_Callback(hObject,~, handles)
% hObject    handle to remove_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% This callback removes the selected session from the list of sessions to be registered

data_struct=handles.data_struct;
if isfield(data_struct,'spatial_footprints')
    % selecting the session to remove:
    number_of_sessions=data_struct.number_of_sessions;
    chosen_session=get(handles.list_of_sessions,'value');
    sessions_to_keep=setdiff(1:number_of_sessions,chosen_session);
    set(handles.list_of_sessions,'value',1)
    number_of_sessions=number_of_sessions-1;
    
    % removing the session from the data:
    if isfield(data_struct,'spatial_footprints_corrected') % if data was aligned
        centroid_locations=data_struct.centroid_locations(sessions_to_keep);
        adjusted_footprints_projections=data_struct.adjusted_footprints_projections(sessions_to_keep);
        spatial_footprints_corrected=data_struct.spatial_footprints_corrected(sessions_to_keep);
        centroid_locations_corrected=data_struct.centroid_locations_corrected(sessions_to_keep);
        footprints_projections_corrected=data_struct.footprints_projections_corrected(sessions_to_keep);
        
        % for variables that are compared to a reference session:
        reference_session_index=data_struct.reference_session_index;
        if reference_session_index==chosen_session
            errordlg('This session was used as a reference for alignment and therefore cannot be removed')
            error('This session was used as a reference for alignment and therefore cannot be removed')
        else
            sessions_without_reference=setdiff(1:number_of_sessions+1,reference_session_index);
            chosen_session_compared_to_reference_index=find(sessions_without_reference==chosen_session);
            sessions_to_keep_compared_to_reference=setdiff(1:number_of_sessions,chosen_session_compared_to_reference_index);
            maximal_cross_correlation=data_struct.maximal_cross_correlation(sessions_to_keep_compared_to_reference);
            alignment_translations=data_struct.alignment_translations(:,sessions_to_keep_compared_to_reference);
        end
        if reference_session_index>chosen_session % index of reference should change
            reference_session_index=reference_session_index-1;
        end
        
        data_struct.centroid_locations=centroid_locations;
        data_struct.adjusted_footprints_projections=adjusted_footprints_projections;
        data_struct.spatial_footprints_corrected=spatial_footprints_corrected;
        data_struct.centroid_locations_corrected=centroid_locations_corrected;
        data_struct.footprints_projections_corrected=footprints_projections_corrected;
        data_struct.maximal_cross_correlation=maximal_cross_correlation;
        data_struct.alignment_translations=alignment_translations;
        data_struct.reference_session_index=reference_session_index;
    end
    spatial_footprints=data_struct.spatial_footprints;
    footprints_projections=data_struct.footprints_projections;
    spatial_footprints=spatial_footprints(sessions_to_keep);
    footprints_projections=footprints_projections(sessions_to_keep);
    file_names=data_struct.file_names(sessions_to_keep);
    sessions_list=cell(1,number_of_sessions);
    for n=1:number_of_sessions
        sessions_list{1,n}=['Session ' num2str(n) ' - ' file_names{1,n}];
    end
    
    data_struct.number_of_sessions=number_of_sessions;
    data_struct.spatial_footprints=spatial_footprints;
    data_struct.footprints_projections=footprints_projections;
    data_struct.sessions_list=sessions_list;
    data_struct.file_names=file_names;
    data_struct.sessions_list=sessions_list;
    set(handles.list_of_sessions,'string',data_struct.sessions_list);
    handles.data_struct=data_struct;
    guidata(hObject, handles)
    msgbox('Finished removing session')
end

% --------------------------------------------------------------------
function load_transformed_data_Callback(hObject,~, handles)
% hObject    handle to load_transformed_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% This callback loads sessions that were already aligned into a reference coordinate system.
% For such data the compute model should be the next step.

msgbox('Please choose the file containing the aligned data structure: ')
pause(3)
[file_name,file_path]=uigetfile('*.mat','MultiSelect','off');
disp('Loading aligned data')
aligned_data_struct=load(fullfile(file_path,file_name));
if ~isstruct(aligned_data_struct)
    errordlg('This file does not contain data with the required format')
    error('This file does not contain data with the required format')
elseif ~isfield(aligned_data_struct,'aligned_data_struct')
    errordlg('This file does not contain data with the required format')
    error('This file does not contain data with the required format')
else
    % loading the aligned data:
    msgbox('Please select the folder in which the results will be saved')
    pause(3)
    results_directory=uigetdir(file_path); % the directory which the final results will be saved
    data_struct=aligned_data_struct.aligned_data_struct;
    data_struct.results_directory=results_directory;
    figures_directory=fullfile(results_directory,'Figures');
    data_struct.figures_directory=figures_directory;
    if exist(figures_directory,'dir')~=7
        mkdir(figures_directory);
    end
   
    % plotting the aligned data:
    footprints_projections_corrected=data_struct.footprints_projections_corrected;
    overlapping_FOV=data_struct.overlapping_FOV;
    if get(handles.figures_visibility_on,'Value');
        figures_visibility='On';
    else
        figures_visibility='Off';
    end
    plot_all_sessions_projections(footprints_projections_corrected,figures_directory,figures_visibility)
    number_of_sessions=length(footprints_projections_corrected);
    if number_of_sessions>2
        RGB_indexes=[1 2 3];
    else
        RGB_indexes=[1 2];
    end
    axes(handles.axes1);
    plot_RGB_overlay(footprints_projections_corrected,RGB_indexes,overlapping_FOV)
    figure('units','normalized','outerposition',[0.325 0.25 0.35 0.5],'Visible',figures_visibility)
    plot_RGB_overlay(footprints_projections_corrected,RGB_indexes,overlapping_FOV)

    % loading configurations to GUI:
    set(handles.microns_per_pixel,'string',num2str(round(100*data_struct.microns_per_pixel)/100));
    set(handles.reference_session_index,'string',num2str(data_struct.reference_session_index))
    set(handles.list_of_sessions,'string',data_struct.sessions_list)    
    if strcmp(data_struct.alignment_type,'Translations')
        set(handles.translations,'Value',1);
    elseif strcmp(data_struct.alignment_type,'Translations and Rotations')
        set(handles.translations_rotations,'Value',1);
    else
        set(handles.non_rigid,'Value',1);
    end
end

handles.data_struct=data_struct;
guidata(hObject, handles)
disp('Done')
msgbox('Finished loading aligned sessions')


% --- Executes on button press in load_modeled_data.
function load_modeled_data_Callback(hObject, eventdata, handles)
% hObject    handle to load_modeled_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This callback loads data that were already modeled.
% For such data the initial registration should be the next step.

msgbox('Please choose the file containing the modeled data structure: ')
pause(3)
[file_name,file_path]=uigetfile('*.mat','MultiSelect','off');
disp('Loading modeled data')
modeled_data_struct=load(fullfile(file_path,file_name));
if ~isstruct(modeled_data_struct)
    errordlg('This file does not contain data with the required format')
    error('This file does not contain data with the required format')
elseif ~isfield(modeled_data_struct,'modeled_data_struct')
    errordlg('This file does not contain data with the required format')
    error('This file does not contain data with the required format')
else
    % loading the aligned data:
    msgbox('Please select the folder in which the results will be saved')
    pause(3)
    results_directory=uigetdir(file_path); % the directory which the final results will be saved
    data_struct=modeled_data_struct.modeled_data_struct;
    data_struct.results_directory=results_directory;
    figures_directory=fullfile(results_directory,'Figures');
    data_struct.figures_directory=figures_directory;
    if exist(figures_directory,'dir')~=7
        mkdir(figures_directory);
    end
   
    % plotting the data:
    footprints_projections_corrected=data_struct.footprints_projections_corrected;
    overlapping_FOV=data_struct.overlapping_FOV;
    if get(handles.figures_visibility_on,'Value');
        figures_visibility='On';
    else
        figures_visibility='Off';
    end
    plot_all_sessions_projections(footprints_projections_corrected,figures_directory,figures_visibility)
    number_of_sessions=length(footprints_projections_corrected);
    if number_of_sessions>2
        RGB_indexes=[1 2 3];
    else
        RGB_indexes=[1 2];
    end
    axes(handles.axes1);
    plot_RGB_overlay(footprints_projections_corrected,RGB_indexes,overlapping_FOV)
    figure('units','normalized','outerposition',[0.325 0.25 0.35 0.5],'Visible',figures_visibility)
    plot_RGB_overlay(footprints_projections_corrected,RGB_indexes,overlapping_FOV)

    % loading configurations to GUI:
    set(handles.microns_per_pixel,'string',num2str(round(100*data_struct.microns_per_pixel)/100));
    set(handles.reference_session_index,'string',num2str(data_struct.reference_session_index))
    set(handles.list_of_sessions,'string',data_struct.sessions_list)
    set(handles.model_maximal_distance,'string',num2str(data_struct.maximal_distance));   
    if strcmp(data_struct.alignment_type,'Translations')
        set(handles.translations,'Value',1);
    elseif strcmp(data_struct.alignment_type,'Translations and Rotations')
        set(handles.translations_rotations,'Value',1);
    else
        set(handles.non_rigid,'Value',1);
    end
end

handles.data_struct=data_struct;
guidata(hObject, handles)
disp('Done')
msgbox('Finished loading modeled data')


% --- Executes on button press in transform_sessions.
function transform_sessions_Callback(hObject,~, handles)
% hObject    handle to transform_sessions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% Stage 2: Aligning all the sessions to a reference coordinate system using
% rigid-body transformation

% This callback performs rigid-body transfomration to all the sessions
% according to a chosen reference ssseion. This stage:
% 1) Corrects sessions for translation/rotations and transforming the
% spatial footprints into a single coordinate frame
% 2) Matches the sizes of all the spatial footprints from the different sessions
% to the intersction of the different FOVs
% 3) Evaluate whether or not it is suitable for
% longitudinal analysis

use_parallel_processing=true; % either true or false

data_struct=handles.data_struct;
number_of_sessions=data_struct.number_of_sessions;
spatial_footprints=data_struct.spatial_footprints;
results_directory=data_struct.results_directory;
figures_directory=data_struct.figures_directory;
microns_per_pixel=data_struct.microns_per_pixel;

% defining the aligned data structure:
aligned_data_struct=struct;
aligned_data_struct.number_of_sessions=number_of_sessions;
aligned_data_struct.spatial_footprints=spatial_footprints;
aligned_data_struct.results_directory=results_directory;
aligned_data_struct.figures_directory=figures_directory;
aligned_data_struct.microns_per_pixel=microns_per_pixel;
aligned_data_struct.footprints_projections=data_struct.footprints_projections;
aligned_data_struct.sessions_list=data_struct.sessions_list;
aligned_data_struct.file_names=data_struct.file_names;

if get(handles.figures_visibility_on,'Value');
    figures_visibility='On';
else
    figures_visibility='Off';
end

% Defining the parameters for image alignment:
translations_value=get(handles.translations,'Value');
rotations_value=get(handles.translations_rotations,'Value');
if translations_value==1
    alignment_type='Translations';
elseif rotations_value==1
    alignment_type='Translations and Rotations';
else
    alignment_type='Non-rigid';
end

if strcmp(alignment_type,'Translations and Rotations')
    maximal_rotation=str2num(get(handles.maximal_rotation,'string'));
end
if strcmp(alignment_type,'Non-rigid')
    transformation_smoothness=str2num(get(handles.transformation_smoothness,'string'));
    if transformation_smoothness>3 || transformation_smoothness<0.5
        errordlg('FOV smoothing parameter should be between 0.5-3')
        error('FOV smoothing parameter should be between 0.5-3')
    end
end

reference_session_index=str2num(get(handles.reference_session_index,'string'));
reference_valid=1;
if isempty(reference_session_index) || reference_session_index<1 || reference_session_index>number_of_sessions
    reference_valid=0;
end
while reference_valid==0
    set(handles.reference_session_index,'value',0)
    msgbox('Please insert a valid reference session number and press enter ');
    waitfor(handles.reference_session_index,'value',1);
    reference_session_index=str2num(get(handles.reference_session_index,'string'));
    if ~isempty(reference_session_index) && reference_session_index>=1 && reference_session_index<=number_of_sessions
        reference_valid=1;
    end
end

% Preparing the data for alignment:
disp('Stage 2 - Aligning sessions')
[normalized_spatial_footprints]=normalize_spatial_footprints(spatial_footprints);
[adjusted_spatial_footprints,adjusted_FOV,adjusted_x_size,adjusted_y_size,adjustment_zero_padding]=...
    adjust_FOV_size(normalized_spatial_footprints);
[adjusted_footprints_projections]=compute_footprints_projections(adjusted_spatial_footprints);
[centroid_locations]=compute_centroid_locations(adjusted_spatial_footprints,microns_per_pixel);
[centroid_projections]=compute_centroids_projections(centroid_locations,adjusted_spatial_footprints);

% Aligning the cells according to the tranlations/rotations that maximize their similarity:
sufficient_correlation_centroids=0.2; % smaller correlation imply no similarity between sessions
sufficient_correlation_footprints=0.2; % smaller correlation imply no similarity between sessions
if strcmp(alignment_type,'Translations and Rotations')
    [spatial_footprints_corrected,centroid_locations_corrected,footprints_projections_corrected,centroid_projections_corrected,maximal_cross_correlation,alignment_translations,overlapping_FOV]=...
        align_images(adjusted_spatial_footprints,centroid_locations,adjusted_footprints_projections,centroid_projections,adjusted_FOV,microns_per_pixel,reference_session_index,alignment_type,sufficient_correlation_centroids,sufficient_correlation_footprints,use_parallel_processing,maximal_rotation);
elseif strcmp(alignment_type,'Non-rigid')
    [spatial_footprints_corrected,centroid_locations_corrected,footprints_projections_corrected,centroid_projections_corrected,maximal_cross_correlation,alignment_translations,overlapping_FOV,displacement_fields]=...
        align_images(adjusted_spatial_footprints,centroid_locations,adjusted_footprints_projections,centroid_projections,adjusted_FOV,microns_per_pixel,reference_session_index,alignment_type,sufficient_correlation_centroids,sufficient_correlation_footprints,use_parallel_processing,transformation_smoothness);
else
    [spatial_footprints_corrected,centroid_locations_corrected,footprints_projections_corrected,centroid_projections_corrected,maximal_cross_correlation,alignment_translations,overlapping_FOV]=...
        align_images(adjusted_spatial_footprints,centroid_locations,adjusted_footprints_projections,centroid_projections,adjusted_FOV,microns_per_pixel,reference_session_index,alignment_type,sufficient_correlation_centroids,sufficient_correlation_footprints,use_parallel_processing);
end

% Evaluating data quality:
[all_projections_correlations,number_of_cells_per_session]=...
    evaluate_data_quality(spatial_footprints_corrected,centroid_projections_corrected,footprints_projections_corrected,maximal_cross_correlation,alignment_translations,reference_session_index,sufficient_correlation_footprints,alignment_type);

% plotting alignment results:
if strcmp(alignment_type,'Non-rigid')
    plot_alignment_results(adjusted_spatial_footprints,centroid_locations,spatial_footprints_corrected,centroid_locations_corrected,adjusted_footprints_projections,footprints_projections_corrected,reference_session_index,all_projections_correlations,maximal_cross_correlation,alignment_translations,overlapping_FOV,alignment_type,number_of_cells_per_session,figures_directory,figures_visibility,displacement_fields)
else
    plot_alignment_results(adjusted_spatial_footprints,centroid_locations,spatial_footprints_corrected,centroid_locations_corrected,adjusted_footprints_projections,footprints_projections_corrected,reference_session_index,all_projections_correlations,maximal_cross_correlation,alignment_translations,overlapping_FOV,alignment_type,number_of_cells_per_session,figures_directory,figures_visibility)
end
if number_of_sessions>2
    RGB_indexes=[1 2 3];
else
    RGB_indexes=[1 2];
end
axes(handles.axes1);
plot_RGB_overlay(footprints_projections_corrected,RGB_indexes,overlapping_FOV)

% saving the results into the data struct for the GUI
data_struct.reference_session_index=reference_session_index;
data_struct.alignment_type=alignment_type;
data_struct.centroid_locations=centroid_locations;
data_struct.spatial_footprints_corrected=spatial_footprints_corrected;
data_struct.centroid_locations_corrected=centroid_locations_corrected;
data_struct.adjusted_footprints_projections=adjusted_footprints_projections;
data_struct.footprints_projections_corrected=footprints_projections_corrected;
data_struct.adjusted_x_size=adjusted_x_size;
data_struct.adjusted_y_size=adjusted_y_size;
data_struct.overlapping_FOV=overlapping_FOV;
data_struct.maximal_cross_correlation=maximal_cross_correlation;
data_struct.alignment_translations=alignment_translations;
data_struct.adjustment_zero_padding=adjustment_zero_padding;

% saving the results into the aligned data structure:
aligned_data_struct.reference_session_index=reference_session_index;
aligned_data_struct.alignment_type=alignment_type;
aligned_data_struct.centroid_locations=centroid_locations;
aligned_data_struct.spatial_footprints_corrected=spatial_footprints_corrected;
aligned_data_struct.centroid_locations_corrected=centroid_locations_corrected;
aligned_data_struct.adjusted_footprints_projections=adjusted_footprints_projections;
aligned_data_struct.footprints_projections_corrected=footprints_projections_corrected;
aligned_data_struct.adjusted_x_size=adjusted_x_size;
aligned_data_struct.adjusted_y_size=adjusted_y_size;
aligned_data_struct.overlapping_FOV=overlapping_FOV;
aligned_data_struct.maximal_cross_correlation=maximal_cross_correlation;
aligned_data_struct.alignment_translations=alignment_translations;
aligned_data_struct.adjustment_zero_padding=adjustment_zero_padding;

handles.data_struct=data_struct;
disp('Saving the aligned data structure')
save(fullfile(results_directory,'aligned_data_struct.mat'),'aligned_data_struct','-v7.3')
guidata(hObject,handles)
if use_parallel_processing
    delete(gcp);
end
disp('Done')
msgbox('Finished aligning sessions')


% --- Executes on button press in compute_model.
function compute_model_Callback(hObject,~, handles)
% hObject    handle to compute_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% Stage 3: Computing a probabilistic model of the spatial footprints similarities
% of neighboring cell-pairs from different sessions using the centroid_locations
% distances and spatial correlations

% This callback computes the probability model for the same cells and
% different cells according to either spatial correlations, centroid
% distances, or both measures. The output is all the probabilities of
% neighboring cell-pairs to be the same cell - P_same.

% reseting figures
cla(handles.axes2,'reset')
axes(handles.axes2);
plot([],[])
set(gca,'xtick',[])
set(gca,'ytick',[])
cla(handles.axes3,'reset')
axes(handles.axes3);
plot([],[])
set(gca,'xtick',[])
set(gca,'ytick',[])
cla(handles.axes4,'reset')
axes(handles.axes4);
plot([],[])
set(gca,'xtick',[])
set(gca,'ytick',[])
cla(handles.axes5,'reset')
axes(handles.axes5);
plot([],[])
set(gca,'xtick',[])
set(gca,'ytick',[])
cla(handles.axes6,'reset')
axes(handles.axes6);
plot([],[])
set(gca,'xtick',[])
set(gca,'ytick',[])

data_struct=handles.data_struct;
spatial_footprints_corrected=data_struct.spatial_footprints_corrected;
centroid_locations_corrected=data_struct.centroid_locations_corrected;
microns_per_pixel=data_struct.microns_per_pixel;
results_directory=data_struct.results_directory;
figures_directory=data_struct.figures_directory;

% defining the modeled data structure:
modeled_data_struct=struct;
modeled_data_struct.number_of_sessions=data_struct.number_of_sessions;
modeled_data_struct.spatial_footprints=data_struct.spatial_footprints;
modeled_data_struct.results_directory=results_directory;
modeled_data_struct.figures_directory=figures_directory;
modeled_data_struct.microns_per_pixel=microns_per_pixel;
modeled_data_struct.reference_session_index=data_struct.reference_session_index;
modeled_data_struct.alignment_type=data_struct.alignment_type;
modeled_data_struct.centroid_locations=data_struct.centroid_locations;
modeled_data_struct.spatial_footprints_corrected=spatial_footprints_corrected;
modeled_data_struct.centroid_locations_corrected=centroid_locations_corrected;
modeled_data_struct.adjusted_footprints_projections=data_struct.adjusted_footprints_projections;
modeled_data_struct.footprints_projections_corrected=data_struct.footprints_projections_corrected;
modeled_data_struct.adjusted_x_size=data_struct.adjusted_x_size;
modeled_data_struct.adjusted_y_size=data_struct.adjusted_y_size;
modeled_data_struct.overlapping_FOV=data_struct.overlapping_FOV;
modeled_data_struct.maximal_cross_correlation=data_struct.maximal_cross_correlation;
modeled_data_struct.alignment_translations=data_struct.alignment_translations;
modeled_data_struct.adjustment_zero_padding=data_struct.adjustment_zero_padding;
modeled_data_struct.footprints_projections=data_struct.footprints_projections;
modeled_data_struct.sessions_list=data_struct.sessions_list;
modeled_data_struct.file_names=data_struct.file_names;

if get(handles.figures_visibility_on,'Value')
    figures_visibility='On';
else
    figures_visibility='Off';
end

% Defining the parameters for the probabilstic modeling:
maximal_distance=str2num(get(handles.model_maximal_distance,'string'));
normalized_maximal_distance=maximal_distance/microns_per_pixel;
p_same_certainty_threshold=0.95; % certain cells are those with p_same>threshld or <1-threshold
[number_of_bins,centers_of_bins]=estimate_number_of_bins(spatial_footprints_corrected,normalized_maximal_distance);

disp('Stage 3 - Calculating a probabilistic model of the data')
[all_to_all_indexes,all_to_all_spatial_correlations,all_to_all_centroid_distances,neighbors_spatial_correlations,neighbors_centroid_distances,neighbors_x_displacements,neighbors_y_displacements,NN_spatial_correlations,NNN_spatial_correlations,NN_centroid_distances,NNN_centroid_distances]=...
    compute_data_distribution(spatial_footprints_corrected,centroid_locations_corrected,normalized_maximal_distance);

% saving the results into the data struct for the GUI
data_struct.all_to_all_indexes=all_to_all_indexes;
data_struct.all_to_all_spatial_correlations=all_to_all_spatial_correlations;
data_struct.all_to_all_centroid_distances=all_to_all_centroid_distances;
data_struct.neighbors_spatial_correlations=neighbors_spatial_correlations;
data_struct.neighbors_centroid_distances=neighbors_centroid_distances;
data_struct.neighbors_x_displacements=neighbors_x_displacements;
data_struct.neighbors_y_displacements=neighbors_y_displacements;
data_struct.NN_spatial_correlations=NN_spatial_correlations;
data_struct.NNN_spatial_correlations=NNN_spatial_correlations;
data_struct.NN_centroid_distances=NN_centroid_distances;
data_struct.NNN_centroid_distances=NNN_centroid_distances;

% saving the results into the modeled data structure:
modeled_data_struct.all_to_all_indexes=all_to_all_indexes;
modeled_data_struct.all_to_all_spatial_correlations=all_to_all_spatial_correlations;
modeled_data_struct.all_to_all_centroid_distances=all_to_all_centroid_distances;
modeled_data_struct.neighbors_spatial_correlations=neighbors_spatial_correlations;
modeled_data_struct.neighbors_centroid_distances=neighbors_centroid_distances;
modeled_data_struct.neighbors_x_displacements=neighbors_x_displacements;
modeled_data_struct.neighbors_y_displacements=neighbors_y_displacements;
modeled_data_struct.NN_spatial_correlations=NN_spatial_correlations;
modeled_data_struct.NNN_spatial_correlations=NNN_spatial_correlations;
modeled_data_struct.NN_centroid_distances=NN_centroid_distances;
modeled_data_struct.NNN_centroid_distances=NNN_centroid_distances;

handles.data_struct=data_struct;
guidata(hObject, handles)

% Plotting the (x,y) displacements:
x_y_displacements=plot_x_y_displacements(neighbors_x_displacements,neighbors_y_displacements,microns_per_pixel,normalized_maximal_distance,number_of_bins,centers_of_bins,figures_directory,figures_visibility);
axes(handles.axes2)
plot_x_y_displacements_GUI(x_y_displacements,microns_per_pixel,centers_of_bins,normalized_maximal_distance,number_of_bins)

disp('Calculating a probabilistic model of the data')
% Modeling the distribution of centroid distances:
[centroid_distances_model_parameters,p_same_given_centroid_distance,centroid_distances_distribution,centroid_distances_model_same_cells,centroid_distances_model_different_cells,centroid_distances_model_weighted_sum,MSE_centroid_distances_model,centroid_distance_intersection]=...
    compute_centroid_distances_model(neighbors_centroid_distances,microns_per_pixel,centers_of_bins);

% Modeling the distribution of spatial correlations:
[spatial_correlations_model_parameters,p_same_given_spatial_correlation,spatial_correlations_distribution,spatial_correlations_model_same_cells,spatial_correlations_model_different_cells,spatial_correlations_model_weighted_sum,MSE_spatial_correlations_model,spatial_correlation_intersection]=...
    compute_spatial_correlations_model(neighbors_spatial_correlations,centers_of_bins);

% estimating registration accuracy:
[p_same_centers_of_bins,uncertain_fraction_centroid_distances,cdf_p_same_centroid_distances,false_positive_per_distance_threshold,true_positive_per_distance_threshold,uncertain_fraction_spatial_correlations,cdf_p_same_spatial_correlations,false_positive_per_correlation_threshold,true_positive_per_correlation_threshold]=...
    estimate_registration_accuracy(p_same_certainty_threshold,neighbors_centroid_distances,centroid_distances_model_same_cells,centroid_distances_model_different_cells,p_same_given_centroid_distance,centers_of_bins,neighbors_spatial_correlations,spatial_correlations_model_same_cells,spatial_correlations_model_different_cells,p_same_given_spatial_correlation);
% Checking which model is better according to a defined cost function:
[best_model_string]=choose_best_model(MSE_centroid_distances_model,centroid_distances_model_same_cells,centroid_distances_model_different_cells,p_same_given_centroid_distance,MSE_spatial_correlations_model,spatial_correlations_model_same_cells,spatial_correlations_model_different_cells,p_same_given_spatial_correlation);

% change the initial and final registration according to the best model:
if strcmp(best_model_string,'Spatial correlation')
    set(handles.spatial_correlations,'Value',1);
    set(handles.spatial_correlations_2,'Value',1);
    set(handles.distance_threshold,'enable','off')
    set(handles.correlation_threshold,'enable','on')
else
    set(handles.centroid_distances,'Value',1);
    set(handles.centroid_distances_2,'Value',1);
    set(handles.correlation_threshold,'enable','off')
    set(handles.distance_threshold,'enable','on')
end

% Plotting the probabilistic models and estimated registration accuracy:
plot_models(centroid_distances_model_parameters,NN_centroid_distances,NNN_centroid_distances,centroid_distances_distribution,centroid_distances_model_same_cells,centroid_distances_model_different_cells,centroid_distances_model_weighted_sum,centroid_distance_intersection,centers_of_bins,microns_per_pixel,normalized_maximal_distance,figures_directory,figures_visibility,spatial_correlations_model_parameters,NN_spatial_correlations,NNN_spatial_correlations,spatial_correlations_distribution,spatial_correlations_model_same_cells,spatial_correlations_model_different_cells,spatial_correlations_model_weighted_sum,spatial_correlation_intersection)
plot_estimated_registration_accuracy(p_same_centers_of_bins,p_same_certainty_threshold,p_same_given_centroid_distance,centroid_distances_distribution,cdf_p_same_centroid_distances,uncertain_fraction_centroid_distances,true_positive_per_distance_threshold,false_positive_per_distance_threshold,centers_of_bins,normalized_maximal_distance,microns_per_pixel,figures_directory,figures_visibility,p_same_given_spatial_correlation,spatial_correlations_distribution,cdf_p_same_spatial_correlations,uncertain_fraction_spatial_correlations,true_positive_per_correlation_threshold,false_positive_per_correlation_threshold)
plot_estimated_accuracy_GUI(handles,p_same_centers_of_bins,p_same_certainty_threshold,p_same_given_centroid_distance,centroid_distances_distribution,cdf_p_same_centroid_distances,true_positive_per_distance_threshold,false_positive_per_distance_threshold,centers_of_bins,normalized_maximal_distance,microns_per_pixel,p_same_given_spatial_correlation,spatial_correlations_distribution,cdf_p_same_spatial_correlations,true_positive_per_correlation_threshold,false_positive_per_correlation_threshold)

% Computing the P_same for each neighboring cell-pair according to the different models:
[all_to_all_p_same_centroid_distance_model,all_to_all_p_same_spatial_correlation_model]=...
    compute_p_same(all_to_all_centroid_distances,p_same_given_centroid_distance,centers_of_bins,all_to_all_spatial_correlations,p_same_given_spatial_correlation);

% saving the results into the data struct for GUI:
data_struct.best_model_string=best_model_string;
data_struct.maximal_distance=maximal_distance;
data_struct.number_of_bins=number_of_bins;
data_struct.centers_of_bins=centers_of_bins;

data_struct.false_positive_per_distance_threshold=false_positive_per_distance_threshold;
data_struct.true_positive_per_distance_threshold=true_positive_per_distance_threshold;
data_struct.cdf_p_same_centroid_distances=cdf_p_same_centroid_distances;
data_struct.uncertain_fraction_centroid_distances=uncertain_fraction_centroid_distances;
data_struct.p_same_given_centroid_distance=p_same_given_centroid_distance;
data_struct.neighbors_centroid_distances=neighbors_centroid_distances;
data_struct.MSE_centroid_distances_model=MSE_centroid_distances_model;
data_struct.centroid_distances_model_parameters=centroid_distances_model_parameters;
data_struct.centroid_distances_distribution=centroid_distances_distribution;
data_struct.centroid_distance_intersection=centroid_distance_intersection;
data_struct.all_to_all_p_same_centroid_distance_model=all_to_all_p_same_centroid_distance_model;

% saving the results into the modeled data structure:
modeled_data_struct.best_model_string=best_model_string;
modeled_data_struct.maximal_distance=maximal_distance;
modeled_data_struct.number_of_bins=number_of_bins;
modeled_data_struct.centers_of_bins=centers_of_bins;

modeled_data_struct.false_positive_per_distance_threshold=false_positive_per_distance_threshold;
modeled_data_struct.true_positive_per_distance_threshold=true_positive_per_distance_threshold;
modeled_data_struct.cdf_p_same_centroid_distances=cdf_p_same_centroid_distances;
modeled_data_struct.uncertain_fraction_centroid_distances=uncertain_fraction_centroid_distances;
modeled_data_struct.p_same_given_centroid_distance=p_same_given_centroid_distance;
modeled_data_struct.neighbors_centroid_distances=neighbors_centroid_distances;
modeled_data_struct.MSE_centroid_distances_model=MSE_centroid_distances_model;
modeled_data_struct.centroid_distances_model_parameters=centroid_distances_model_parameters;
modeled_data_struct.centroid_distances_distribution=centroid_distances_distribution;
modeled_data_struct.centroid_distance_intersection=centroid_distance_intersection;
modeled_data_struct.all_to_all_p_same_centroid_distance_model=all_to_all_p_same_centroid_distance_model;

data_struct.false_positive_per_correlation_threshold=false_positive_per_correlation_threshold;
data_struct.true_positive_per_correlation_threshold=true_positive_per_correlation_threshold;
data_struct.cdf_p_same_spatial_correlations=cdf_p_same_spatial_correlations;
data_struct.uncertain_fraction_spatial_correlations=uncertain_fraction_spatial_correlations;
data_struct.all_to_all_spatial_correlations=all_to_all_spatial_correlations;
data_struct.MSE_spatial_correlations_model=MSE_spatial_correlations_model;
data_struct.spatial_correlations_model_parameters=spatial_correlations_model_parameters;
data_struct.p_same_given_spatial_correlation=p_same_given_spatial_correlation;
data_struct.spatial_correlations_distribution=spatial_correlations_distribution;
data_struct.spatial_correlation_intersection=spatial_correlation_intersection;
data_struct.all_to_all_p_same_spatial_correlation_model=all_to_all_p_same_spatial_correlation_model;

% saving the results into the modeled data structure:
modeled_data_struct.false_positive_per_correlation_threshold=false_positive_per_correlation_threshold;
modeled_data_struct.true_positive_per_correlation_threshold=true_positive_per_correlation_threshold;
modeled_data_struct.cdf_p_same_spatial_correlations=cdf_p_same_spatial_correlations;
modeled_data_struct.uncertain_fraction_spatial_correlations=uncertain_fraction_spatial_correlations;
modeled_data_struct.all_to_all_spatial_correlations=all_to_all_spatial_correlations;
modeled_data_struct.MSE_spatial_correlations_model=MSE_spatial_correlations_model;
modeled_data_struct.spatial_correlations_model_parameters=spatial_correlations_model_parameters;
modeled_data_struct.p_same_given_spatial_correlation=p_same_given_spatial_correlation;
modeled_data_struct.spatial_correlations_distribution=spatial_correlations_distribution;
modeled_data_struct.spatial_correlation_intersection=spatial_correlation_intersection;
modeled_data_struct.all_to_all_p_same_spatial_correlation_model=all_to_all_p_same_spatial_correlation_model;

% setting the intersection point as the threshold
set(handles.correlation_threshold,'string',num2str(spatial_correlation_intersection))
set(handles.distance_threshold,'string',num2str(centroid_distance_intersection))
handles.data_struct=data_struct;
disp('Saving the modeled data structure')
save(fullfile(results_directory,'modeled_data_struct.mat'),'modeled_data_struct','-v7.3')
guidata(hObject, handles)
disp('Done')
msgbox(['Finished computing probabilistic model - The ' best_model_string ' model is best suited for the data'])


% --- Executes on button press in register_cells_initial.
function register_cells_initial_Callback(hObject,~,handles)
% hObject    handle to register_cells_initial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% Stage 4: Obtaining an initial cell registration according to an optimized
% registration threshold.

% This callback performs initial cell registration according to either
% spatial correlations or centroid distances:
data_struct=handles.data_struct;

microns_per_pixel=data_struct.microns_per_pixel;
spatial_footprints_corrected=data_struct.spatial_footprints_corrected;
centroid_locations_corrected=data_struct.centroid_locations_corrected;
normalized_maximal_distance=data_struct.maximal_distance/microns_per_pixel;
number_of_bins=data_struct.number_of_bins;
figures_directory=data_struct.figures_directory;

if get(handles.figures_visibility_on,'Value');
    figures_visibility='On';
else
    figures_visibility='Off';
end

% Computing the initial registration according to a simple threshold:
if get(handles.spatial_correlations,'Value')==1 % if spatial correlations are used    
        initial_registration_type='Spatial correlation';
        initial_threshold=str2num(get(handles.correlation_threshold,'string'));
        [cell_to_index_map,registered_cells_spatial_correlations,non_registered_cells_spatial_correlations]=...
            initial_registration_spatial_correlations(normalized_maximal_distance,initial_threshold,spatial_footprints_corrected,centroid_locations_corrected);
        plot_initial_registration(cell_to_index_map,number_of_bins,spatial_footprints_corrected,initial_registration_type,figures_directory,figures_visibility,registered_cells_spatial_correlations,non_registered_cells_spatial_correlations)    
else
    initial_registration_type='Centroid distances';
    initial_threshold=str2num(get(handles.distance_threshold,'string'));
    centroid_distances_distribution_threshold=initial_threshold/microns_per_pixel;
    [cell_to_index_map,registered_cells_centroid_distances,non_registered_cells_centroid_distances]=...
        initial_registration_centroid_distances(normalized_maximal_distance,centroid_distances_distribution_threshold,centroid_locations_corrected);
    plot_initial_registration(cell_to_index_map,number_of_bins,spatial_footprints_corrected,initial_registration_type,figures_directory,figures_visibility,registered_cells_centroid_distances,non_registered_cells_centroid_distances,microns_per_pixel,normalized_maximal_distance)
end

disp([num2str(size(cell_to_index_map,1)) ' cells were found'])
disp('Done')

data_struct.initial_registration_type=initial_registration_type;
data_struct.cell_to_index_map=cell_to_index_map;
data_struct.initial_threshold=initial_threshold;
data_struct.initial_registration_type=initial_registration_type;

handles.data_struct=data_struct;
guidata(hObject, handles)
msgbox(['Finished performing initial cell registration - ' num2str(size(cell_to_index_map,1)) ' were found'])


% --- Executes on button press in register_cells_final.
function register_cells_final_Callback(hObject,~, handles)
% hObject    handle to register_cells_final (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% Stage 5: Obtaining the final cell registration based on a correlation clustering algorithm.

% This callback performs the final cell registration according to the
% probability model for same cells and different cells. P_same can be
% either according to centroid distances, spatial correlations or both:

data_struct=handles.data_struct;

if ~isfield(data_struct,'cell_to_index_map')
    errordlg('Final registration cannot be performed before initial registration')
    error('Final registration cannot be performed before initial registration')
end

if  get(handles.spatial_correlations_2,'Value')==1  
    if isfield(data_struct,'all_to_all_p_same_spatial_correlation_model')
        all_to_all_p_same_spatial_correlation_model=data_struct.all_to_all_p_same_spatial_correlation_model;
        all_to_all_indexes=data_struct.all_to_all_indexes;
        all_to_all_spatial_correlations=data_struct.all_to_all_spatial_correlations;
    else
        errordlg('Please compute the spatial correlations probability model before performing final cell registration')
        error('Please compute the spatial correlations probability model before performing final cell registration')
    end
elseif get(handles.centroid_distances_2,'Value')==1
    if isfield(data_struct,'all_to_all_p_same_centroid_distance_model')
        all_to_all_p_same_centroid_distance_model=data_struct.all_to_all_p_same_centroid_distance_model;
        all_to_all_indexes=data_struct.all_to_all_indexes;
        all_to_all_centroid_distances=data_struct.all_to_all_centroid_distances;
    end
end

results_directory=data_struct.results_directory;
figures_directory=data_struct.figures_directory;
overlapping_FOV=data_struct.overlapping_FOV;
cell_to_index_map=data_struct.cell_to_index_map;
centroid_locations_corrected=data_struct.centroid_locations_corrected;
spatial_footprints_corrected=data_struct.spatial_footprints_corrected;
number_of_sessions=data_struct.number_of_sessions;
microns_per_pixel=data_struct.microns_per_pixel;
maximal_distance=data_struct.maximal_distance;
normalized_maximal_distance=maximal_distance/microns_per_pixel;

if get(handles.figures_visibility_on,'Value');
    figures_visibility='On';
else
    figures_visibility='Off';
end

if get(handles.spatial_correlations_2,'Value')==1;
    model_type='Spatial correlation';
else
    model_type='Centroid distance';
end

transform_data=false;
if get(handles.use_model,'Value')==1;
    p_same_threshold=str2num(get(handles.decision_thresh,'string'));
    final_threshold=p_same_threshold;
    registration_approach='Probabilistic';
else
    registration_approach='Simple threshold';
    if strcmp(model_type,'Spatial correlation')
        final_threshold=str2num(get(handles.simple_correlation_threshold,'string'));
    elseif strcmp(model_type,'Centroid distance')
        final_threshold=str2num(get(handles.simple_distance_threshold,'string'));
        centroid_distances_distribution_threshold=(maximal_distance-final_threshold)/maximal_distance;
        transform_data=true;
    end
end

data_struct.final_threshold=final_threshold;
data_struct.registration_approach=registration_approach;
data_struct.model_type=model_type;

% Registering the cells with the clustering algorithm:
disp('Stage 5 - Performing final registration')
if strcmp(registration_approach,'Probabilistic')
    if strcmp(model_type,'Spatial correlation')
        [optimal_cell_to_index_map,registered_cells_centroids,cell_scores,cell_scores_positive,cell_scores_negative,cell_scores_exclusive,p_same_registered_pairs]=...
            cluster_cells(cell_to_index_map,all_to_all_p_same_spatial_correlation_model,all_to_all_indexes,normalized_maximal_distance,p_same_threshold,centroid_locations_corrected,registration_approach,transform_data);
    elseif strcmp(model_type,'Centroid distance')
        [optimal_cell_to_index_map,registered_cells_centroids,cell_scores,cell_scores_positive,cell_scores_negative,cell_scores_exclusive,p_same_registered_pairs]=...
            cluster_cells(cell_to_index_map,all_to_all_p_same_centroid_distance_model,all_to_all_indexes,normalized_maximal_distance,p_same_threshold,centroid_locations_corrected,registration_approach,transform_data);
    end
    plot_cell_scores(cell_scores_positive,cell_scores_negative,cell_scores_exclusive,cell_scores,p_same_registered_pairs,figures_directory,figures_visibility)
elseif strcmp(registration_approach,'Simple threshold')
    if strcmp(model_type,'Spatial correlation')
        [optimal_cell_to_index_map,registered_cells_centroids]=...
            cluster_cells(cell_to_index_map,all_to_all_spatial_correlations,all_to_all_indexes,normalized_maximal_distance,final_threshold,centroid_locations_corrected,registration_approach,transform_data);
    elseif strcmp(model_type,'Centroid distance')
        [optimal_cell_to_index_map,registered_cells_centroids]=...
            cluster_cells(cell_to_index_map,all_to_all_centroid_distances,all_to_all_indexes,normalized_maximal_distance,centroid_distances_distribution_threshold,centroid_locations_corrected,registration_approach,transform_data);
    end
end
[is_in_overlapping_FOV]=check_if_in_overlapping_FOV(registered_cells_centroids,overlapping_FOV);


% Plotting the registration results with the cell maps from all sessions:
plot_all_registered_projections(spatial_footprints_corrected,optimal_cell_to_index_map,figures_directory,figures_visibility)

% saving the clustering results:
disp('Saving the results')
cell_registered_struct=struct;
cell_registered_struct.cell_to_index_map=optimal_cell_to_index_map;
if strcmp(registration_approach,'Probabilistic');
    cell_registered_struct.cell_scores=cell_scores';
    cell_registered_struct.true_positive_scores=cell_scores_positive';
    cell_registered_struct.true_negative_scores=cell_scores_negative';
    cell_registered_struct.exclusivity_scores=cell_scores_exclusive';
    cell_registered_struct.p_same_registered_pairs=p_same_registered_pairs';
end
cell_registered_struct.is_cell_in_overlapping_FOV=is_in_overlapping_FOV';
cell_registered_struct.registered_cells_centroids=registered_cells_centroids';
cell_registered_struct.centroid_locations_corrected=centroid_locations_corrected';
cell_registered_struct.spatial_footprints_corrected=spatial_footprints_corrected';
cell_registered_struct.alignment_x_translations=data_struct.alignment_translations(1,:);
cell_registered_struct.alignment_y_translations=data_struct.alignment_translations(2,:);
if strcmp(data_struct.alignment_type,'Translations and Rotations')
    cell_registered_struct.alignment_rotations=data_struct.alignment_translations(3,:);
end
cell_registered_struct.adjustment_x_zero_padding=data_struct.adjustment_zero_padding(1,:);
cell_registered_struct.adjustment_y_zero_padding=data_struct.adjustment_zero_padding(2,:);
save(fullfile(results_directory,['cellRegistered_' datestr(clock,'yyyymmdd_HHMMss') '.mat']),'cell_registered_struct','-v7.3')

% Saving a log file with all the chosen parameters:
comments=get(handles.comments,'string');
file_names=data_struct.file_names;
adjusted_x_size=data_struct.adjusted_x_size;
adjusted_y_size=data_struct.adjusted_y_size;
alignment_type=data_struct.alignment_type;
reference_session_index=data_struct.reference_session_index;
number_of_bins=data_struct.number_of_bins;
initial_registration_type=data_struct.initial_registration_type;
initial_threshold=data_struct.initial_threshold;
if strcmp(registration_approach,'Probabilistic')
    if strcmp(model_type,'Spatial correlation')
        uncertain_fraction_spatial_correlations=data_struct.uncertain_fraction_spatial_correlations;
        false_positive_per_correlation_threshold=data_struct.false_positive_per_correlation_threshold;
        true_positive_per_correlation_threshold=data_struct.true_positive_per_correlation_threshold;
        MSE_spatial_correlations_model=data_struct.MSE_spatial_correlations_model;
        save_log_file(results_directory,file_names,microns_per_pixel,adjusted_x_size,adjusted_y_size,alignment_type,reference_session_index,maximal_distance,number_of_bins,initial_registration_type,initial_threshold,registration_approach,model_type,final_threshold,optimal_cell_to_index_map,cell_registered_struct,comments,uncertain_fraction_spatial_correlations,false_positive_per_correlation_threshold,true_positive_per_correlation_threshold,MSE_spatial_correlations_model)
    elseif strcmp(model_type,'Centroid distance')
        uncertain_fraction_centroid_distances=data_struct.uncertain_fraction_centroid_distances;
        false_positive_per_distance_threshold=data_struct.false_positive_per_distance_threshold;
        true_positive_per_distance_threshold=data_struct.true_positive_per_distance_threshold;
        MSE_centroid_distances_model=data_struct.MSE_centroid_distances_model;
        save_log_file(results_directory,file_names,microns_per_pixel,adjusted_x_size,adjusted_y_size,alignment_type,reference_session_index,maximal_distance,number_of_bins,initial_registration_type,initial_threshold,registration_approach,model_type,final_threshold,optimal_cell_to_index_map,cell_registered_struct,comments,uncertain_fraction_centroid_distances,false_positive_per_distance_threshold,true_positive_per_distance_threshold,MSE_centroid_distances_model)
    end
elseif strcmp(registration_approach,'Simple threshold')
    if strcmp(model_type,'Spatial correlation')
        save_log_file(results_directory,file_names,microns_per_pixel,adjusted_x_size,adjusted_y_size,alignment_type,reference_session_index,maximal_distance,number_of_bins,initial_registration_type,initial_threshold,registration_approach,model_type,final_threshold,optimal_cell_to_index_map,cell_registered_struct,comments)
    elseif strcmp(model_type,'Centroid distance')
        save_log_file(results_directory,file_names,microns_per_pixel,adjusted_x_size,adjusted_y_size,alignment_type,reference_session_index,maximal_distance,number_of_bins,initial_registration_type,initial_threshold,registration_approach,model_type,final_threshold,optimal_cell_to_index_map,cell_registered_struct,comments)
    end
end
disp([num2str(size(optimal_cell_to_index_map,1)) ' cells were found'])
disp('End of cell registration procedure')

handles.data_struct=data_struct;
guidata(hObject, handles)
msgbox(['Finished performing final cell registration - ' num2str(size(optimal_cell_to_index_map,1)) ' were found'])


% --- Executes on button press in reset.
function reset_Callback(hObject,~, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% defining data struct:
guidata(hObject, handles);
data_struct=struct;
data_struct.sessions_list=[];

% Reseting Figures and GUI parameters:
cla(handles.axes1,'reset')
axes(handles.axes1);
plot([],[])
set(gca,'xtick',[])
set(gca,'ytick',[])
cla(handles.axes2,'reset')
axes(handles.axes2);
plot([],[])
set(gca,'xtick',[])
set(gca,'ytick',[])
cla(handles.axes3,'reset')
axes(handles.axes3);
plot([],[])
set(gca,'xtick',[])
set(gca,'ytick',[])
cla(handles.axes4,'reset')
axes(handles.axes4);
plot([],[])
set(gca,'xtick',[])
set(gca,'ytick',[])
cla(handles.axes5,'reset')
axes(handles.axes5);
plot([],[])
set(gca,'xtick',[])
set(gca,'ytick',[])
cla(handles.axes6,'reset')
axes(handles.axes6);
plot([],[])
set(gca,'xtick',[])
set(gca,'ytick',[])

set(handles.list_of_sessions,'value',1)
set(handles.list_of_sessions,'string',[]);
set(handles.red_session,'string','1')
set(handles.green_session,'string','2')
set(handles.blue_session,'string','3')
set(handles.decision_thresh,'string','0.5')
set(handles.initial_p_same_slider,'value',0.5);
set(handles.initial_p_same_threshold,'string','0.5');
set(handles.final_p_same_slider,'value',0.5);
set(handles.model_maximal_distance,'string','14')
set(handles.distance_threshold,'string','5')
set(handles.correlation_threshold,'string','0.65')
set(handles.simple_distance_threshold,'string','5')
set(handles.simple_correlation_threshold,'string','0.65')
set(handles.figures_visibility_on,'Value',1);
set(handles.translations_rotations,'Value',1);
set(handles.spatial_correlations_2,'Value',1);
set(handles.use_model,'Value',1);
set(handles.spatial_correlations,'Value',1);
set(handles.microns_per_pixel,'value',0)
set(handles.microns_per_pixel,'string',[])
set(handles.microns_per_pixel,'backgroundColor',[1 1 1]);
set(handles.reference_session_index,'string','1')
set(handles.maximal_rotation','string','30')
set(handles.maximal_rotation','enable','on')
set(handles.transformation_smoothness,'string','2')
set(handles.transformation_smoothness,'enable','off')
set(handles.distance_threshold,'enable','off')
set(handles.correlation_threshold,'enable','on')
set(handles.simple_distance_threshold,'enable','off')
set(handles.simple_correlation_threshold,'enable','off')
set(handles.comments,'string',[])

handles.data_struct=data_struct;
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function correlation_threshold_CreateFcn(~,hObject,~)
% hObject    handle to correlation_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function distance_threshold_CreateFcn(hObject,~,~)
% hObject    handle to distance_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function list_of_sessions_CreateFcn(hObject, ~,~)
% hObject    handle to list_of_sessions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function maximal_distance_CreateFcn(hObject,~,~)
% hObject    handle to maximal_distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function number_of_bins_CreateFcn(hObject,~,~)
% hObject    handle to number_of_bins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function compute_model_CreateFcn(~,~,~)
% hObject    handle to compute_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function transformation_type_CreateFcn(~,~,~)
% hObject    handle to transformation_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function microns_per_pixel_Callback(~,~, handles)
% hObject    handle to microns_per_pixel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of microns_per_pixel as text
%        str2num(get(hObject,'String')) returns contents of microns_per_pixel as a double
microns_per_pixel=get(handles.microns_per_pixel,'string');
if ~isempty(microns_per_pixel)
    set(handles.microns_per_pixel,'value',1)
end


% --- Executes during object creation, after setting all properties.
function microns_per_pixel_CreateFcn(hObject,~,~)
% hObject    handle to microns_per_pixel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function reference_session_index_Callback(~,~, handles)
% hObject    handle to reference_session_index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of reference_session_index as text
%        str2num(get(hObject,'String')) returns contents of reference_session_index as a double

reference_session_index=get(handles.reference_session_index,'string');
if ~isempty(reference_session_index)
    set(handles.reference_session_index,'value',1)
end

% --- Executes during object creation, after setting all properties.
function reference_session_index_CreateFcn(hObject, ~, ~)
% hObject    handle to reference_session_index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function model_maximal_distance_CreateFcn(hObject,~,~)
% hObject    handle to model_maximal_distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function cluster_maximal_distance_CreateFcn(hObject,~,~)
% hObject    handle to cluster_maximal_distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function max_iterations_CreateFcn(hObject,~,~)
% hObject    handle to max_iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function decision_thresh_CreateFcn(hObject,~,~)
% hObject    handle to decision_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function load_Callback(~,~,~)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)


% --- Executes on selection change in list_of_sessions.
function list_of_sessions_Callback(hObject, eventdata, handles)
% hObject    handle to list_of_sessions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list_of_sessions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_of_sessions


function model_maximal_distance_Callback(hObject, eventdata, handles)
% hObject    handle to model_maximal_distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of model_maximal_distance as text
%        str2num(get(hObject,'String')) returns contents of model_maximal_distance as a double


function number_of_bins_Callback(hObject, eventdata, handles)
% hObject    handle to number_of_bins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_of_bins as text
%        str2num(get(hObject,'String')) returns contents of number_of_bins as a double


% --- Executes when selected object is changed in transformation_type.
function transformation_type_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in transformation_type
% eventdata  struct with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    struct with handles and user data (see GUIDATA)

if get(handles.translations,'Value')==1;
    set(handles.maximal_rotation','enable','off')
else
    set(handles.maximal_rotation','enable','on')
    set(handles.maximal_rotation','string','30')
end


function maximal_rotation_Callback(hObject, eventdata, handles)
% hObject    handle to maximal_rotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maximal_rotation as text
%        str2num(get(hObject,'String')) returns contents of maximal_rotation as a double


% --- Executes during object creation, after setting all properties.
function maximal_rotation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maximal_rotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when selected object is changed in initial_register_select.
function initial_register_select_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in initial_register_select
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

if get(handles.spatial_correlations,'Value')==1;
    set(handles.distance_threshold,'enable','off')
    set(handles.correlation_threshold,'enable','on')
else
    set(handles.correlation_threshold,'enable','off')
    set(handles.distance_threshold,'enable','on')
end


function comments_Callback(hObject, eventdata, handles)
% hObject    handle to comments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of comments as text
%        str2num(get(hObject,'String')) returns contents of comments as a double


% --- Executes during object creation, after setting all properties.
function comments_CreateFcn(hObject, eventdata, handles)
% hObject    handle to comments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in display_rgb.
function display_rgb_Callback(hObject, eventdata, handles)
% hObject    handle to display_rgb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

data_struct=handles.data_struct;

if ~isfield(data_struct,'footprints_projections_corrected')
    errordlg('RGB overlay cannot be displayed before transformation is performed')
else
    spatial_footprints_projections=data_struct.footprints_projections_corrected;
    overlapping_FOV=data_struct.overlapping_FOV;
    
    if get(handles.figures_visibility_on,'Value');
        figures_visibility='On';
    else
        figures_visibility='Off';
    end
    number_of_sessions=data_struct.number_of_sessions;
    red_session=str2num(get(handles.red_session,'string'));
    green_session=str2num(get(handles.green_session,'string'));
    if number_of_sessions>2
        blue_session=str2num(get(handles.blue_session,'string'));
        RGB_indexes=[red_session green_session blue_session];
    else
        RGB_indexes=[red_session green_session];
    end
    axes(handles.axes1);
    plot_RGB_overlay(spatial_footprints_projections,RGB_indexes,overlapping_FOV)
    figure('units','normalized','outerposition',[0.325 0.25 0.35 0.5],'Visible',figures_visibility)
    plot_RGB_overlay(spatial_footprints_projections,RGB_indexes,overlapping_FOV)
end

handles.data_struct=data_struct;
guidata(hObject, handles)


function red_session_Callback(hObject, eventdata, handles)
% hObject    handle to red_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of red_session as text
%        str2num(get(hObject,'String')) returns contents of red_session as a double


% --- Executes during object creation, after setting all properties.
function red_session_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function green_session_Callback(hObject, eventdata, handles)
% hObject    handle to green_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of green_session as text
%        str2num(get(hObject,'String')) returns contents of green_session as a double


% --- Executes during object creation, after setting all properties.
function green_session_CreateFcn(hObject, eventdata, handles)
% hObject    handle to green_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function blue_session_Callback(hObject, eventdata, handles)
% hObject    handle to blue_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of blue_session as text
%        str2num(get(hObject,'String')) returns contents of blue_session as a double


% --- Executes during object creation, after setting all properties.
function blue_session_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blue_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function initial_p_same_slider_Callback(hObject, eventdata, handles)
% hObject    handle to initial_p_same_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

data_struct=handles.data_struct;
initial_p_same_slider_value=round(100*get(handles.initial_p_same_slider,'value'))/100;
set(handles.initial_p_same_threshold,'string',num2str(initial_p_same_slider_value));
centers_of_bins=data_struct.centers_of_bins;

if get(handles.spatial_correlations,'Value')==1
    spatial_correlations_distribution=data_struct.spatial_correlations_distribution;
    p_same_given_spatial_correlation=data_struct.p_same_given_spatial_correlation;
    [~,p_same_ind]=min(abs(initial_p_same_slider_value-(p_same_given_spatial_correlation)));
    spatial_correlation_threshold=round(100*centers_of_bins{2}(p_same_ind))/100;
    set(handles.correlation_threshold,'string',num2str(spatial_correlation_threshold));
    axes(handles.axes4)
    cla(handles.axes4,'reset')
    plot_p_same_spatial_correlation_slider(spatial_correlations_distribution,p_same_given_spatial_correlation,spatial_correlation_threshold,centers_of_bins)
else
    centroid_distances_distribution=data_struct.centroid_distances_distribution;
    p_same_given_centroid_distance=data_struct.p_same_given_centroid_distance;
    microns_per_pixel=data_struct.microns_per_pixel;
    maximal_distance=data_struct.maximal_distance/microns_per_pixel;
    [~,p_same_ind]=min(abs(initial_p_same_slider_value-(p_same_given_centroid_distance)));
    centroid_distance_threshold=round(microns_per_pixel*100*centers_of_bins{1}(p_same_ind))/100;
    set(handles.distance_threshold,'string',num2str(centroid_distance_threshold));
    axes(handles.axes3)
    cla(handles.axes3,'reset')
    plot_p_same_centroid_distance_slider(centroid_distances_distribution,p_same_given_centroid_distance,centroid_distance_threshold,centers_of_bins,maximal_distance,microns_per_pixel)
end

handles.data_struct=data_struct;
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function initial_p_same_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initial_p_same_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function initial_p_same_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to initial_all_to_all_centroid_distances (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of initial_all_to_all_centroid_distances as text
%        str2num(get(hObject,'String')) returns contents of initial_all_to_all_centroid_distances as a double

data_struct=handles.data_struct;
initial_p_same_slider_value=str2num(get(handles.initial_p_same_threshold,'string'));
set(handles.initial_p_same_slider,'value',initial_p_same_slider_value);
centers_of_bins=data_struct.centers_of_bins;

if get(handles.spatial_correlations,'Value')==1
    spatial_correlations_distribution=data_struct.spatial_correlations_distribution;
    p_same_given_spatial_correlation=data_struct.p_same_given_spatial_correlation;
    [~,p_same_ind]=min(abs(initial_p_same_slider_value-(p_same_given_spatial_correlation)));
    spatial_correlation_threshold=round(100*centers_of_bins{2}(p_same_ind))/100;
    set(handles.correlation_threshold,'string',num2str(spatial_correlation_threshold));
    axes(handles.axes4)
    cla(handles.axes4,'reset')
    plot_p_same_spatial_correlation_slider(spatial_correlations_distribution,p_same_given_spatial_correlation,spatial_correlation_threshold,centers_of_bins)
else
    centroid_distances_distribution=data_struct.centroid_distances_distribution;
    p_same_given_centroid_distance=data_struct.p_same_given_centroid_distance;
    microns_per_pixel=data_struct.microns_per_pixel;
    maximal_distance=data_struct.maximal_distance/microns_per_pixel;
    [~,p_same_ind]=min(abs(initial_p_same_slider_value-(p_same_given_centroid_distance)));
    centroid_distance_threshold=round(microns_per_pixel*100*centers_of_bins{1}(p_same_ind))/100;
    set(handles.distance_threshold,'string',num2str(centroid_distance_threshold));
    axes(handles.axes3)
    cla(handles.axes3,'reset')
    plot_p_same_centroid_distance_slider(centroid_distances_distribution,p_same_given_centroid_distance,centroid_distance_threshold,centers_of_bins,maximal_distance,microns_per_pixel)
end

handles.data_struct=data_struct;
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function initial_p_same_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initial_all_to_all_centroid_distances (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in default_initial_registration.
function default_initial_registration_Callback(hObject, eventdata, handles)
% hObject    handle to default_initial_registration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

initial_p_same_slider_value=0.5;
data_struct=handles.data_struct;
set(handles.initial_p_same_slider,'value',initial_p_same_slider_value);
set(handles.initial_p_same_threshold,'string',num2str(initial_p_same_slider_value));
centers_of_bins=data_struct.centers_of_bins;

spatial_correlations_distribution=data_struct.spatial_correlations_distribution;
p_same_given_spatial_correlation=data_struct.p_same_given_spatial_correlation;
[~,p_same_ind]=min(abs(initial_p_same_slider_value-(p_same_given_spatial_correlation)));
spatial_correlation_threshold=round(100*centers_of_bins{2}(p_same_ind))/100;
set(handles.correlation_threshold,'string',num2str(spatial_correlation_threshold));
axes(handles.axes4)
cla(handles.axes4,'reset')
plot_p_same_spatial_correlation_slider(spatial_correlations_distribution,p_same_given_spatial_correlation,spatial_correlation_threshold,centers_of_bins)

centroid_distances_distribution=data_struct.centroid_distances_distribution;
p_same_given_centroid_distance=data_struct.p_same_given_centroid_distance;
microns_per_pixel=data_struct.microns_per_pixel;
maximal_distance=data_struct.maximal_distance/microns_per_pixel;
[~,p_same_ind]=min(abs(initial_p_same_slider_value-(p_same_given_centroid_distance)));
centroid_distance_threshold=round(microns_per_pixel*100*centers_of_bins{1}(p_same_ind))/100;
set(handles.distance_threshold,'string',num2str(centroid_distance_threshold));
axes(handles.axes3)
cla(handles.axes3,'reset')
plot_p_same_centroid_distance_slider(centroid_distances_distribution,p_same_given_centroid_distance,centroid_distance_threshold,centers_of_bins,maximal_distance,microns_per_pixel)

handles.data_struct=data_struct;
guidata(hObject, handles)


% --- Executes on slider movement.
function final_p_same_slider_Callback(hObject, eventdata, handles)
% hObject    handle to final_p_same_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

if get(handles.use_model,'Value')==1;
    data_struct=handles.data_struct;
    final_p_same_slider_value=round(100*get(handles.final_p_same_slider,'value'))/100;
    set(handles.decision_thresh,'string',num2str(final_p_same_slider_value));
    centers_of_bins=data_struct.centers_of_bins;
    
    if get(handles.spatial_correlations_2,'Value')==1
        spatial_correlations_distribution=data_struct.spatial_correlations_distribution;
        p_same_given_spatial_correlation=data_struct.p_same_given_spatial_correlation;
        [~,p_same_ind]=min(abs(final_p_same_slider_value-(p_same_given_spatial_correlation)));
        spatial_correlation_threshold=round(100*centers_of_bins{2}(p_same_ind))/100;
        set(handles.correlation_threshold,'string',num2str(spatial_correlation_threshold));
        axes(handles.axes4)
        cla(handles.axes4,'reset')
        plot_p_same_spatial_correlation_slider(spatial_correlations_distribution,p_same_given_spatial_correlation,spatial_correlation_threshold,centers_of_bins)
    else
        centroid_distances_distribution=data_struct.centroid_distances_distribution;
        p_same_given_centroid_distance=data_struct.p_same_given_centroid_distance;
        microns_per_pixel=data_struct.microns_per_pixel;
        maximal_distance=data_struct.maximal_distance/microns_per_pixel;
        [~,p_same_ind]=min(abs(final_p_same_slider_value-(p_same_given_centroid_distance)));
        centroid_distance_threshold=round(microns_per_pixel*100*centers_of_bins{1}(p_same_ind))/100;
        set(handles.distance_threshold,'string',num2str(centroid_distance_threshold));
        axes(handles.axes3)
        cla(handles.axes3,'reset')
        plot_p_same_centroid_distance_slider(centroid_distances_distribution,p_same_given_centroid_distance,centroid_distance_threshold,centers_of_bins,maximal_distance,microns_per_pixel)
    end
    
    handles.data_struct=data_struct;
    guidata(hObject, handles)
end


% --- Executes during object creation, after setting all properties.
function final_p_same_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to final_p_same_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in default_final_registration.
function default_final_registration_Callback(hObject, eventdata, handles)
% hObject    handle to default_final_registration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

if get(handles.use_model,'Value')==1;
    final_p_same_slider_value=0.5;
    data_struct=handles.data_struct;
    set(handles.final_p_same_slider,'value',final_p_same_slider_value);
    set(handles.decision_thresh,'string',num2str(final_p_same_slider_value));
    optimal_threshold=data_struct.spatial_correlation_intersection;
    set(handles.simple_correlation_threshold,'string',optimal_threshold);
    centers_of_bins=data_struct.centers_of_bins;
    
    spatial_correlations_distribution=data_struct.spatial_correlations_distribution;
    p_same_given_spatial_correlation=data_struct.p_same_given_spatial_correlation;
    [~,p_same_ind]=min(abs(final_p_same_slider_value-(p_same_given_spatial_correlation)));
    spatial_correlation_threshold=round(100*centers_of_bins{2}(p_same_ind))/100;
    set(handles.correlation_threshold,'string',num2str(spatial_correlation_threshold));
    axes(handles.axes4)
    cla(handles.axes4,'reset')
    plot_p_same_spatial_correlation_slider(spatial_correlations_distribution,p_same_given_spatial_correlation,spatial_correlation_threshold,centers_of_bins)
    
    centroid_distances_distribution=data_struct.centroid_distances_distribution;
    p_same_given_centroid_distance=data_struct.p_same_given_centroid_distance;
    microns_per_pixel=data_struct.microns_per_pixel;
    maximal_distance=data_struct.maximal_distance/microns_per_pixel;
    [~,p_same_ind]=min(abs(final_p_same_slider_value-(p_same_given_centroid_distance)));
    centroid_distance_threshold=round(microns_per_pixel*100*centers_of_bins{1}(p_same_ind))/100;
    set(handles.distance_threshold,'string',num2str(centroid_distance_threshold));
    axes(handles.axes3)
    cla(handles.axes3,'reset')
    plot_p_same_centroid_distance_slider(centroid_distances_distribution,p_same_given_centroid_distance,centroid_distance_threshold,centers_of_bins,maximal_distance,microns_per_pixel)
    
    handles.data_struct=data_struct;
    guidata(hObject, handles)
end

function decision_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to decision_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of decision_thresh as text
%        str2num(get(hObject,'String')) returns contents of decision_thresh as a double

if get(handles.use_model,'Value')==1;
    data_struct=handles.data_struct;
    final_p_same_slider_value=str2num(get(handles.decision_thresh,'string'));
    set(handles.final_p_same_slider,'value',final_p_same_slider_value);
    centers_of_bins=data_struct.centers_of_bins;
    
    if get(handles.spatial_correlations_2,'Value')==1
        spatial_correlations_distribution=data_struct.spatial_correlations_distribution;
        p_same_given_spatial_correlation=data_struct.p_same_given_spatial_correlation;
        [~,p_same_ind]=min(abs(final_p_same_slider_value-(p_same_given_spatial_correlation)));
        spatial_correlation_threshold=round(100*centers_of_bins{2}(p_same_ind))/100;
        set(handles.correlation_threshold,'string',num2str(spatial_correlation_threshold));
        axes(handles.axes4)
        cla(handles.axes4,'reset')
        plot_p_same_spatial_correlation_slider(spatial_correlations_distribution,p_same_given_spatial_correlation,spatial_correlation_threshold,centers_of_bins)
    else
        centroid_distances_distribution=data_struct.centroid_distances_distribution;
        p_same_given_centroid_distance=data_struct.p_same_given_centroid_distance;
        microns_per_pixel=data_struct.microns_per_pixel;
        maximal_distance=data_struct.maximal_distance/microns_per_pixel;
        [~,p_same_ind]=min(abs(final_p_same_slider_value-(p_same_given_centroid_distance)));
        centroid_distance_threshold=round(microns_per_pixel*100*centers_of_bins{1}(p_same_ind))/100;
        set(handles.distance_threshold,'string',num2str(centroid_distance_threshold));
        axes(handles.axes3)
        cla(handles.axes3,'reset')
        plot_p_same_centroid_distance_slider(centroid_distances_distribution,p_same_given_centroid_distance,centroid_distance_threshold,centers_of_bins,maximal_distance,microns_per_pixel)
    end
    
    handles.data_struct=data_struct;
    guidata(hObject, handles)
end

% --- Executes on button press in use_simple_threshold.
function use_simple_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to use_simple_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_simple_threshold


% --- Executes on button press in use_model.
function use_model_Callback(hObject, eventdata, handles)
% hObject    handle to use_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_model


% --- Executes during object creation, after setting all properties.
function Registration_approach_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Registration_approach (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when selected object is changed in Registration_approach.
function Registration_approach_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in transformation_type
% eventdata  struct with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    struct with handles and user data (see GUIDATA)

if get(handles.use_model,'Value')==1;
    set(handles.use_simple_threshold,'value',0)
    set(handles.decision_thresh,'enable','on')
    set(handles.simple_correlation_threshold,'enable','off')
    set(handles.simple_distance_threshold,'enable','off')
elseif get(handles.use_simple_threshold,'Value')==1;
    set(handles.use_model,'value',0)
    set(handles.decision_thresh,'enable','off')
    if get(handles.spatial_correlations_2,'value')==1
        set(handles.simple_correlation_threshold,'enable','on')
        set(handles.simple_distance_threshold,'enable','off')
    elseif get(handles.centroid_distances_2,'value')==1
        set(handles.simple_correlation_threshold,'enable','off')
        set(handles.simple_distance_threshold,'enable','on')
    end
end


function simple_correlation_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to simple_correlation_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of simple_correlation_threshold as text
%        str2num(get(hObject,'String')) returns contents of simple_correlation_threshold as a double


% --- Executes during object creation, after setting all properties.
function simple_correlation_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to simple_correlation_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when selected object is changed in probability_model_select.
function probability_model_select_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in probability_model_select
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)


if get(handles.use_model,'Value')==1;
    set(handles.use_simple_threshold,'value',0)
    set(handles.decision_thresh,'enable','on')
    set(handles.simple_correlation_threshold,'enable','off')
    set(handles.simple_distance_threshold,'enable','off')
elseif get(handles.use_simple_threshold,'Value')==1;
    set(handles.use_model,'value',0)
    set(handles.decision_thresh,'enable','off')
    if get(handles.spatial_correlations_2,'value')==1
        set(handles.simple_correlation_threshold,'enable','on')
        set(handles.simple_distance_threshold,'enable','off')
    elseif get(handles.centroid_distances_2,'value')==1
        set(handles.simple_correlation_threshold,'enable','off')
        set(handles.simple_distance_threshold,'enable','on')
    end
end


function simple_distance_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to simple_distance_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    struct with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of simple_distance_threshold as text
%        str2num(get(hObject,'String')) returns contents of simple_distance_threshold as a double


% --- Executes during object creation, after setting all properties.
function simple_distance_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to simple_distance_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function distance_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to distance_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of distance_threshold as text
%        str2double(get(hObject,'String')) returns contents of distance_threshold as a double


function correlation_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to correlation_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of correlation_threshold as text
%        str2double(get(hObject,'String')) returns contents of correlation_threshold as a double


% --- Executes on button press in non_rigid.
function non_rigid_Callback(hObject, eventdata, handles)
% hObject    handle to non_rigid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of non_rigid

if get(handles.non_rigid,'Value')==1
    set(handles.maximal_rotation,'enable','off')
    set(handles.transformation_smoothness,'enable','on')
end


function transformation_smoothness_Callback(hObject, eventdata, handles)
% hObject    handle to transformation_smoothness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of transformation_smoothness as text
%        str2double(get(hObject,'String')) returns contents of transformation_smoothness as a double

transformation_smoothness=str2num(get(handles.transformation_smoothness,'string'));
if transformation_smoothness>3 || transformation_smoothness<0.5
    errordlg('FOV smoothing parameter should be between 0.5-3')
    error('FOV smoothing parameter should be between 0.5-3')
end


% --- Executes during object creation, after setting all properties.
function transformation_smoothness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to transformation_smoothness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in translations.
function translations_Callback(hObject, eventdata, handles)
% hObject    handle to translations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of translations

if get(handles.translations,'Value')==1
    set(handles.maximal_rotation,'enable','off')
    set(handles.transformation_smoothness,'enable','off')
end


% --- Executes on button press in translations_rotations.
function translations_rotations_Callback(hObject, eventdata, handles)
% hObject    handle to translations_rotations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of translations_rotations

if get(handles.translations_rotations,'Value')==1
    set(handles.maximal_rotation,'enable','on')
    set(handles.transformation_smoothness,'enable','off')
end
