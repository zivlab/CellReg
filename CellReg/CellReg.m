function varargout = CellReg(varargin)
% This GUI is an implementation of a probabilistic approach for the
% identification of the same neurons (cell registration) across multiple sessions 
% in Ca2+ imaging data, developed by Sheintuch et al., 2016.

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

% Last Modified by GUIDE v2.5 25-Apr-2017 11:50:59

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
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CellReg (see VARARGIN)

% Choose default command line output for CellReg
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CellReg wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Defining the main path of the cell registration scripts:
scriptName=mfilename('fullpath');
[currentpath, ~, ~]=fileparts(scriptName);
cd(currentpath);
addpath(currentpath);

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
sessions_list=cell(1,1);
sessions_list{1}='List of sessions:';
set(handles.list_of_sessions,'string',sessions_list);
set(handles.red_session,'string','1')
set(handles.green_session,'string','2')
set(handles.blue_session,'string','3')
set(handles.decision_thresh,'string','0.5')
set(handles.initial_p_same_slider,'value',0.5);
set(handles.initial_p_same_threshold,'string','0.5');
set(handles.final_p_same_slider,'value',0.5);
set(handles.model_maximal_distance,'string','12')
set(handles.distance_threshold,'string','5')
set(handles.correlation_threshold,'string','0.65')
set(handles.one_photon,'Value',1);
set(handles.translations_rotations,'Value',1);
set(handles.spatial_correlations_2,'Value',1);
set(handles.spatial_correlations,'Value',1);
set(handles.use_joint_model,'Value',0);
set(handles.use_model,'Value',1);
set(handles.x_scale,'string',[])
set(handles.y_scale,'string',[])
set(handles.frame_rate,'string',[])
set(handles.frame_rate,'value',0);
set(handles.pixel_to_mic,'string',[])
set(handles.x_scale,'value',0)
set(handles.y_scale,'value',0)
set(handles.pixel_to_mic,'value',0)
set(handles.frame_rate,'value',0)
set(handles.reference_session,'string','1')
set(handles.maximal_rotation','string','30')
set(handles.maximal_rotation','enable','on')
set(handles.distance_threshold,'enable','off')
set(handles.correlation_threshold,'enable','on')
set(handles.comments,'string',[])
handles.data_struct=data_struct;
guidata(hObject, handles)


% --- Outputs from this function are returned to the command line.
function varargout = CellReg_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function load_new_data_Callback(hObject,~, handles)
% hObject    handle to load_new_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Stage 1: Loading the spatial footprints of cellular activity from the different
% sessions.

% This callback loads a new data set which includes several sessions with
% their spatial footprints, centroid locations (optional), and events (optional). A
% single folder should be selected with all the mat files with number:
% example: "finalFiltersMat_1", "finalFiltersMat_2", "finalEventsMat_1", "finalEventsMat_2"

data_struct=handles.data_struct;
msgbox('Please choose the folder containing the spatial footprints from all the sessions: ')
pause(3)
sessions_dir=uigetdir();
cd(sessions_dir)
list=dir(sessions_dir);
a=1:numel(list);
namesOfFiles = arrayfun(@(x) list(x).name,a,'UniformOutput',false);
tifFilesIdx = find(cellfun (@(x) ~isempty(strfind(x,'Filters')) ,namesOfFiles));
if isempty(tifFilesIdx)
    errordlg('This folder does not contain data with the required format')
else
    num_sessions=length(tifFilesIdx);
    tifFilesIdx_order=zeros(1,length(tifFilesIdx));

    button_cent = questdlg('Does your data include centroid locations?','title','Yes','No','Yes');
    
    button = questdlg('Does your data include the events?','title','Yes','No','Yes');
    
    all_centroids=cell(1,num_sessions);
    all_filters=cell(1,num_sessions);
    if strcmp(button,'Yes')
        all_events=cell(1,num_sessions);
    end
    for n=1:num_sessions
        is_data=0;
        tifFilesIdx = find(cellfun (@(x) ~isempty(strfind(x,['FiltersMat_' num2str(n)])) ,namesOfFiles));
        if isempty(tifFilesIdx)
            tifFilesIdx = find(cellfun (@(x) ~isempty(strfind(x,['Filters_' num2str(n)])) ,namesOfFiles));
            if isempty(tifFilesIdx)
                errordlg('This folder does not contain spatial footprints with the required format')
            else
                is_data=1;                
                tifFilesIdx_order(n)=tifFilesIdx(1);
            end
        else
            is_data=1;
            tifFilesIdx_order(n)=tifFilesIdx(1);
        end
        if is_data==1;
            data_string=namesOfFiles(tifFilesIdx);
            temp_data=load(data_string{1});
            if isstruct(temp_data)
                field_name=fieldnames(temp_data);
                all_filters{n}=getfield(temp_data,field_name{1});
            else
                all_filters{n}=temp_data;
            end
        end
        if strcmp(button_cent,'Yes')
            is_data=0;
            tifFilesIdx = find(cellfun (@(x) ~isempty(strfind(x,['LocationsMat_' num2str(n)])) ,namesOfFiles));
            if isempty(tifFilesIdx)
                tifFilesIdx = find(cellfun (@(x) ~isempty(strfind(x,['Locations_' num2str(n)])) ,namesOfFiles));
                if isempty(tifFilesIdx)
                    errordlg('This folder does not contain centroid locations with the required format')
                else
                    is_data=1;
                end
            else
                is_data=1;
            end
            if is_data==1;
                data_string=namesOfFiles(tifFilesIdx);
                temp_data=load(data_string{1});
                if isstruct(temp_data)
                    field_name=fieldnames(temp_data);
                    all_centroids{n}=getfield(temp_data,field_name{1});
                else
                    all_centroids{n}=temp_data;
                end
                if size(all_centroids{n},1)~=size(all_filters{n},1);
                    errordlg(['There is an inconsistency in the number of centroids and spatial footprints in session number ' num2str(n)]);
                end
            end
        end
        
        if strcmp(button,'Yes')
            frame_rate=str2num(get(handles.frame_rate,'string'));
            if isempty(frame_rate)
                msgbox('Please insert the frame rate and press enter ')
                waitfor(handles.frame_rate,'value',1);
                frame_rate=str2double(get(handles.frame_rate,'string'));
            end
            data_struct.frame_rate=frame_rate;
            is_data=0;
            tifFilesIdx = find(cellfun (@(x) ~isempty(strfind(x,['EventsMat_' num2str(n)])) ,namesOfFiles));
            if isempty(tifFilesIdx)
                tifFilesIdx = find(cellfun (@(x) ~isempty(strfind(x,['Events_' num2str(n)])) ,namesOfFiles));
                if isempty(tifFilesIdx)
                    errordlg('This folder does not contain events with the required format')
                else
                    is_data=1;
                end
            else
                is_data=1;
            end
            if is_data==1;
                data_string=namesOfFiles(tifFilesIdx);
                temp_data=load(data_string{1});
                if isstruct(temp_data)
                    field_name=fieldnames(temp_data);
                    all_events{n}=getfield(temp_data,field_name{1});
                else
                    all_events{n}=temp_data;
                end
                if size(all_events{n},2)~=size(all_filters{n},1);
                    errordlg(['There is an inconsistency in the number of cell-events and spatial footprints in session number ' num2str(n)]);
                end
            end
        end
    end

    center(1)=size(all_filters{n},2)/2;
    center(2)=size(all_filters{n},3)/2;
    
    pixel_to_mic=str2num(get(handles.pixel_to_mic,'string'));
    if isempty(pixel_to_mic)
        msgbox('Please insert the pixel size in microns and press enter ')
        waitfor(handles.pixel_to_mic,'value',1);
        pixel_to_mic=str2double(get(handles.pixel_to_mic,'string'));
    end
    image_size=[size(all_filters{1},2) size(all_filters{1},3)];
    x_scale=image_size(2)*pixel_to_mic;
    y_scale=image_size(1)*pixel_to_mic;
    
    set(handles.pixel_to_mic,'string',num2str(round(100*pixel_to_mic)/100));
    set(handles.x_scale,'string',num2str(x_scale));
    set(handles.y_scale,'string',num2str(y_scale));
    
    msgbox('Please select the folder in which the results will be saved')
    pause(3)
    results_dir = uigetdir; % the directory which the final results will be saved
    data_struct.results_dir=results_dir;
    cd(results_dir)
    mkdir('Figures');
    figures_dir=fullfile(results_dir,'Figures');
    data_struct.figures_dir=figures_dir;
    
    msgbox('Please wait a few seconds while the spatial footprints are being loaded');    
    data_struct.button_events=button;
    data_struct.button_cent=button_cent;
    data_struct.pixel_to_mic=pixel_to_mic;
    data_struct.x_scale=x_scale;
    data_struct.y_scale=y_scale;
    if strcmp(button_cent,'Yes')
        data_struct.all_centroids=all_centroids;
    end
    data_struct.all_filters=all_filters;
    if strcmp(button,'Yes')
        data_struct.all_events=all_events;
    end
    data_struct.center=center;
    data_struct.image_size=image_size;
    data_struct.num_sessions=num_sessions;
    data_struct.sessions_dir=sessions_dir;
    data_struct.sessions_list=cell(1,num_sessions+1);
    data_struct.sessions_list{1}='List of sessions:';
    tifFilesIdx = find(cellfun (@(x) ~isempty(strfind(x,'Filters')) ,namesOfFiles));
    for n=1:num_sessions
        this_file_name=namesOfFiles{tifFilesIdx_order(n)};
        data_struct.sessions_list{n+1}=['Session ' num2str(n) ' - ' this_file_name];
    end
    
    % Calculating all footprints projections:
    all_projections=cell(1,num_sessions);
    for n=1:num_sessions
        this_session_filters=all_filters{n};
        num_filters=size(this_session_filters,1);
        normalized_filters=zeros(size(this_session_filters));
        for k=1:num_filters
            normalized_filters(k,:,:)=this_session_filters(k,:,:)/max(max(this_session_filters(k,:,:)));
        end
        all_projections{n}=squeeze(sum(normalized_filters,1));
    end
    
    subx=4;
    suby=ceil(num_sessions/subx);
    if num_sessions>4
        figure('units','normalized','outerposition',[0.1 0.04 0.8 0.96])
        for n=1:num_sessions
            subplot(suby,subx,n)
            imagesc(all_projections{n},[0 2])
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            colormap('gray')
            this_file_name=namesOfFiles{tifFilesIdx_order(n)};
            title(['Session ' num2str(n) ' - ' strrep(this_file_name,'_','\_')],'fontsize',14,'fontweight','bold')
        end
    else
        figure('units','normalized','outerposition',[0.1 0.2 0.8 0.5])
        for n=1:num_sessions
            subplot(1,num_sessions,n)
            imagesc(all_projections{n},[0 2])
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            colormap('gray')
            this_file_name=namesOfFiles{tifFilesIdx_order(n)};
            title(['Session ' num2str(n) ' - ' strrep(this_file_name,'_','\_')],'fontsize',14,'fontweight','bold')
        end
    end    
    data_struct.button_cent=button_cent;
    set(handles.list_of_sessions,'string',data_struct.sessions_list);
    handles.data_struct=data_struct;
    guidata(hObject, handles)
    msgbox('Finished loading sessions')
end


% --- Executes on button press in add_session.
function add_session_Callback(hObject,~, handles)
% hObject    handle to add_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This callback adds another session to the list of sessions to be
% registered. The folder containg the filters, centroids, and events
% (optional) should be selected

data_struct=handles.data_struct;
if isfield(data_struct,'num_sessions')
    all_centroids=data_struct.all_centroids;
    all_filters=data_struct.all_filters;
    if strcmp(data_struct.button_events,'Yes')
        all_events=data_struct.all_events;
    end
    msgbox('Please choose the folder containing the spatial footprints for this session: ')
    pause(1)
    sessions_dir=uigetdir();
    cd(sessions_dir)
    list=dir(sessions_dir);
    a=1:numel(list);
    namesOfFiles = arrayfun(@(x) list(x).name,a,'UniformOutput',false);
    msgbox('Please wait a few seconds while the spatial footprints are being loaded');    
    tifFilesIdx = find(cellfun (@(x) ~isempty(strfind(x,'Filters')) ,namesOfFiles));
    if length(tifFilesIdx)>1
        errordlg('This folder contains data from more than one session')
    elseif isempty(tifFilesIdx)
        errordlg('This folder does not contain data with the required format')
    else
        
        tifFilesIdx = cellfun (@(x) ~isempty(strfind(x,'Filters')) ,namesOfFiles);
        if isempty(tifFilesIdx)
            errordlg('This folder does not contain spatial footprints with the required format')
        else
            data_string=namesOfFiles(tifFilesIdx);
            temp_data=load(data_string{1});
            if isstruct(temp_data)
                field_name=fieldnames(temp_data);
                all_filters{data_struct.num_sessions+1}=getfield(temp_data,field_name{1});
            else
                all_filters{data_struct.num_sessions+1}=temp_data;
            end
        end
        
        if strcmp(data_struct.button_cent,'Yes')
            tifFilesIdx = cellfun (@(x) ~isempty(strfind(x,'Locations')) ,namesOfFiles);
            if isempty(tifFilesIdx)
                errordlg('This folder does not contain centroid locations with the required format')
            else
                locations_string=namesOfFiles(tifFilesIdx);                
                temp_data=load(locations_string{1});
                if isstruct(temp_data)
                    field_name=fieldnames(temp_data);
                    all_centroids{data_struct.num_sessions+1}=getfield(temp_data,field_name{1});
                else
                    all_centroids{data_struct.num_sessions+1}=temp_data;
                end
                if size(all_centroids{data_struct.num_sessions+1},1)~=size(all_filters{data_struct.num_sessions+1},1);
                    errordlg(['There is an inconsistency in the number of centroids and spatial footprints in session number ' num2str(n)]);
                end
            end
        end
        
        % Calculating all footprints projections:
        this_session_filters=all_filters{data_struct.num_sessions+1};
        num_filters=size(this_session_filters,1);
        normalized_filters=zeros(size(this_session_filters));
        for k=1:num_filters
            normalized_filters(k,:,:)=this_session_filters(k,:,:)/max(max(this_session_filters(k,:,:)));
        end
        all_projections=squeeze(sum(normalized_filters,1));        
        
        data_struct.all_centroids=all_centroids;
        data_struct.all_filters=all_filters;
        if strcmp(data_struct.button_events,'Yes')
            tifFilesIdx = cellfun (@(x) ~isempty(strfind(x,'Events')) ,namesOfFiles);
            if isempty(tifFilesIdx)
                errordlg('This folder does not contain events with the required format')
            else
                events_string=namesOfFiles(tifFilesIdx);                
                temp_data=load(events_string{1});
                if isstruct(temp_data)
                    field_name=fieldnames(temp_data);
                    all_events{data_struct.num_sessions+1}=getfield(temp_data,field_name{1});
                else
                    all_events{data_struct.num_sessions+1}=temp_data;
                end
                if size(all_events{data_struct.num_sessions+1},2)~=size(all_filters{data_struct.num_sessions+1},1);
                    errordlg(['There is an inconsistency in the number of cell-events and spatial footprints in session number ' num2str(n)]);
                end   
                data_struct.all_events=all_events;
            end
        end
        
        data_struct.num_sessions=data_struct.num_sessions+1;
        sessions_list=data_struct.sessions_list;        
        directory_ind=find(sessions_dir=='\');
        directory_name=sessions_dir(directory_ind(end)+1:end);
        if ~isempty(strfind(directory_name,'finalResults'))
            directory_name=sessions_dir(directory_ind(end-1)+1:directory_ind(end)-1);
        end
        sessions_list{data_struct.num_sessions+1}=['Session ' num2str(data_struct.num_sessions) ' - ' directory_name];        
        data_struct.sessions_list=sessions_list;
        set(handles.list_of_sessions,'string',data_struct.sessions_list);
        
        figure
        imagesc(all_projections,[0 2])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        colormap('gray')
        title(['Session ' num2str(data_struct.num_sessions) ' - ' strrep(directory_name,'_','\_')],'fontsize',18,'fontweight','bold')
        set(gcf,'PaperPositionMode','auto')                
    end
else
    msgbox('Please choose the folder containing the spatial footprints for this session: ')
    pause(3)
    sessions_dir=uigetdir();
    cd(sessions_dir)
    list=dir(sessions_dir);
    a=1:numel(list);
    namesOfFiles = arrayfun(@(x) list(x).name,a,'UniformOutput',false);
    tifFilesIdx = find(cellfun (@(x) ~isempty(strfind(x,'Filters')) ,namesOfFiles));
    if length(tifFilesIdx)>1
        errordlg('This folder contains data from more than one session')
    elseif isempty(tifFilesIdx)
        errordlg('This folder does not contain data with the required format')
    else
        num_sessions=1;
        button_cent = questdlg('Does your data include centroid locations?','title','Yes','No','Yes');
        button = questdlg('Does your data include the events?','title','Yes','No','Yes');
        
        all_centroids=cell(1,num_sessions);
        all_filters=cell(1,num_sessions);
        if strcmp(button,'Yes')
            all_events=cell(1,num_sessions);
            frame_rate=str2num(get(handles.frame_rate,'string'));
            if isempty(frame_rate)
                msgbox('Please insert the frame rate and press enter ')
                waitfor(handles.frame_rate,'value',1);
                frame_rate=str2double(get(handles.frame_rate,'string'));
            end
            data_struct.frame_rate=frame_rate;
        end
        msgbox('Please wait a few seconds while the spatial footprints are being loaded');       

        tifFilesIdx = cellfun (@(x) ~isempty(strfind(x,'Filters')) ,namesOfFiles);
        if isempty(tifFilesIdx)
            errordlg('This folder does not contain spatial footprints with the required format')
        else
            data_string=namesOfFiles(tifFilesIdx);            
            temp_data=load(data_string{1});
            if isstruct(temp_data)
                field_name=fieldnames(temp_data);
                all_filters{1}=getfield(temp_data,field_name{1});
            else
                all_filters{1}=temp_data;
            end            
        end
        
        if strcmp(button_cent,'Yes')
            tifFilesIdx = cellfun (@(x) ~isempty(strfind(x,'Locations')) ,namesOfFiles);
            locations_string=namesOfFiles(tifFilesIdx);
            if isempty(tifFilesIdx)
                errordlg('This folder does not contain centroid locations with the required format')
            else                
                temp_data=load(locations_string{1});
                if isstruct(temp_data)
                    field_name=fieldnames(temp_data);
                    all_centroids{1}=getfield(temp_data,field_name{1});
                else
                    all_centroids{1}=temp_data;
                end
                if size(all_centroids{1},1)~=size(all_filters{1},1);
                    errordlg(['There is an inconsistency in the number of centroids and spatial footprints in session number ' num2str(n)]);
                end                
            end
        end
        
        % Calculating all footprints projections:
        this_session_filters=all_filters{1};
        num_filters=size(this_session_filters,1);
        normalized_filters=zeros(size(this_session_filters));
        for k=1:num_filters
            normalized_filters(k,:,:)=this_session_filters(k,:,:)/max(max(this_session_filters(k,:,:)));
        end
        all_projections=squeeze(sum(normalized_filters,1));        
        
        if strcmp(button,'Yes')
            tifFilesIdx = cellfun (@(x) ~isempty(strfind(x,'Events')) ,namesOfFiles);
            if isempty(tifFilesIdx)
                errordlg('This folder does not contain events with the required format')
            else
                events_string=namesOfFiles(tifFilesIdx);
                temp_data=load(events_string{1});
                if isstruct(temp_data)
                    field_name=fieldnames(temp_data);
                    all_events{1}=getfield(temp_data,field_name{1});
                else
                    all_events{1}=temp_data;
                end
                if size(all_events{1},2)~=size(all_filters{1},1);
                    errordlg(['There is an inconsistency in the number of cell-events and spatial footprints in session number ' num2str(n)]);
                end         
            end
        end
        
        center(1)=size(all_filters{1},2)/2;
        center(2)=size(all_filters{1},3)/2;
            
        pixel_to_mic=str2num(get(handles.pixel_to_mic,'string'));
        if isempty(pixel_to_mic)
            msgbox('Please insert the pixel size in microns and press enter ')
            waitfor(handles.pixel_to_mic,'value',1);
            pixel_to_mic=str2double(get(handles.pixel_to_mic,'string'));
        end
        image_size=[size(all_filters{1},2) size(all_filters{1},3)];
        x_scale=image_size(2)*pixel_to_mic;
        y_scale=image_size(1)*pixel_to_mic;
                
        set(handles.pixel_to_mic,'string',num2str(round(100*pixel_to_mic)/100));
        set(handles.x_scale,'string',num2str(x_scale));
        set(handles.y_scale,'string',num2str(y_scale));
        
        msgbox('Please select the folder in which the results will be saved')
        pause(3)
        results_dir = uigetdir; % the directory which the final results will be saved
        data_struct.results_dir=results_dir;
        cd(results_dir)
        mkdir('Figures');
        figures_dir=fullfile(results_dir,'Figures');
        data_struct.figures_dir=figures_dir;
        
        data_struct.button_events=button;
        data_struct.pixel_to_mic=pixel_to_mic;
        data_struct.x_scale=x_scale;
        data_struct.y_scale=y_scale;
        data_struct.all_centroids=all_centroids;
        data_struct.all_filters=all_filters;
        if strcmp(button,'Yes')
            data_struct.all_events=all_events;
        end
        data_struct.button_cent=button_cent;
        data_struct.center=center;
        data_struct.image_size=image_size;
        data_struct.num_sessions=num_sessions;
        data_struct.sessions_dir=sessions_dir;
        data_struct.sessions_list=cell(1,num_sessions+1);
        data_struct.sessions_list{1}='List of sessions:';
        directory_ind=find(sessions_dir=='\');
        directory_name=sessions_dir(directory_ind(end)+1:end);
        if ~isempty(strfind(directory_name,'finalResults'))
            directory_name=sessions_dir(directory_ind(end-1)+1:directory_ind(end)-1);
        end        
        data_struct.sessions_list{num_sessions+1}=['Session 1 - ' directory_name];
        set(handles.list_of_sessions,'string',data_struct.sessions_list);
        
        figure
        imagesc(all_projections,[0 2])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        colormap('gray')
        title(['Session ' num2str(data_struct.num_sessions) ' - ' strrep(directory_name,'_','\_')],'fontsize',18,'fontweight','bold')
        set(gcf,'PaperPositionMode','auto')        
    end
end
handles.data_struct=data_struct;
guidata(hObject, handles)
msgbox('Finished loading session')


% --- Executes on button press in remove_session.
function remove_session_Callback(hObject,~, handles)
% hObject    handle to remove_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This callback removes the selected session from the list of sessions to be registered

data_struct=handles.data_struct;
if isfield(data_struct,'num_sessions')
    num_sessions=data_struct.num_sessions;
    all_centroids=data_struct.all_centroids;
    all_filters=data_struct.all_filters;
    if strcmp(data_struct.button_events,'Yes')
        all_events=data_struct.all_events;
    end
    chosen_session=get(handles.list_of_sessions,'value')-1;
    if chosen_session==0
        errordlg('Please choose a session to remove')
    end
    sessions_to_keep=setdiff(1:num_sessions,chosen_session);
    set(handles.list_of_sessions,'value',1)
    temp_centroids=all_centroids(sessions_to_keep);
    all_centroids=temp_centroids;
    temp_filters=all_filters(sessions_to_keep);
    all_filters=temp_filters;    
    if strcmp(data_struct.button_events,'Yes')
        temp_events=all_events(sessions_to_keep);
        all_events=temp_events;
    end  
    
    if isfield(data_struct,'all_projections_corrected')
        temp_projections=data_struct.all_projections_corrected(sessions_to_keep);
        data_struct.all_projections_corrected=temp_projections;
        temp_centroids_corrected=data_struct.all_centroids_corrected(sessions_to_keep);
        data_struct.all_centroids_corrected=temp_centroids_corrected;
        temp_filters_corrected=data_struct.all_filters_corrected(sessions_to_keep);
        data_struct.all_filters_corrected=temp_filters_corrected;
        temp_projections=data_struct.all_projections(sessions_to_keep);
        data_struct.all_projections=temp_projections;
    end
    num_sessions=num_sessions-1;
    data_struct.all_centroids=all_centroids;
    data_struct.all_filters=all_filters;
    data_struct.num_sessions=num_sessions;
    if strcmp(data_struct.button_events,'Yes')
        data_struct.all_events=all_events;
    end
    original_names=data_struct.sessions_list;
    namesIdx = find(cellfun (@(x) ~isempty(strfind(x,'Session')) ,original_names));
    data_struct.sessions_list=cell(1,data_struct.num_sessions+1);
    data_struct.sessions_list{1}='List of sessions:';    
    for n=1:num_sessions
        this_string=original_names{namesIdx(sessions_to_keep(n))};              
        if n<10
            this_file_name=this_string(13:end);
        else
            this_file_name=this_string(14:end);
        end
        data_struct.sessions_list{n+1}=['Session ' num2str(n) ' - ' this_file_name];
    end
    
    set(handles.list_of_sessions,'string',data_struct.sessions_list);
    handles.data_struct=data_struct;
    guidata(hObject, handles)
    if chosen_session>0
        msgbox('Finished removing session')
    end
end


% --------------------------------------------------------------------
function load_transformed_data_Callback(hObject,~, handles)
% hObject    handle to load_transformed_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This callback loads sessions that were already transformed into a reference coordinate system.
% For such data the compute model should be the next step.

msgbox('Please choose the folder containing the transformed sessions: ')
pause(3)
sessions_dir=uigetdir();
cd(sessions_dir)
list=dir(sessions_dir);
msgbox('Please select the folder in which the results will be saved')
pause(3)
results_dir=uigetdir; % the directory which the final results will be saved
data_struct.results_dir=results_dir;
a=1:numel(list);
namesOfFiles = arrayfun(@(x) list(x).name,a,'UniformOutput',false);
tifFilesIdx = cellfun (@(x) ~isempty(strfind(x,'transformed_data_struct')) ,namesOfFiles);
fullfilename=fullfile(sessions_dir,namesOfFiles(tifFilesIdx));
msgbox('Please wait a few seconds while the transformed data is being loaded');
if ~isempty(fullfilename)
    load(fullfilename{1});
else
    tifFilesIdx = cellfun (@(x) ~isempty(strfind(x,'registered_data_struct')) ,namesOfFiles);
    fullfilename=fullfile(sessions_dir,namesOfFiles(tifFilesIdx));
    if ~isempty(fullfilename)
        load(fullfilename{1});
    else
        errordlg('This folder does not contain data with the required format')
    end
end
data_struct.results_dir=results_dir;

all_projections=data_struct.all_projections;
if isfield(data_struct,'all_projections_corrected')
    all_projections_corrected=data_struct.all_projections_corrected;
else
    all_projections_corrected=data_struct.all_projections;
end
num_sessions=data_struct.num_sessions;
max_x=data_struct.max_x;
max_y=data_struct.max_y;

rgb_ind=[1 2 3];
all_projections_corrected_rgb=zeros(max_y,max_x,3);
all_projections_corrected_rgb(:,:,1)=all_projections_corrected{rgb_ind(1)};
all_projections_corrected_rgb(:,:,2)=all_projections_corrected{rgb_ind(2)};
if num_sessions>2
    all_projections_corrected_rgb(:,:,3)=all_projections_corrected{rgb_ind(3)};
end
all_projections_corrected_rgb(all_projections_corrected_rgb>1)=1;
axes(handles.axes1);
imagesc(all_projections_corrected_rgb)
title('RGB overlay','FontWeight','Bold','fontsize',16)
set(gca,'xtick',[])
set(gca,'ytick',[])

subx=4;
suby=ceil(num_sessions/subx);
if num_sessions>4
    figure('units','normalized','outerposition',[0.1 0.04 0.8 0.96])
    for n=1:num_sessions
        subplot(suby,subx,n)
        imagesc(all_projections{n},[0 2])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        colormap('gray')
        this_session_name=data_struct.sessions_list{n+1};
        title(strrep(this_session_name,'_','\_'),'fontsize',14,'fontweight','bold')
    end
else
    figure('units','normalized','outerposition',[0.1 0.2 0.8 0.5])
    for n=1:num_sessions
        subplot(1,num_sessions,n)
        imagesc(all_projections{n},[0 2])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        colormap('gray')
        this_session_name=data_struct.sessions_list{n+1};
        title(strrep(this_session_name,'_','\_'),'fontsize',14,'fontweight','bold')
    end
end
    
set(handles.pixel_to_mic,'string',num2str(round(100*data_struct.pixel_to_mic)/100));
set(handles.x_scale,'string',num2str(data_struct.x_scale));
set(handles.y_scale,'string',num2str(data_struct.y_scale));
set(handles.reference_session,'string',num2str(data_struct.reference_session))
set(handles.list_of_sessions,'string',data_struct.sessions_list)

if isfield(data_struct,'frame_rate')
    set(handles.frame_rate,'string',num2str(data_struct.frame_rate))
end

handles.data_struct=data_struct;
guidata(hObject, handles)
msgbox('Finished loading transformed sessions')


% --------------------------------------------------------------------
function load_modeled_data_Callback(hObject,~, handles)
% hObject    handle to load_modeled_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This callback loads sessions that were already modeled. For such data
% the register cells initial should be the next step

msgbox('Please choose the folder containing the modeled data: ')
pause(3)
sessions_dir=uigetdir();
cd(sessions_dir)
list=dir(sessions_dir);
msgbox('Please select the folder in which the results will be saved')
pause(3)
results_dir = uigetdir; % the directory which the final results will be saved
data_struct.results_dir=results_dir;
a=1:numel(list);
namesOfFiles = arrayfun(@(x) list(x).name,a,'UniformOutput',false);
tifFilesIdx = cellfun (@(x) ~isempty(strfind(x,'model_data_struct')) ,namesOfFiles);
fullfilename=fullfile(sessions_dir,namesOfFiles(tifFilesIdx));
msgbox('Please wait a few seconds while the modeled data is being loaded');
if ~isempty(fullfilename);
    load(fullfilename{1})
else
    errordlg('This folder does not contain data with the required format')
end

data_struct.results_dir=results_dir;
set(handles.pixel_to_mic,'string',num2str(round(100*data_struct.pixel_to_mic)/100));
set(handles.x_scale,'string',num2str(data_struct.x_scale));
set(handles.y_scale,'string',num2str(data_struct.y_scale));
set(handles.reference_session,'string',num2str(data_struct.reference_session))
if isfield(data_struct,'thresh_dist_from_intersection')
    set(handles.distance_threshold,'string',num2str(data_struct.thresh_dist_from_intersection))
end
if isfield(data_struct,'thresh_corr_from_intersection')
    set(handles.correlation_threshold,'string',num2str(data_struct.thresh_corr_from_intersection))
end

if isfield(data_struct,'frame_rate')
    set(handles.frame_rate,'string',num2str(data_struct.frame_rate))
end

all_projections=data_struct.all_projections;
if isfield(data_struct,'all_projections_corrected')
    all_projections_corrected=data_struct.all_projections_corrected;
else
    all_projections_corrected=data_struct.all_projections;
end
num_sessions=data_struct.num_sessions;
max_x=data_struct.max_x;
max_y=data_struct.max_y;

rgb_ind=[1 2 3];
all_projections_corrected_rgb=zeros(max_y,max_x,3);
all_projections_corrected_rgb(:,:,1)=all_projections_corrected{rgb_ind(1)};
all_projections_corrected_rgb(:,:,2)=all_projections_corrected{rgb_ind(2)};
if num_sessions>2
    all_projections_corrected_rgb(:,:,3)=all_projections_corrected{rgb_ind(3)};
end

if isfield(data_struct,'overlapping_matrix')
    overlapping_matrix=data_struct.overlapping_matrix;
    all_projections_corrected_rgb=all_projections_corrected_rgb+0.25*repmat(overlapping_matrix,1,1,3);
end
all_projections_corrected_rgb(all_projections_corrected_rgb>1)=1;

axes(handles.axes1);
imagesc(all_projections_corrected_rgb)
title('RGB overlay','FontWeight','Bold','fontsize',16)
set(gca,'xtick',[])
set(gca,'ytick',[])

subx=4;
suby=ceil(num_sessions/subx);
if num_sessions>4
    figure('units','normalized','outerposition',[0.1 0.04 0.8 0.96])
    for n=1:num_sessions
        subplot(suby,subx,n)
        imagesc(all_projections{n},[0 2])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        colormap('gray')
        this_session_name=data_struct.sessions_list{n+1};
        title(this_session_name,'fontsize',14,'fontweight','bold')
    end
else
    figure('units','normalized','outerposition',[0.1 0.2 0.8 0.5])
    for n=1:num_sessions
        subplot(1,num_sessions,n)
        imagesc(all_projections{n},[0 2])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        colormap('gray')
        this_session_name=data_struct.sessions_list{n+1};
        title(this_session_name,'fontsize',14,'fontweight','bold')
    end
end
set(handles.list_of_sessions,'string',data_struct.sessions_list)

initial_p_same_slider_value=0.5;
set(handles.initial_p_same_slider,'value',initial_p_same_slider_value);
set(handles.initial_p_same_threshold,'string',num2str(initial_p_same_slider_value));

ctrs=data_struct.ctrs;
pixel_to_mic=data_struct.pixel_to_mic;
n_corr=data_struct.n_corr;
normalized_distance=data_struct.normalized_distance;
temp_dist_thresh=str2double(get(handles.model_maximal_distance,'string'))/pixel_to_mic;

p_value_corr=data_struct.p_value_corr;
[~,p_same_ind_corr]=min(abs(initial_p_same_slider_value-(1-p_value_corr)));
p_value_dist=data_struct.p_value_dist;
[~,p_same_ind_dist]=min(abs(initial_p_same_slider_value-(1-p_value_dist)));

corr_thresh=round(100*ctrs{2}(p_same_ind_corr))/100;
set(handles.correlation_threshold,'string',num2str(corr_thresh));
dist_thresh=round(pixel_to_mic*100*ctrs{1}(p_same_ind_dist))/100;
set(handles.distance_threshold,'string',num2str(dist_thresh));

axes(handles.axes4)
plot([],[])
set(gca,'xtick',[])
set(gca,'ytick',[])
cla(handles.axes4,'reset')
start_x=0;
end_x=1;
y_vec=repmat(n_corr,[2 1]);
y_vec=y_vec(:);
x_vec=(ctrs{2}(2:end)+ctrs{2}(1:end-1))/2;
x_vec=repmat(x_vec,[2 1]);
x_vec=[start_x; x_vec(:); end_x];
for run_bins=1:length(x_vec)/2
    current_color=p_value_corr(run_bins)*[1 1 1];
    patch(x_vec([1 1 2 2]+2*(run_bins-1)),[0 [1 1]*n_corr(run_bins) 0],current_color,'EdgeColor',current_color)
    hold on
end
plot(x_vec,y_vec,'k-','linewidth',2)
xlabel('Spatial correlation','fontsize',12,'fontweight','bold')
x_label=linspace(0,1,6);
x=linspace(0,1,6);
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
xlim([0 1])
[~,ind_005]=min(abs(0.05-(1-p_value_corr)));
p_005=ctrs{2}(ind_005);
[~,ind_05]=min(abs(0.5-(1-p_value_corr)));
p_05=ctrs{2}(ind_05);
[~,ind_095]=min(abs(0.95-(1-p_value_corr)));
p_095=ctrs{2}(ind_095);
hold on
plot([p_005 p_005],[0 max(n_corr)],'--','linewidth',3,'color','k')
hold on
plot([p_05 p_05],[0 max(n_corr)],'--','linewidth',3,'color','k')
hold on
plot([p_095 p_095],[0 max(n_corr)],'--','linewidth',3,'color','k')
hold on
plot([corr_thresh corr_thresh],[0 max(n_corr)],'linewidth',2,'color','r')
text(p_005,1.1*max(n_corr),'0.05','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
text(p_095,1.1*max(n_corr),'0.95','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
text(p_05,1.1*max(n_corr),'0.5','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
axes(handles.axes3)
plot([],[])
set(gca,'xtick',[])
set(gca,'ytick',[])
cla(handles.axes3,'reset')
start_x=0;
end_x=temp_dist_thresh*pixel_to_mic;
y_vec=repmat(normalized_distance,[2 1]);
y_vec=y_vec(:);
x_vec=(pixel_to_mic*ctrs{1}(2:end)+pixel_to_mic*ctrs{1}(1:end-1))/2;
x_vec=repmat(x_vec,[2 1]);
x_vec=[start_x; x_vec(:); end_x];
for run_bins=1:length(x_vec)/2
    current_color=p_value_dist(run_bins)*[1 1 1];
    patch(x_vec([1 1 2 2]+2*(run_bins-1)),[0 [1 1]*normalized_distance(run_bins) 0],current_color,'EdgeColor',current_color)
    hold on
end
plot(x_vec,y_vec,'k-','linewidth',2)
xlim([0 pixel_to_mic*temp_dist_thresh])
xlabel('Centroids distance (\mum)','FontWeight','Bold','fontsize',12)
x_label=0:3:pixel_to_mic*temp_dist_thresh;
x=0:3:pixel_to_mic*temp_dist_thresh;
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
[~,ind_005]=min(abs(0.05-(1-p_value_dist)));
p_005=pixel_to_mic*ctrs{1}(ind_005);
[~,ind_05]=min(abs(0.5-(1-p_value_dist)));
p_05=pixel_to_mic*ctrs{1}(ind_05);
[~,ind_095]=min(abs(0.95-(1-p_value_dist)));
p_095=pixel_to_mic*ctrs{1}(ind_095);
hold on
plot([p_005 p_005],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
hold on
plot([p_05 p_05],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
hold on
plot([p_095 p_095],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
hold on
plot([dist_thresh dist_thresh],[0 max(normalized_distance)],'linewidth',2,'color','r')
text(p_005,1.1*max(normalized_distance),'0.05','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
text(p_095,1.1*max(normalized_distance),'0.95','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
text(p_05,1.1*max(normalized_distance),'0.5','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')

if isfield(data_struct,'log_grid')
    log_grid=data_struct.log_grid;
    p_value=data_struct.p_value;
    num_bins=data_struct.number_of_bins;
    axes(handles.axes5)
    imagesc(log_grid)
    colormap('jet')
    y=round(linspace(1,num_bins,5));
    y_label=round(10*linspace(pixel_to_mic*temp_dist_thresh,0,5))/10;
    x=round(linspace(1,num_bins,6));
    x_label=linspace(0,1,6);
    set(gca,'YTick',y)
    set(gca,'YTickLabel',y_label,'fontsize',14,'FontWeight','Bold')
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',14,'FontWeight','Bold')
    xlabel('Spatial correlation','FontWeight','Bold','fontsize',14)
    ylabel('Centroids distance (\mum)','FontWeight','Bold','fontsize',14)
    freezeColors
    hold on
    contour_values=0.05:0.05:0.95;
    contour(p_value,contour_values,'linewidth',1)
    colormap('gray')
    hold on
    contour(p_value,[0.05 0.5 0.95],'linewidth',3)
    hold on
    contour(p_value,[0.495 0.505],'linewidth',2,'color','r')
    set(gca,'fontsize',14)
end

handles.data_struct=data_struct;
guidata(hObject, handles)
msgbox('Finished loading modeled data')


% --- Executes on button press in transform_sessions.
function transform_sessions_Callback(hObject,~, handles)
% hObject    handle to transform_sessions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Stage 2: Transforming all the sessions to a reference coordinate system using
% rigid-body transformation

% This callback performs rigid-body transfomration to all the sessions
% according to a chosen reference ssseion. This stage includes:
% 1) Corrects sessions for translation/rotations and transforming the
% spatial footprints into a single coordinate frame
% 2) Matches the sizes of all the spatial footprints from the different sessions
% to the intersction of the different FOVs
% 3) Tests the data to evaluate whether or not it is suitable for
% longitudinal analysis

data_struct=handles.data_struct;
button_cent=data_struct.button_cent;
if strcmp(button_cent,'Yes');
    all_centroids=data_struct.all_centroids;
end

button_sub_pixel='Yes';

all_filters=data_struct.all_filters;
center=data_struct.center;
num_sessions=data_struct.num_sessions;
pixel_to_mic=data_struct.pixel_to_mic;

reference_session=str2num(get(handles.reference_session,'string'));
reference_valid=1;
if isempty(reference_session) || reference_session<1 || reference_session>num_sessions
    reference_valid=0;
end
while reference_valid==0
    set(handles.reference_session,'value',0)
    msgbox('Please insert a valid reference session number and press enter ');
    waitfor(handles.reference_session,'value',1);
    reference_session=str2num(get(handles.reference_session,'string'));
    if ~isempty(reference_session) && reference_session>=1 && reference_session<=num_sessions
        reference_valid=1;
    end
end
data_struct.reference_session=reference_session;

if num_sessions<2
    errordlg('Please choose more than one session to perform cell registration ')
else
    trans_val=get(handles.translations,'Value');
    if trans_val==1
        button_0='No';
    else
        button_0='Yes';
    end
    
    % Normalizing all footprints (to sum up to 1) and put threshold on pixels
    pixel_thresh=0.5;
    h=waitbar(0,'Normalizing spatial footprints');
    for n=1:num_sessions
        waitbar((n)/(num_sessions),h,['Normalizing spatial footprints for session number ' num2str(n) '/' num2str(num_sessions)])
        this_session_filters=all_filters{n};
        num_filters=size(this_session_filters,1);
        for k=1:num_filters
            temp_filter=squeeze(this_session_filters(k,:,:));
            max_filter=max(temp_filter(:));
            temp_filter(temp_filter<pixel_thresh*max_filter)=0;
            temp_filter=temp_filter./sum(sum(temp_filter));
            all_filters{n}(k,:,:)=temp_filter;
        end
    end
    close(h);
    
    % Calculating all centroids
    if strcmp(button_cent,'No');
        all_centroids=cell(1,num_sessions);
        properties ='all';
        h=waitbar(0,'Calculating centroids locations');
        for n=1:num_sessions
            waitbar((n)/(num_sessions),h,['Calculating centroids locations for session number ' num2str(n) '/' num2str(num_sessions)])
            this_session_filters=all_filters{n};
            num_filters=size(this_session_filters,1);
            all_centroids{n}=zeros(num_filters,2);
            for k=1:num_filters
                temp_filter=squeeze(this_session_filters(k,:,:));
                level = graythresh(temp_filter); % Using otso method
                bw = im2bw(temp_filter ,level);
                filtersFeatures = regionprops(bw,properties);
                numOfROIs = numel(filtersFeatures);
                if numOfROIs>1
                    areas =cell2mat( {filtersFeatures.Area});
                    [~,idxOfLargestArea] = max(areas);
                    filtersFeatures=filtersFeatures(idxOfLargestArea);
                end
                all_centroids{n}(k,:)=filtersFeatures.Centroid;
            end
        end
        close(h);
    end
    
    % Calculating all footprints projections:
    all_projections=cell(1,num_sessions);
    reference_session_size=[size(all_filters{reference_session},2) size(all_filters{reference_session},3)];
    for n=1:num_sessions
        this_session_filters=all_filters{n};
        if size(this_session_filters,2)>1.1*reference_session_size(1)
            warndlg(['The y axis FOV in session number ' num2str(n) ' is significantly larger than the reference'])
        elseif size(this_session_filters,2)<0.9*reference_session_size(1)
            warndlg(['The y axis FOV in session number ' num2str(n) ' is significantly smaller than the reference'])
        end
        if size(this_session_filters,3)>1.1*reference_session_size(2)
            warndlg(['The x axis FOV in session number ' num2str(n) ' is significantly larger than the reference'])
        elseif size(this_session_filters,3)<0.9*reference_session_size(2)
            warndlg(['The x axis FOV in session number ' num2str(n) ' is significantly smaller than the reference'])
        end
        
        num_filters=size(this_session_filters,1);
        normalized_filters=zeros(size(this_session_filters));
        for k=1:num_filters
            normalized_filters(k,:,:)=this_session_filters(k,:,:)/max(max(this_session_filters(k,:,:)));
        end
        all_projections{n}=squeeze(sum(normalized_filters,1));
    end
    
    subx=4;
    suby=ceil(num_sessions/subx);
    if num_sessions>4
        figure('units','normalized','outerposition',[0.1 0.04 0.8 0.96])
        for n=1:num_sessions
            subplot(suby,subx,n)
            imagesc(all_projections{n},[0 2])
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            colormap('gray')
            this_session_name=data_struct.sessions_list{n+1};
            title(strrep(this_session_name,'_','\_'),'fontsize',14,'fontweight','bold')
        end
    else
        figure('units','normalized','outerposition',[0.1 0.2 0.8 0.5])
        for n=1:num_sessions
            subplot(1,num_sessions,n)
            imagesc(all_projections{n},[0 2])
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            colormap('gray')
            this_session_name=data_struct.sessions_list{n+1};
            title(strrep(this_session_name,'_','\_'),'fontsize',14,'fontweight','bold')
        end
    end
    figures_dir=data_struct.figures_dir;
    cd(figures_dir)
    set(gcf,'PaperPositionMode','auto')
    savefig('Stage 1 - spatial footprints')
    saveas(gcf,'Stage 1 - spatial footprints','tif')
    
    % Calculating all centroids projections:
    centroid_projections=cell(1,num_sessions);
    for n=1:num_sessions
        this_session_centroids=all_centroids{n};
        this_session_filters=all_filters{n};
        num_filters=size(this_session_filters,1);
        normalized_centroids=zeros(size(this_session_filters));
        for k=1:num_filters
            if round(this_session_centroids(k,2))>1.5 && round(this_session_centroids(k,1))>1.5 && round(this_session_centroids(k,2))<size(normalized_centroids,2)-1 && round(this_session_centroids(k,1))<size(normalized_centroids,3)-1
                normalized_centroids(k,round(this_session_centroids(k,2))-1:round(this_session_centroids(k,2))+1,round(this_session_centroids(k,1))-1:round(this_session_centroids(k,1))+1)=1/4;
                normalized_centroids(k,round(this_session_centroids(k,2))-1:round(this_session_centroids(k,2))+1,round(this_session_centroids(k,1)))=1/2;
                normalized_centroids(k,round(this_session_centroids(k,2)),round(this_session_centroids(k,1))-1:round(this_session_centroids(k,1))+1)=1/2;
                normalized_centroids(k,round(this_session_centroids(k,2)),round(this_session_centroids(k,1)))=1;
            elseif round(this_session_centroids(k,2))>0 && round(this_session_centroids(k,1))>0 && round(this_session_centroids(k,2))<size(normalized_centroids,2) && round(this_session_centroids(k,1))<size(normalized_centroids,3)
                normalized_centroids(k,round(this_session_centroids(k,2)),round(this_session_centroids(k,1)))=1;
            end
        end
        centroid_projections{n}=squeeze(sum(normalized_centroids,1));
    end
    
    % adjusting the image sizes:
    max_x=0;
    max_y=0;
    
    overlapping_cell_all_sessions=cell(1,num_sessions);
    for n=1:num_sessions
        this_session_projection=all_projections{n};
        overlapping_cell_all_sessions{n}=ones(size(this_session_projection,1),size(this_session_projection,2));
        if n==1
            max_x=size(this_session_projection,2);
            max_y=size(this_session_projection,1);
        else
            max_x=max(max_x,size(this_session_projection,2));
            max_y=max(max_y,size(this_session_projection,1));
        end
    end
    
    overlapping_matrix_all_sessions=zeros(max_y,max_x,num_sessions);
    for n=1:num_sessions
        projection_temp=all_projections{n};
        overlapping_matrix_temp=overlapping_cell_all_sessions{n};
        centroid_projections_temp=centroid_projections{n};
        new_projections=zeros(max_y,max_x);
        new_overlapping_matrix=zeros(max_y,max_x);
        new_centroid_projections=zeros(max_y,max_x);
        new_projections(1:size(projection_temp,1),1:size(projection_temp,2))=projection_temp;
        new_overlapping_matrix(1:size(projection_temp,1),1:size(projection_temp,2))=overlapping_matrix_temp;
        new_centroid_projections(1:size(centroid_projections_temp,1),1:size(centroid_projections_temp,2))=centroid_projections_temp;
        all_projections{n}=new_projections;
        overlapping_matrix_all_sessions(:,:,n)=new_overlapping_matrix;
        centroid_projections{n}=new_centroid_projections;
        all_filters_temp=all_filters{n};
        num_filters=size(all_filters_temp,1);
        new_filters=zeros(num_filters,max_y,max_x);
        for k=1:num_filters
            filter_temp=squeeze(all_filters_temp(k,:,:));
            new_filter=zeros(max_y,max_x);
            new_filter(1:size(filter_temp,1),1:size(filter_temp,2))=filter_temp;
            new_filters(k,:,:)=new_filter;
        end
        all_filters{n}=new_filters;
    end
    
    % Finding the translations/rotations for across-session registration:
    all_centroids_corrected=all_centroids;
    all_projections_corrected=all_projections;
    all_filters_corrected=all_filters;
    centroid_projections_corrected=centroid_projections;
    
    x_ind_vec=zeros(1,num_sessions-1);
    y_ind_vec=zeros(1,num_sessions-1);
    
    if strcmp(button_0,'Yes') % correcting for rotations as well
        max_rotation=str2num(get(handles.maximal_rotation,'string'));
        rotation_vec=zeros(1,num_sessions-1);
        possible_rotations=-max_rotation:0.5:max_rotation;
        all_rotated_projections=cell(1,num_sessions);
        all_centroids_rotated=cell(1,num_sessions);
        all_rotated_projections{reference_session}=all_projections_corrected{reference_session};
        all_centroids_rotated{reference_session}=centroid_projections_corrected{reference_session};
    end
    
    max_correlation_vec=zeros(1,num_sessions-1);
    if strcmp(button_0,'Yes')
        best_rotation_vec=zeros(1,num_sessions-1);
    end
    best_x_translation_vec=zeros(1,num_sessions-1);
    best_y_translation_vec=zeros(1,num_sessions-1);
    h2=waitbar(0,'Transforming sessions','Units', 'normalized', 'Position', [0.4 0.5 0.2 0.07]);
    cross_corr_partial_vec=cell(1,num_sessions-1);
    registration_order=setdiff(1:num_sessions,reference_session);
    overlapping_matrix=ones(max_y,max_x);
    overlapping_matrix=overlapping_matrix.*overlapping_matrix_all_sessions(:,:,reference_session);
    for n=1:num_sessions-1
        overlapping_matrix_temp=overlapping_matrix_all_sessions(:,:,n);
        waitbar((n)/(num_sessions-1),h2,['Transforming session number ' num2str(n) '/' num2str(num_sessions-1)])
        if strcmp(button_0,'Yes') % correcting for rotations as well
            max_correlation=0;
            best_rotation=0;
            h1=waitbar(0,'Checking for rotations','Units', 'normalized', 'Position',[0.4 0.4 0.2 0.07]);
            for k=1:length(possible_rotations)
                waitbar((k)/length(possible_rotations),h1,'Checking for rotations')
                rotated_image=rotate_image_interp(centroid_projections_corrected{registration_order(n)},possible_rotations(k),[0 0],center);
                cross_corr=normxcorr2(centroid_projections_corrected{reference_session},rotated_image);
                if max(max(cross_corr))>max_correlation
                    best_rotation=possible_rotations(k);
                    max_correlation=max(max(cross_corr));
                end
            end
            close(h1);
            rotation_vec(n)=best_rotation;
            if abs(best_rotation)>0
                if abs(best_rotation)>10
                    warndlg([num2str(best_rotation) ' degrees rotation was found for session number ' num2str(registration_order(n))])
                end
                rotated_projections=rotate_image_interp(all_projections_corrected{registration_order(n)}',-best_rotation,[0 0],center);
                overlapping_matrix_temp=rotate_image_interp(overlapping_matrix_all_sessions(:,:,n)',-best_rotation,[0 0],center)';
                all_rotated_projections{registration_order(n)}=rotated_projections';
                theta=best_rotation*pi/180;
                transformation=[cos(theta) -sin(theta) ; sin(theta) cos(theta)]';
                trans_inv=transformation^-1;
                centroids_temp=trans_inv*(all_centroids_corrected{registration_order(n)}-repmat([center(1) ;center(2)],1,size(all_centroids_corrected{registration_order(n)},1))')'+repmat([center(1) ;center(2)],1,size(all_centroids_corrected{registration_order(n)},1));
            else
                all_rotated_projections{registration_order(n)}=all_projections_corrected{registration_order(n)};
                centroids_temp=all_centroids_corrected{registration_order(n)}';
            end
            all_centroids_corrected{registration_order(n)}=centroids_temp';
            this_session_centroids=all_centroids_corrected{registration_order(n)};
            num_filters=size(this_session_centroids,1);
            this_session_filters=all_filters_corrected{registration_order(n)};
            normalized_centroids=zeros(size(this_session_filters));
            for k=1:num_filters
                if round(this_session_centroids(k,2))>1.5 && round(this_session_centroids(k,1))>1.5 && round(this_session_centroids(k,2))<size(normalized_centroids,2)-1 && round(this_session_centroids(k,1))<size(normalized_centroids,3)-1
                    normalized_centroids(k,round(this_session_centroids(k,2))-1:round(this_session_centroids(k,2))+1,round(this_session_centroids(k,1))-1:round(this_session_centroids(k,1))+1)=1/4;
                    normalized_centroids(k,round(this_session_centroids(k,2))-1:round(this_session_centroids(k,2))+1,round(this_session_centroids(k,1)))=1/2;
                    normalized_centroids(k,round(this_session_centroids(k,2)),round(this_session_centroids(k,1))-1:round(this_session_centroids(k,1))+1)=1/2;
                    normalized_centroids(k,round(this_session_centroids(k,2)),round(this_session_centroids(k,1)))=1;
                elseif round(this_session_centroids(k,2))>0.5 && round(this_session_centroids(k,1))>0.5 && round(this_session_centroids(k,2))<size(normalized_centroids,2) && round(this_session_centroids(k,1))<size(normalized_centroids,3)
                    normalized_centroids(k,round(this_session_centroids(k,2)),round(this_session_centroids(k,1)))=1;
                end
            end
            all_centroids_rotated{registration_order(n)}=squeeze(sum(normalized_centroids,1));
            
            if abs(best_rotation)>0 && max_correlation>0.15;                
                h = waitbar(0,'Rotating spatial footprints','Units', 'normalized', 'Position',[0.4 0.4 0.2 0.07]);
                for m=1:size(all_centroids_corrected{registration_order(n)},1);
                    waitbar((m)/size(all_centroids_corrected{registration_order(n)},1),h,['Rotating spatial footprint number ' num2str(m) '/' num2str(size(all_centroids_corrected{registration_order(n)},1))])
                    temp_filter=squeeze(all_filters_corrected{registration_order(n)}(m,:,:));
                    temp_centroid=all_centroids{registration_order(n)}(m,:);
                    temp_filter_rotated=rotate_cell(temp_filter',-best_rotation,[0 0],center,temp_centroid,pixel_to_mic);
                    all_filters_corrected{registration_order(n)}(m,:,:)=temp_filter_rotated';
                end
                close(h);
            end
        end
        if strcmp(button_0,'Yes')
            cross_corr_cent=normxcorr2(all_centroids_rotated{reference_session},all_centroids_rotated{registration_order(n)});
        else
            cross_corr_cent=normxcorr2(centroid_projections_corrected{reference_session},centroid_projections_corrected{registration_order(n)});
        end
        cross_corr_size=size(cross_corr_cent);
        cross_corr_partial=cross_corr_cent(round(cross_corr_size(1)/2-cross_corr_size(1)/6):round(cross_corr_size(1)/2+cross_corr_size(1)/6)...
            ,round(cross_corr_size(2)/2-cross_corr_size(2)/6):round(cross_corr_size(2)/2+cross_corr_size(2)/6));
        [max_correlation_vec(n),x_ind]=max(max(cross_corr_partial));
        if strcmp(button_0,'Yes')
            best_rotation_vec(n)=best_rotation;
        end
        [~,y_ind]=max(cross_corr_partial(:,x_ind));
        if strcmp(button_sub_pixel,'Yes')
            temp_corr_x=cross_corr_partial(y_ind-1:y_ind+1,x_ind-5:x_ind+5);
            [~,mu_x_1]=gaussfit(-5:5,temp_corr_x(1,:)./sum(temp_corr_x(1,:)),1,0);
            [~,mu_x_2]=gaussfit(-5:5,temp_corr_x(2,:)./sum(temp_corr_x(2,:)),1,0);
            [~,mu_x_3]=gaussfit(-5:5,temp_corr_x(3,:)./sum(temp_corr_x(2,:)),1,0);
            sub_x=mean([mu_x_1 , mu_x_2 , mu_x_2 , mu_x_3]);
            if abs(sub_x)>1
                warndlg(['X axis sub-pixel correction was ' num2str(round(100*sub_x)/100) ' for session number ' num2str(registration_order(n))])
            end
            temp_corr_y=cross_corr_partial(y_ind-5:y_ind+5,x_ind-1:x_ind+1);
            [~,mu_y_1]=gaussfit(-5:5,temp_corr_y(:,1)./sum(temp_corr_y(:,1)),1,0);
            [~,mu_y_2]=gaussfit(-5:5,temp_corr_y(:,2)./sum(temp_corr_y(:,2)),1,0);
            [~,mu_y_3]=gaussfit(-5:5,temp_corr_y(:,3)./sum(temp_corr_y(:,3)),1,0);
            sub_y=mean([mu_y_1 , mu_y_2 , mu_y_2 , mu_y_3]);
            if abs(sub_y)>1
                warndlg(['Y axis sub-pixel correction was ' num2str(round(100*sub_y)/100) ' for session number ' num2str(registration_order(n))])
            end
            x_ind=x_ind+round(cross_corr_size(2)/2-cross_corr_size(2)/6)-1;
            y_ind=y_ind+round(cross_corr_size(1)/2-cross_corr_size(1)/6)-1;
            x_ind_sub=x_ind+sub_x;
            y_ind_sub=y_ind+sub_y;
            x_ind_vec(n)=x_ind_sub;
            y_ind_vec(n)=y_ind_sub;
            best_x_translation_vec(n)=(x_ind_sub-max_x);
            best_y_translation_vec(n)=(y_ind_sub-max_y);
        else
            x_ind=x_ind+round(cross_corr_size(2)/2-cross_corr_size(2)/6)-1;
            y_ind=y_ind+round(cross_corr_size(1)/2-cross_corr_size(1)/6)-1;
            x_ind_vec(n)=x_ind;
            y_ind_vec(n)=y_ind;
            best_x_translation_vec(n)=(x_ind-max_x);
            best_y_translation_vec(n)=(y_ind-max_y);
        end
        
        if pixel_to_mic*abs(x_ind-max_x)>50
            warndlg([num2str(pixel_to_mic*abs(x_ind-max_x)) ' micron x translation was found for session number ' num2str(registration_order(n))])
        end
        if pixel_to_mic*abs(y_ind-max_y)>50
            warndlg([num2str(pixel_to_mic*abs(y_ind-max_y)) ' micron y translation was found for session number ' num2str(registration_order(n))])
        end
        
        if max_correlation_vec(n)<0.15
            if strcmp(button_0,'Yes')
                warndlg(['No appropriate translations/rotations were found for session number ' num2str(registration_order(n))])
            else
                warndlg(['No appropriate translations were found for session number ' num2str(registration_order(n)) ' - consider using rotations as well'])
            end
        end
        if strcmp(button_0,'Yes')
            projections_temp=all_rotated_projections{registration_order(n)};
            centroid_projections_temp=all_centroids_rotated{registration_order(n)};
            centroids_temp=all_centroids_corrected{registration_order(n)};
        else
            projections_temp=all_projections{registration_order(n)};
            centroid_projections_temp=centroid_projections_corrected{registration_order(n)};
            centroids_temp=all_centroids{registration_order(n)};
        end
        
        if strcmp(button_sub_pixel,'Yes')
            centroids_temp(:,1)=centroids_temp(:,1)-(x_ind_sub-max_x);
            centroids_temp(:,2)=centroids_temp(:,2)-(y_ind_sub-max_y);
        else
            centroids_temp(:,1)=centroids_temp(:,1)-(x_ind-max_x);
            centroids_temp(:,2)=centroids_temp(:,2)-(y_ind-max_y);
        end
        all_centroids_corrected{registration_order(n)}=centroids_temp;
        
        new_projections=zeros(max_y,max_x);
        new_centroid_projections=zeros(max_y,max_x);
        new_overlapping_matrix_temp=zeros(max_y,max_x);
        if y_ind-max_y>=0 && x_ind-max_x>=0
            new_projections(1:end-(y_ind-max_y),1:end-(x_ind-max_x))=projections_temp(1+(y_ind-max_y):end,1+(x_ind-max_x):end);
            new_centroid_projections(1:end-(y_ind-max_y),1:end-(x_ind-max_x))=centroid_projections_temp(1+(y_ind-max_y):end,1+(x_ind-max_x):end);
            new_overlapping_matrix_temp(1:end-(y_ind-max_y),1:end-(x_ind-max_x))=overlapping_matrix_temp(1+(y_ind-max_y):end,1+(x_ind-max_x):end);
        elseif y_ind-max_y>=0 && x_ind-max_x<0
            new_projections(1:end-(y_ind-max_y),1-(x_ind-max_x):end)=projections_temp(1+(y_ind-max_y):end,1:end+(x_ind-max_x));
            new_centroid_projections(1:end-(y_ind-max_y),1-(x_ind-max_x):end)=centroid_projections_temp(1+(y_ind-max_y):end,1:end+(x_ind-max_x));
            new_overlapping_matrix_temp(1:end-(y_ind-max_y),1-(x_ind-max_x):end)=overlapping_matrix_temp(1+(y_ind-max_y):end,1:end+(x_ind-max_x));
        elseif y_ind-max_y<0 && x_ind-max_x>=0
            new_projections(1-(y_ind-max_y):end,1:end-(x_ind-max_x))=projections_temp(1:end+(y_ind-max_y),1+(x_ind-max_x):end);
            new_centroid_projections(1-(y_ind-max_y):end,1:end-(x_ind-max_x))=centroid_projections_temp(1:end+(y_ind-max_y),1+(x_ind-max_x):end);
            new_overlapping_matrix_temp(1-(y_ind-max_y):end,1:end-(x_ind-max_x))=overlapping_matrix_temp(1:end+(y_ind-max_y),1+(x_ind-max_x):end);
        elseif y_ind-max_y<0 && x_ind-max_x<0
            new_projections(1-(y_ind-max_y):end,1-(x_ind-max_x):end)=projections_temp(1:end+(y_ind-max_y),1:end+(x_ind-max_x));
            new_centroid_projections(1-(y_ind-max_y):end,1-(x_ind-max_x):end)=centroid_projections_temp(1:end+(y_ind-max_y),1:end+(x_ind-max_x));
            new_overlapping_matrix_temp(1-(y_ind-max_y):end,1-(x_ind-max_x):end)=overlapping_matrix_temp(1:end+(y_ind-max_y),1:end+(x_ind-max_x));
        end
        if strcmp(button_sub_pixel,'Yes')
            temp_projections=translate_projections(new_projections',[sub_y sub_x]);
            temp_centroid_projections=translate_projections(new_centroid_projections',[sub_y sub_x]);
            all_projections_corrected{registration_order(n)}=temp_projections';
            centroid_projections_corrected{registration_order(n)}=temp_centroid_projections';
        else
            all_projections_corrected{registration_order(n)}=new_projections;
            centroid_projections_corrected{registration_order(n)}=new_centroid_projections;
        end
        
        if strcmp(button_0,'Yes')
            all_filters_temp=all_filters_corrected{registration_order(n)};
        else
            all_filters_temp=all_filters{registration_order(n)};
        end
        num_filters=size(all_filters_temp,1);
        h = waitbar(0,'Translating spatial footprints','Units', 'normalized', 'Position',[0.4 0.4 0.2 0.07]);
        for k=1:num_filters
            waitbar((k)/size(all_centroids_corrected{registration_order(n)},1),h,['Translating spatial footprint number ' num2str(k) '/' num2str(size(all_centroids_corrected{registration_order(n)},1))])
            filters_temp=squeeze(all_filters_temp(k,:,:));
            new_filter=zeros(max_y,max_x);
            if y_ind-max_y>=0 && x_ind-max_x>=0
                new_filter(1:end-(y_ind-max_y),1:end-(x_ind-max_x))=filters_temp(1+(y_ind-max_y):end,1+(x_ind-max_x):end);
            elseif y_ind-max_y>=0 && x_ind-max_x<0
                new_filter(1:end-(y_ind-max_y),1-(x_ind-max_x):end)=filters_temp(1+(y_ind-max_y):end,1:end+(x_ind-max_x));
            elseif y_ind-max_y<0 && x_ind-max_x>=0
                new_filter(1-(y_ind-max_y):end,1:end-(x_ind-max_x))=filters_temp(1:end+(y_ind-max_y),1+(x_ind-max_x):end);
            elseif y_ind-max_y<0 && x_ind-max_x<0
                new_filter(1-(y_ind-max_y):end,1-(x_ind-max_x):end)=filters_temp(1:end+(y_ind-max_y),1:end+(x_ind-max_x));
            end
            if strcmp(button_sub_pixel,'Yes')
                temp_centroid=all_centroids_corrected{registration_order(n)}(k,:);
                temp_filter=translate_cell(new_filter',[sub_y sub_x],temp_centroid,pixel_to_mic);
                all_filters_temp(k,:,:)=temp_filter';
            else
                all_filters_temp(k,:,:)=new_filter;
            end
        end
        close(h);
        all_filters_corrected{registration_order(n)}=all_filters_temp;
        cross_corr_partial_vec{n}=cross_corr_partial;
        overlapping_matrix=overlapping_matrix.*new_overlapping_matrix_temp;
    end
    close(h2);
    
    % Visualizing 3 sessions by RGB coloring before and after registraion:
    rgb_ind=[1 2 3];
    all_projections_rgb=zeros(max_y,max_x,3);
    all_projections_rgb(:,:,1)=all_projections{rgb_ind(1)};
    all_projections_rgb(:,:,2)=all_projections{rgb_ind(2)};
    if num_sessions>2
        all_projections_rgb(:,:,3)=all_projections{rgb_ind(3)};
    end
    all_projections_rgb(all_projections_rgb>1)=1;
    
    all_projections_corrected_rgb=zeros(max_y,max_x,3);
    all_projections_corrected_rgb(:,:,1)=all_projections_corrected{rgb_ind(1)};
    all_projections_corrected_rgb(:,:,2)=all_projections_corrected{rgb_ind(2)};
    if num_sessions>2
        all_projections_corrected_rgb(:,:,3)=all_projections_corrected{rgb_ind(3)};
    end
    all_projections_corrected_rgb=all_projections_corrected_rgb+0.25*repmat(overlapping_matrix,1,1,3);
    all_projections_corrected_rgb(all_projections_corrected_rgb>1)=1;            
    
    figure('units','normalized','outerposition',[0 0.04 1 0.96])
    subplot(2,2,1)
    imshow(all_projections_rgb)
    if num_sessions>2
        title('RGB overlay (sessions 1-3): Pre-transformation','FontWeight','Bold','fontsize',18)
    else
        title('RGB overlay: Pre-transformation','FontWeight','Bold','fontsize',18)
    end
    subplot(2,2,2)
    imshow(all_projections_corrected_rgb)
    if num_sessions>3
        title('RGB overlay (sessions 1-3): Post-transformation','FontWeight','Bold','fontsize',18)
    else
        title('RGB overlay: Post-transformation','FontWeight','Bold','fontsize',18)
    end
    legend_strings=cell(1,num_sessions);
    if num_sessions>16
        color_vec=rand(num_sessions,3);
    else
        color_vec=[1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 0; 0 0 0; 0.4 0.7 0.4; 1 0.5 0; 0.8 0.8 0.8;...
            0.2 0.1 0.5; 0.4 0.3 0.2; 0.7 0 0.3; 0.6 0.9 0.2 ; 0.1 0.2 0.3; 0.3 0.8 0.3];
    end
    
    subplot(2,2,3)
    for n=1:num_sessions
        centroids=all_centroids{n};
        h=scatter(centroids(:,1),max_y-centroids(:,2),10);
        set(h,'MarkerFaceColor',color_vec(n,:));
        hold on
        legend_strings{n}=['S. ' num2str(n)];
    end
    title('Centroid locations: Pre-transformation','FontWeight','Bold','FontSize',22)
    xlim([0 size(all_filters{1},3)])
    ylim([0 size(all_filters{1},2)])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    legend(legend_strings)
    legend('boxoff')
    
    subplot(2,2,4)
    for n=1:num_sessions
        centroids=all_centroids_corrected{n};
        h=scatter(centroids(:,1),max_y-centroids(:,2),15);
        set(h,'MarkerFaceColor',color_vec(n,:));
        hold on
    end
    title('Centroid locations: Post-transformation','FontWeight','Bold','FontSize',22)
    xlim([0 size(all_filters_corrected{1},3)])
    ylim([0 size(all_filters_corrected{1},2)])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    legend(legend_strings)
    legend('boxoff')
    set(gcf,'PaperPositionMode','auto')
    savefig('Stage 2 - pre vs post transformation')
    saveas(gcf,'Stage 2 - pre vs post transformation','tif')
    
    axes(handles.axes1);
    imagesc(all_projections_corrected_rgb)
    title('RGB overlay','FontWeight','Bold','fontsize',16) 
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    
    data_struct.all_centroids_corrected=all_centroids_corrected;
    data_struct.all_filters_corrected=all_filters_corrected;
    data_struct.max_x=max_x;
    data_struct.max_y=max_y;
    data_struct.all_projections=all_projections;
    data_struct.all_projections_corrected=all_projections_corrected;
    data_struct.overlapping_matrix=overlapping_matrix;
            
    results_dir=data_struct.results_dir;
    cd(results_dir);
    handles.data_struct=data_struct;
    save('transformed_data_struct.mat','data_struct','-v7.3')
    
    % Evaluating the data if it is suitable for longitudinal analysis:
    if num_sessions>2
        all_projections_correlations=zeros(num_sessions,num_sessions);
        for n=1:num_sessions
            all_projections_correlations(n,n)=1;
            for k=n+1:num_sessions
                all_projections_correlations(n,k)=corr2(centroid_projections_corrected{n},centroid_projections_corrected{k});
                all_projections_correlations(k,n)=all_projections_correlations(n,k);
            end
        end
        
        figure('units','normalized','outerposition',[0.2 0.3 0.6 0.4])
        subplot(1,2,1)
        imagesc(all_projections_correlations)
        colormap('jet')
        colorbar
        caxis([0 1])
        x_label=cell(1,num_sessions);
        if num_sessions<8
            for n=1:num_sessions
                x_label{n}=num2str(n);
            end
        else
            for n=1:3:num_sessions
                x_label{n}=num2str(n);
            end
        end
        x=1:num_sessions;
        y_label=cell(1,num_sessions);
        if num_sessions<8
            for n=1:num_sessions
                y_label{n}=num2str(n);
            end
        else
            for n=1:3:num_sessions
                y_label{n}=num2str(n);
            end
        end
        y=1:num_sessions;
        hold on
        plot([0.5 0.5 num_sessions+0.5 num_sessions+0.5 0.5],[reference_session-0.5 reference_session+0.5 reference_session+0.5 reference_session-0.5 reference_session-0.5],'linewidth',3,'color','k')
        set(gca,'fontsize',16)
        set(gca,'XTick',x)
        set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
        set(gca,'YTick',y)
        set(gca,'YTickLabel',y_label,'fontsize',12,'fontweight','bold')
        title('Maximal Correlation across sessions','fontweight','bold')
        set(gca,'fontsize',14)
        subplot(1,2,2)
        for n=1:num_sessions
            ind_to_plot=setdiff(1:num_sessions,n);
            plot(1:num_sessions-1,all_projections_correlations(n,ind_to_plot),'*--','linewidth',2,'markersize',8,'color',color_vec(n,:));
            hold on
        end
        ylim([0 1])
        xlim([0 num_sessions])
        legend_strings=cell(1,num_sessions);
        for n=1:num_sessions
            legend_strings{n}=['Ref. S. - ' num2str(n)];
        end
        legend(legend_strings)
        legend('boxoff')
        if num_sessions<8
            for n=1:num_sessions
                x_label{n}=num2str(n);
            end
        else
            for n=1:3:num_sessions
                x_label{n}=num2str(n);
            end
        end
        x=1:num_sessions-1;
        set(gca,'XTick',x)
        set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
        set(gca,'fontsize',14)
        xlabel('Session number','FontWeight','Bold','fontsize',16)
        ylabel('Maximal Correlation','FontWeight','Bold','fontsize',16)
        set(gcf,'PaperPositionMode','auto')
        cd(figures_dir);
        savefig('Stage 2 - abnormalities test - Correlations')
        saveas(gcf,'Stage 2 - abnormalities test - Correlations','tif')
        max_correlation_vec=all_projections_correlations(reference_session,registration_order);
    end   
    cd(results_dir);
    data_struct.max_correlation_vec=max_correlation_vec;

    figure('units','normalized','outerposition',[0 0.04 1 0.96])
    if strcmp(button_0,'Yes')
        subplot(2,3,2)
        plot(registration_order,best_rotation_vec,'*','linewidth',2,'markersize',8,'color','b')
        xlim([0 num_sessions+1])
        ylim([-30 30])
        hold on
        plot([0 num_sessions+1],[-10 -10],'--','color','k','linewidth',2)
        hold on
        plot([0 num_sessions+1],[-20 -20],'--','color','k','linewidth',2)
        hold on
        plot([0 num_sessions+1],[10 10],'--','color','k','linewidth',2)
        hold on
        plot([0 num_sessions+1],[20 20],'--','color','k','linewidth',2)
        if num_sessions<8
            x_label=cell(1,num_sessions);
            for n=1:num_sessions
                x_label{n}=num2str(n);
            end
        else
            x_label=cell(1,length(1:3:num_sessions));
            for n=1:3:num_sessions
                x_label{n}=num2str(n);
            end
        end
        x_label{reference_session}='Ref.';
        x=1:num_sessions;
        set(gca,'XTick',x)
        set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
        y_label=-30:10:30;
        y=-30:10:30;
        set(gca,'YTick',y)
        set(gca,'YTickLabel',y_label,'fontsize',12,'fontweight','bold')
        xlabel('Session number','FontWeight','Bold','fontsize',16)
        ylabel('Rotation (degrees)','FontWeight','Bold','fontsize',16)
        set(gca,'fontsize',16)
    end
    
    subplot(2,3,1)
    for n=1:num_sessions-1
        plot(pixel_to_mic*best_x_translation_vec(n),pixel_to_mic*best_y_translation_vec(n),'*','markersize',8,'linewidth',2,'color',color_vec(n,:));
        hold on
    end
    legend(legend_strings{registration_order},'Location','NorthWest')
    legend('boxoff')
    x_for_circle=0:0.01:2*pi;
    plot(25*sin(x_for_circle),25*cos(x_for_circle),'--','color','k','linewidth',1);
    hold on
    plot(50*sin(x_for_circle),50*cos(x_for_circle),'--','color','k','linewidth',1);
    hold on
    plot(75*sin(x_for_circle),75*cos(x_for_circle),'--','color','k','linewidth',1);
    hold on
    plot(100*sin(x_for_circle),100*cos(x_for_circle),'--','color','k','linewidth',1);
    ylim([-100 100])
    xlim([-100 100])
    x_label=-100:50:100;
    x=-100:50:100;
    y=-100:50:100;
    y_label=-100:50:100;
    set(gca,'YTick',y)
    set(gca,'YTickLabel',y_label,'fontsize',14,'fontweight','bold')
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
    xlabel('Y translation (\mum)','FontWeight','Bold','fontsize',16)
    ylabel('X translation (\mum)','FontWeight','Bold','fontsize',16)
    
    num_cells_per_session=zeros(1,num_sessions);
    for n=1:num_sessions
        num_cells_per_session(n)=size(all_filters{n},1);
    end
    mean_num_cells=mean(num_cells_per_session);
    is_high_num_cells=zeros(1,num_sessions);
    is_low_num_cells=zeros(1,num_sessions);
    for n=1:num_sessions
        if num_cells_per_session(n)>1.5*mean_num_cells
            is_high_num_cells(n)=1;
        end
        if num_cells_per_session(n)<2/3*mean_num_cells
            is_low_num_cells(n)=1;
        end
    end
    if isfield(data_struct,'all_events')
        all_events=data_struct.all_events;
        frame_rate=data_struct.frame_rate;
        events_rate_per_session=zeros(1,num_sessions);
        std_events_rate_per_session=zeros(1,num_sessions);
        amp_per_session=zeros(1,num_sessions);
        std_amp_per_session=zeros(1,num_sessions);
        for n=1:num_sessions
            events_rate_per_session(n)=frame_rate*sum(sum(all_events{n}>0))/size(all_events{n},1)/size(all_events{n},2);
            std_events_rate_per_session(n)=std(frame_rate*sum(all_events{n}>0)./size(all_events{n},1));
            all_amps_temp=sum(all_events{n})./sum(all_events{n}>0);
            amp_per_session(n)=sum(sum(all_events{n}))/sum(sum(all_events{n}>0));
            std_amp_per_session(n)=std(all_amps_temp(~isnan(all_amps_temp)));
        end
        average_event_amp=mean(amp_per_session);
        average_events_rate=mean(events_rate_per_session); % should be changed to fs*
        is_high_event_amps=zeros(1,num_sessions);
        is_low_event_amps=zeros(1,num_sessions);
        is_high_event_rates=zeros(1,num_sessions);
        is_low_event_rates=zeros(1,num_sessions);
        for n=1:num_sessions
            if amp_per_session(n)>1.5*average_event_amp
                is_high_event_amps(n)=1;
                warndlg(['Extermely high event amplitudes are observed in session number ' num2str(n)])
            end
            if amp_per_session(n)<2/3*average_event_amp
                is_low_event_amps(n)=1;
                warndlg(['Extermely low event amplitudes are observed in session number ' num2str(n)])
            end
            if events_rate_per_session(n)>1.2*average_events_rate
                is_high_event_rates(n)=1;
            end
            if events_rate_per_session(n)<5/6*average_events_rate
                is_low_event_rates(n)=1;
            end
        end
        subplot(2,3,5)
        shadedErrorBar(1:num_sessions,events_rate_per_session,std_events_rate_per_session,'--*b')
        hold on
        plot([1 num_sessions],[average_events_rate average_events_rate],'--','linewidth',2,'color','r')
        xlim([0 num_sessions+1])
        ylim([0 max([2*average_events_rate,max(events_rate_per_session)])])
        xlabel('Session number','fontsize',16,'fontweight','bold')
        ylabel('Event rate (1/sec)','fontsize',16,'fontweight','bold')
        if num_sessions<8
            x_label=1:num_sessions;
        else
            x_label=cell(1,num_sessions);
            for n=1:3:num_sessions
                x_label{n}=num2str(n);
            end
        end
        x=1:num_sessions;
        set(gca,'XTick',x)
        set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
        set(gca,'fontsize',16)
        subplot(2,3,6)
        shadedErrorBar(1:num_sessions,amp_per_session,std_amp_per_session,'--*b')
        hold on
        plot([1 num_sessions],[average_event_amp average_event_amp],'--','linewidth',2,'color','r')
        xlim([0 num_sessions+1])
        ylim([0 max([2*average_event_amp,max(amp_per_session)])])
        xlabel('Session number','fontsize',16,'fontweight','bold')
        ylabel('Event amplitude (dF/F_0)','fontsize',16,'fontweight','bold')
        if num_sessions<8
            x_label=1:num_sessions;
        else
            x_label=cell(1,num_sessions);
            for n=1:3:num_sessions
                x_label{n}=num2str(n);
            end
        end
        x=1:num_sessions;
        set(gca,'XTick',x)
        set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
        set(gca,'fontsize',16)
    end
    subplot(2,3,4)
    plot(1:num_sessions,num_cells_per_session,'*--','markersize',8,'linewidth',2,'color','b')
    hold on
    plot([1 num_sessions],[mean_num_cells mean_num_cells],'--','linewidth',2,'color','r')
    xlim([0 num_sessions+1])
    ylim([0 max([2*mean_num_cells,max(num_cells_per_session)])])
    xlabel('Session number','fontsize',16,'fontweight','bold')
    ylabel('Number of cells','fontsize',16,'fontweight','bold')
    if num_sessions<8
        x_label=1:num_sessions;
    else
        x_label=cell(1,num_sessions);
        for n=1:3:num_sessions
            x_label{n}=num2str(n);
        end
    end
    x=1:num_sessions;
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
    set(gca,'fontsize',16)
    subplot(2,3,3)
    plot(registration_order,max_correlation_vec,'*--','markersize',8,'linewidth',2,'color','b')
    hold on
    plot([1 num_sessions],[mean(max_correlation_vec) mean(max_correlation_vec)],'--','linewidth',2,'color','r')
    xlim([0 num_sessions+1])
    ylim([0 max([2*mean(max_correlation_vec),max(max_correlation_vec)])])
    x_label=cell(1,num_sessions);
    if num_sessions<8
        for n=1:num_sessions
            x_label{n}=num2str(n);
        end
    else        
        for n=1:3:num_sessions
            x_label{n}=num2str(n);
        end
    end
    x_label{reference_session}='Ref.';
    x=1:num_sessions;
    set(gca,'fontsize',16)
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
    xlabel('Session number','FontWeight','Bold','fontsize',16)
    ylabel('Max. projection correlation','FontWeight','Bold','fontsize',16)
    set(gca,'fontsize',16)
    set(gcf,'PaperPositionMode','auto')
    cd(figures_dir);
    savefig('Stage 2 - abnormalities test - general')
    saveas(gcf,'Stage 2 - abnormalities test - general','tif')    
    
    for n=1:num_sessions 
        if is_high_num_cells(n)==1;
            warndlg(['Extermely high number of cells is observed in session number ' num2str(n)])
        end
        if is_low_num_cells(n)==1;
            warndlg(['Extermely low number of cells is observed in session number ' num2str(n)])
        end
        if isfield(data_struct,'all_events')
            if is_high_event_rates(n)==1;
                warndlg(['Extermely high event rates are observed in session number ' num2str(n)])
            end
            if is_low_event_rates(n)==1;
                warndlg(['Extermely low event rates are observed in session number ' num2str(n)])
            end
            if is_high_event_amps(n)==1;
                warndlg(['Extermely high event amplitudes are observed in session number ' num2str(n)])
            end
            if is_low_event_amps(n)==1;
                warndlg(['Extermely low event amplitudes are observed in session number ' num2str(n)])
            end
        end
    end
    
    cd(results_dir);
    handles.data_struct=data_struct;
    guidata(hObject, handles)
    msgbox('Finished transforming sessions')
end

% --- Executes on button press in compute_model.
function compute_model_Callback(hObject,~, handles)
% hObject    handle to compute_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Stage 3: Computing a probabilistic model of the spatial footprints similarities
% of neighboring cell-pairs from different sessions using the centroids
% distances and spatial correlations

% This callback computes the probability model for the same cells and
% different cells according to either spatial correlations, centroid
% distances, or both measures. The output is all the probabilities of
% neighboring cell-pairs to be the same cell - P(same).

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

cell_memory_allocation=10000;
data_struct=handles.data_struct;
all_filters_corrected=data_struct.all_filters_corrected;
all_centroids_corrected=data_struct.all_centroids_corrected;
num_sessions=data_struct.num_sessions;
pixel_to_mic=str2double(get(handles.pixel_to_mic,'string'));

if get(handles.use_joint_model,'value')==1
    errordlg('The joint model is not implemented in this version. Please uncheck this checkbox')
    error_variable=non_existing_variable;
end

if get(handles.use_joint_model,'value')==1 & get(handles.two_photon,'value')==1
    errordlg('The joint model is only applied for 1-photon imaging data');
    error_variable=non_existing_variable;
end

if num_sessions>7
    num_bins=100;
elseif num_sessions>5
    num_bins=80;
    elseif num_sessions>3
        num_bins=60;
else
    num_bins=50;
end
data_struct.number_of_bins=num_bins;

temp_dist_thresh=str2double(get(handles.model_maximal_distance,'string'))/pixel_to_mic;

compute_model_again=0;
if ~isfield(data_struct,'all_to_all_p_value_dist_multi')
    compute_model_again=1;
else
    max_measured_distance=max(data_struct.all_neighbor_distances);
    if max_measured_distance<0.99*temp_dist_thresh
        compute_model_again=1;
    end
end

cd(data_struct.results_dir);
mkdir('Figures');
figures_dir=fullfile(data_struct.results_dir,'Figures');
data_struct.figures_dir=figures_dir;

if compute_model_again % If the distributions were not estiamted yet
    % Computing correlations and distances across days:
    neighbor_count=0;
    NN_count=0;
    NNN_count=0;
    all_to_all_matrix=cell(1,num_sessions);
    all_to_all_correlation=cell(1,num_sessions);
    all_to_all_matrix_multi=cell(1,num_sessions);
    all_to_all_correlation_multi=cell(1,num_sessions);
    all_to_all_distance_multi=cell(1,num_sessions);
    all_neighbor_correlations=zeros(1,num_sessions^2*cell_memory_allocation);
    all_neighbor_distances=zeros(1,num_sessions^2*cell_memory_allocation);
    all_neighbor_distances_x=zeros(1,num_sessions^2*cell_memory_allocation);
    all_neighbor_distances_y=zeros(1,num_sessions^2*cell_memory_allocation);
    NN_correlations=zeros(1,num_sessions^2*cell_memory_allocation);
    NNN_correlations=zeros(1,num_sessions^2*cell_memory_allocation);
    NN_distances=zeros(1,num_sessions^2*cell_memory_allocation);
    NNN_distances=zeros(1,num_sessions^2*cell_memory_allocation);
    h = waitbar(0,'Computing correlations and distances','Units', 'normalized', 'Position',[0.4 0.5 0.2 0.07]);
    for n=1:num_sessions
        waitbar((n)/num_sessions,h,['Computing correlations and distances - session number ' num2str(n) '/' num2str(num_sessions)])
        new_filters=all_filters_corrected{n};
        new_centroids=all_centroids_corrected{n};
        num_cells=size(new_filters,1);
        all_to_all_matrix{n}=zeros(num_cells,num_sessions);
        all_to_all_correlation{n}=zeros(num_cells,num_sessions);
        all_to_all_correlation_multi{n}=cell(num_cells,num_sessions);
        all_to_all_distance_multi{n}=cell(num_cells,num_sessions);
        all_to_all_matrix_multi{n}=cell(num_cells,num_sessions);
        sessions_to_compare=1:num_sessions;
        sessions_to_compare(n)=[];
        h2 = waitbar(0,'Computing correlations and distances','Units', 'normalized', 'Position',[0.4 0.4 0.2 0.07]);
        for k=1:num_cells
            waitbar((k)/num_cells,h2,['Computing correlations and distances - cell number ' num2str(k) '/' num2str(num_cells)]);
            new_filter=squeeze(new_filters(k,:,:));
            for m=1:length(sessions_to_compare)
                this_session=sessions_to_compare(m);
                this_session_centroids=all_centroids_corrected{this_session};
                centroid=repmat(new_centroids(k,:),size(this_session_centroids,1),1);
                this_session_filters=all_filters_corrected{this_session};
                distance_vec=sqrt(sum((centroid-this_session_centroids).^2,2));
                diff_temp=centroid-this_session_centroids;
                filters_to_check=find(distance_vec<temp_dist_thresh);
                distance_vec_x=diff_temp(filters_to_check,1);
                distance_vec_y=diff_temp(filters_to_check,2);
                this_distance_vec=distance_vec(filters_to_check);
                if ~isempty(filters_to_check)
                    corr_vec=zeros(1,length(filters_to_check));
                    num_empty_filters=0;
                    for l=1:length(filters_to_check)
                        suspected_filter=squeeze(this_session_filters(filters_to_check(l),:,:));
                        if sum(sum(suspected_filter))==0 || sum(sum(new_filter))==0
                            num_empty_filters=num_empty_filters+1;
                        else
                            neighbor_count=neighbor_count+1;
                            corr_vec(l)=corr2(suspected_filter,new_filter);
                            all_neighbor_correlations(neighbor_count)=corr_vec(l);
                            all_neighbor_distances(neighbor_count)=distance_vec(filters_to_check(l));
                            all_neighbor_distances_x(neighbor_count)=distance_vec_x(l);
                            all_neighbor_distances_y(neighbor_count)=distance_vec_y(l);
                        end
                    end                    
                    if num_empty_filters<length(filters_to_check);
                        NN_count=NN_count+1;
                        NN_correlations(NN_count)=max(corr_vec);
                        NN_distances(NN_count)=min(distance_vec(filters_to_check));
                        if length(filters_to_check)>1 && num_empty_filters<length(filters_to_check)-1
                            NNN_count=NNN_count+length(filters_to_check)-1;
                            [dist_vec_sorted,ind]=sort(distance_vec(filters_to_check));
                            dist_vec_sorted(1)=[];
                            NNN_distances(2+NNN_count-length(filters_to_check):NNN_count)=dist_vec_sorted(1:end);
                            [corr_vec_sorted]=corr_vec(ind);
                            corr_vec_sorted(1)=[];
                            NNN_correlations(2+NNN_count-length(filters_to_check):NNN_count)=corr_vec_sorted(1:end);
                        end
                    end
                    all_to_all_correlation_multi{n}{k,this_session}=corr_vec;
                    [highest_corr,highest_corr_ind]=max(corr_vec);
                    index=filters_to_check(highest_corr_ind);
                    all_to_all_matrix{n}(k,this_session)=index;
                    all_to_all_correlation{n}(k,this_session)=highest_corr;
                    all_to_all_distance_multi{n}{k,this_session}=this_distance_vec;
                    all_to_all_matrix_multi{n}{k,this_session}=filters_to_check;
                end
            end
        end
        close(h2);
    end
    close(h);
    all_neighbor_correlations(neighbor_count+1:end)=[];
    all_neighbor_distances(neighbor_count+1:end)=[];
    all_neighbor_distances_x(neighbor_count+1:end)=[];
    all_neighbor_distances_y(neighbor_count+1:end)=[];
    NN_correlations(NN_count+1:end)=[];
    NNN_correlations(NNN_count+1:end)=[];
    NN_distances(NN_count+1:end)=[];
    NNN_distances(NNN_count+1:end)=[];
    
    % saving all the estimated disributions:
    data_struct.all_to_all_matrix=all_to_all_matrix;
    data_struct.all_to_all_correlation=all_to_all_correlation;
    data_struct.all_to_all_correlation_multi=all_to_all_correlation_multi;
    data_struct.all_neighbor_correlations=all_neighbor_correlations;
    data_struct.NN_correlations=NN_correlations;
    data_struct.NNN_correlations=NNN_correlations;
    data_struct.all_to_all_matrix_multi=all_to_all_matrix_multi;
    data_struct.all_to_all_distance_multi=all_to_all_distance_multi;
    data_struct.all_neighbor_distances=all_neighbor_distances;
    data_struct.all_neighbor_distances_x=all_neighbor_distances_x;
    data_struct.all_neighbor_distances_y=all_neighbor_distances_y;
    data_struct.NN_distances=NN_distances;
    data_struct.NNN_distances=NNN_distances;
    handles.data_struct=data_struct;
    guidata(hObject, handles)
else % if the distributions were already estimated
    all_to_all_correlation_multi=data_struct.all_to_all_correlation_multi;
    all_neighbor_correlations=data_struct.all_neighbor_correlations;
    NN_correlations=data_struct.NN_correlations;
    NNN_correlations=data_struct.NNN_correlations;
    all_to_all_matrix_multi=data_struct.all_to_all_matrix_multi;
    all_to_all_distance_multi=data_struct.all_to_all_distance_multi;
    all_neighbor_distances=data_struct.all_neighbor_distances;
    all_neighbor_distances_x=data_struct.all_neighbor_distances_x;
    all_neighbor_distances_y=data_struct.all_neighbor_distances_y;
    NN_distances=data_struct.NN_distances;
    NNN_distances=data_struct.NNN_distances;
end

if get(handles.one_photon,'value')==1
    imaging_technique='one_photon';
else
    imaging_technique='two_photon';
end
data_struct.imaging_technique=imaging_technique;

all_neighbor_correlations_temp=all_neighbor_correlations;
all_neighbor_correlations_temp(all_neighbor_correlations_temp==0)=10^-10;
all_neighbor_correlations_temp(all_neighbor_correlations_temp==1)=1-10^-10;
all_neighbor_distances_temp=all_neighbor_distances;
all_neighbor_distances_temp(all_neighbor_distances_temp==0)=10^-10;

all_neighbor_correlations_temp(all_neighbor_distances>temp_dist_thresh)=[];
all_neighbor_distances_temp(all_neighbor_distances>temp_dist_thresh)=[];

if strcmp(imaging_technique,'one_photon');
    if sum(all_neighbor_correlations_temp<0)/length(all_neighbor_correlations_temp)>0.05
        errordlg('The cells seem to be smaller than expected. Either the micron/pixel ratio is incorrect or you should lower the maximal distance')
        error_variable=non_existing_variable;
    else
        all_neighbor_distances_temp(all_neighbor_correlations_temp<0)=[];
        all_neighbor_correlations_temp(all_neighbor_correlations_temp<0)=[];
    end
end


ctrs=cell(1,2);
xout_temp=linspace(0,1,2*num_bins+1);
xcorrout=xout_temp(2:2:end);
ctrs{2}=xcorrout;
xout_temp=linspace(0,temp_dist_thresh,2*num_bins+1);
xout=xout_temp(2:2:end);
ctrs{1}=xout;

grid=hist3([all_neighbor_distances_temp;all_neighbor_correlations_temp]',ctrs);
grid=flipud(grid);

% computing the nearest/ other neighbors distributions:
NN_correlations_temp=NN_correlations;
NNN_correlations_temp=NNN_correlations;

NN_distances_temp=NN_distances;
NNN_distances_temp=NNN_distances;
NN_correlations_temp(NN_distances>temp_dist_thresh)=[];
NN_distances_temp(NN_distances>temp_dist_thresh)=[];
NNN_correlations_temp(NNN_distances>temp_dist_thresh)=[];
NNN_distances_temp(NNN_distances>temp_dist_thresh)=[];
if strcmp(imaging_technique,'one_photon');
    NN_distances_temp(NN_correlations_temp<0)=[];
    NN_correlations_temp(NN_correlations_temp<0)=[];
    NNN_distances_temp(NNN_correlations_temp<0)=[];
    NNN_correlations_temp(NNN_correlations_temp<0)=[];
end
if length(NNN_distances_temp)<0.1*length(NN_distances_temp);
    errordlg('There is insufficient number of non-nearest neighboring cells to estiamte the different cells distribution. Either the micron/pixel ratio is incorrect or you should increase the maximal distance')
    error_variable=non_existing_variable;
end

grid_NN=hist3([NN_distances_temp;NN_correlations_temp]',ctrs);
grid_NNN=hist3([NNN_distances_temp;NNN_correlations_temp]',ctrs);
grid_NN=flipud(grid_NN);
grid_NNN=flipud(grid_NNN);

more_than_06_fraction=sum(NN_correlations_temp>0.6==1)/length(NN_correlations_temp);
less_than_06_fraction=sum(NNN_correlations_temp<0.6==1)/length(NNN_correlations_temp);
less_than_7_fraction=sum(NN_distances_temp<7/pixel_to_mic==1)/length(NN_distances_temp);
more_than_7_fraction=sum(NNN_distances_temp>7/pixel_to_mic==1)/length(NNN_distances_temp);

%(x,y) distnaces distributions:
all_neighbor_distances_x_temp=all_neighbor_distances_x;
all_neighbor_distances_y_temp=all_neighbor_distances_y;
all_neighbor_distances_x_temp(abs(all_neighbor_distances_y_temp)>temp_dist_thresh)=[];
all_neighbor_distances_y_temp(abs(all_neighbor_distances_y_temp)>temp_dist_thresh)=[];
all_neighbor_distances_y_temp(abs(all_neighbor_distances_x_temp)>temp_dist_thresh)=[];
all_neighbor_distances_x_temp(abs(all_neighbor_distances_x_temp)>temp_dist_thresh)=[];

ctrs_xy=cell(1,2);
xout_temp_2=linspace(0,temp_dist_thresh,num_bins+1);
xout_2=xout_temp_2(2:2:end);
ctrs_xy{1}=[-flip(xout_2), xout_2];
ctrs_xy{2}=[-flip(xout_2), xout_2];
grid_x_y=hist3([all_neighbor_distances_x_temp;all_neighbor_distances_y_temp]',ctrs_xy);
grid_x_y=flipud(grid_x_y);
grid_x_y=fliplr(grid_x_y);

figure
fig_size_y=15;
fig_size_x=18;
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 5 fig_size_x fig_size_y]);
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[5 0 0 0]);
size_x=0.75;
size_y=0.75;
axes('position',[0.12 0.15 size_x size_y])
imagesc(log(1+grid_x_y)./max(max(log(1+grid_x_y))))
axis square
title('Centroid distances','fontsize',22)
cmap_jet=colormap('jet');
y=round(linspace(1,num_bins,9));
y_label=round(linspace(pixel_to_mic*max(ctrs{1}),-pixel_to_mic*max(ctrs{1}),9));
x=round(linspace(1,num_bins,9));
x_label=round(linspace(-pixel_to_mic*max(ctrs{1}),pixel_to_mic*max(ctrs{1}),9));set(gca,'YTick',y)
set(gca,'YTickLabel',y_label,'fontsize',18)
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',18)
xlabel('x difference (\mum)','FontWeight','Bold','fontsize',22)
ylabel('y difference (\mum)','FontWeight','Bold','fontsize',22)
x_for_circle=0:0.01:2*pi;
hold on
plot(num_bins/2+num_bins/2*4/temp_dist_thresh/pixel_to_mic*sin(x_for_circle),num_bins/2+num_bins/2*4/temp_dist_thresh/pixel_to_mic*cos(x_for_circle),':','color',[1 1 1],'linewidth',4);
hold on
plot(num_bins/2+num_bins/2*8/temp_dist_thresh/pixel_to_mic*sin(x_for_circle),num_bins/2+num_bins/2*8/temp_dist_thresh/pixel_to_mic*cos(x_for_circle),'--','color',[1 1 1],'linewidth',4);
axes('position',[0.855 0.15 0.02 size_y])
num_colors=size(cmap_jet,1);
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim([0 1])
ylim([0 1])
for n=1:num_colors
    hold on
    p=patch([0 1 1 0],[n/num_colors n/num_colors (n-1)/num_colors (n-1)/num_colors],cmap_jet(n,:));
    set(p,'FaceAlpha',1,'EdgeColor','none');
end
text(3.5,0.5,'Number of cell-pairs (log)','fontsize',22,'fontweight','bold','rotation',90,'HorizontalAlignment','Center')
text(1.5,0,'0','fontsize',18,'fontweight','bold','HorizontalAlignment','Left')
text(1.5,1,'Max','fontsize',18,'fontweight','bold','HorizontalAlignment','Left')
set(gca,'fontsize',18)
set(gcf,'PaperPositionMode','auto')
cd(figures_dir);
savefig('Stage 3 - x-y differences')
saveas(gcf,'Stage 3 - x-y differences','tif')
close;

axes(handles.axes2)
imagesc(log(1+grid_x_y)./max(max(log(1+grid_x_y))))
title('Centroid distances','fontsize',16)
colormap('jet');
freezeColors
y=round(linspace(1,num_bins,9));
y_label=round(linspace(pixel_to_mic*max(ctrs{1}),-pixel_to_mic*max(ctrs{1}),9));
x=round(linspace(1,num_bins,9));
x_label=round(linspace(-pixel_to_mic*max(ctrs{1}),pixel_to_mic*max(ctrs{1}),9));
set(gca,'YTick',y)
set(gca,'YTickLabel',y_label,'fontsize',16)
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',16)
xlabel('x difference (\mum)','FontWeight','Bold','fontsize',16)
ylabel('y difference (\mum)','FontWeight','Bold','fontsize',16)
x_for_circle=0:0.01:2*pi;
hold on
plot(num_bins/2+num_bins/2*4/temp_dist_thresh/pixel_to_mic*sin(x_for_circle),num_bins/2+num_bins/2*4/temp_dist_thresh/pixel_to_mic*cos(x_for_circle),':','color',[1 1 1],'linewidth',4);
hold on
plot(num_bins/2+num_bins/2*8/temp_dist_thresh/pixel_to_mic*sin(x_for_circle),num_bins/2+num_bins/2*8/temp_dist_thresh/pixel_to_mic*cos(x_for_circle),'--','color',[1 1 1],'linewidth',4);

% Estimating the distributions of spatial correlations for same cell and different cells:
if strcmp(imaging_technique,'one_photon');
    [n1,~]=hist(all_neighbor_correlations_temp,ctrs{2});
    if sum(n1(ctrs{2}<0.05))>0.04*sum(n1)
        errordlg('The cells seem to be smaller than expected. Either the micron/pixel ratio is incorrect or you should increase the maximal distance')
        error_variable=non_existing_variable;
    end
    
    pdf_betamixture = @(x,p,A1,A2,B1,B2) ...
        p*lognpdf(x,A1,B1) + (1-p)*betapdf(x,A2,B2);
    
    phat_same=lognfit(1-all_neighbor_correlations_temp(all_neighbor_correlations_temp>=0.7));
    phat_diff=betafit(all_neighbor_correlations_temp(all_neighbor_correlations_temp<0.7));
    pStart=0.5;
    AStart=[phat_same(1) phat_diff(1)];
    BStart=[phat_same(2) phat_diff(2)];
    start=[pStart AStart BStart];
    lb = [0 -inf 0 0 0];
    ub = [1 Inf Inf Inf Inf];
    options = statset('MaxIter',1000, 'MaxFunEvals',2000);
    paramEsts=mle(1-all_neighbor_correlations_temp,'pdf',pdf_betamixture, 'start',start, 'lower',lb, 'upper',ub,'options',options);
    p_corr_marg_same=lognpdf(1-ctrs{2},paramEsts(2),paramEsts(4));
    p_corr_marg_same=p_corr_marg_same./sum(p_corr_marg_same)*(num_bins/((ctrs{2}(2)-ctrs{2}(1))+(ctrs{2}(end)-ctrs{2}(1))));
    smoothing_func=sigmf(ctrs{2},[20 min(ctrs{2})+0.25]);
    p_corr_marg_same=p_corr_marg_same.*smoothing_func;
    p_corr_marg_same(1:round(num_bins/10:end))=0;
    p_corr_marg_diff=betapdf(1-ctrs{2},paramEsts(3),paramEsts(5));
    p_corr_marg_diff=p_corr_marg_diff./sum(p_corr_marg_diff)*(num_bins/((ctrs{2}(2)-ctrs{2}(1))+(ctrs{2}(end)-ctrs{2}(1))));
    p_corr_marg_both=paramEsts(1)*p_corr_marg_same+(1-paramEsts(1))*p_corr_marg_diff;
    [n1,~]=hist(all_neighbor_correlations_temp,ctrs{2});
    n_corr=n1./sum(n1)*(num_bins/((ctrs{2}(2)-ctrs{2}(1))+(ctrs{2}(end)-ctrs{2}(1))));
    ind_for_thresh=find(ctrs{2}>0.3 & ctrs{2}<0.95);
    [~,ind_intersection]=min(abs(paramEsts(1)*p_corr_marg_same(ind_for_thresh)-(1-paramEsts(1))*p_corr_marg_diff(ind_for_thresh)));
    thresh_corr_from_intersection=round(100*ctrs{2}(ind_intersection+ind_for_thresh(1)-1))/100;
    set(handles.correlation_threshold,'string',num2str(thresh_corr_from_intersection))
    data_struct.thresh_corr_from_intersection=thresh_corr_from_intersection;
    MSE_corr_model=sum(abs(((n_corr-p_corr_marg_both))*((ctrs{2}(2)-ctrs{2}(1))+(ctrs{2}(end)-ctrs{2}(1)))/num_bins))/2;
    p_same_corr_model=paramEsts(1);
end

% Estimating the distributions of centroid distances for same cell and different cells:
[n1,~]=hist(all_neighbor_distances_temp,ctrs{1});
n1=n1./sum(n1)*(num_bins/(pixel_to_mic*(ctrs{1}(2)-ctrs{1}(1))+pixel_to_mic*(ctrs{1}(end)-ctrs{1}(1))));

max_dist_to_fit=9;
data_to_fit=all_neighbor_distances_temp(pixel_to_mic*all_neighbor_distances_temp<max_dist_to_fit);
parmhat=lognfit(data_to_fit);

if strcmp(imaging_technique,'one_photon');
    optimal_delta=paramEsts(1);
else
    optimal_delta=length(data_to_fit)/length(all_neighbor_distances_temp);
end

p_0=optimal_delta;
c_0=6;
a_0=1;
b_0=(n1(end)-n1(round(num_bins/2)))/(pixel_to_mic*(ctrs{1}(end)-ctrs{1}(round(num_bins/2))));
b_0=b_0/(1-p_0);
params_init=[p_0 parmhat a_0 c_0 b_0];
F = @(x,xdata)...
    x(1)*(1./(xdata.*x(3).*sqrt(2*pi)).*exp(-(log(xdata)-x(2)).^2./(2*x(3)^2))...
    + (1-x(1)).*x(6).*xdata./(1+exp(-x(4).*(xdata-x(5)))));

lb = [0 -inf 0 0 0 0];
ub = [1 Inf Inf Inf Inf inf];
options = statset('MaxIter',1000, 'MaxFunEvals',2000);
params_final=lsqcurvefit(F,params_init,pixel_to_mic*ctrs{1},n1,lb,ub,options);

p_dist_marg_same=lognpdf(pixel_to_mic*ctrs{1},params_final(2),params_final(3));
p_dist_marg_same=p_dist_marg_same./sum(p_dist_marg_same)*(num_bins/(pixel_to_mic*(ctrs{1}(2)-ctrs{1}(1))+pixel_to_mic*(ctrs{1}(end)-ctrs{1}(1))));
smoothing_func=1-sigmf(pixel_to_mic*ctrs{1},[1.5 pixel_to_mic*temp_dist_thresh-3]);
p_dist_marg_same=p_dist_marg_same.*smoothing_func;
p_dist_marg_diff=params_final(6)*pixel_to_mic*ctrs{1}./(1+exp(-params_final(4)*(pixel_to_mic*ctrs{1}-params_final(5))));
p_dist_marg_diff=p_dist_marg_diff./sum(p_dist_marg_diff)*(num_bins/(pixel_to_mic*(ctrs{1}(2)-ctrs{1}(1))+pixel_to_mic*(ctrs{1}(end)-ctrs{1}(1))));
p_dist_both=params_final(1)*p_dist_marg_same+(1-params_final(1))*p_dist_marg_diff;
normalized_distance=n1;
ind_for_thresh=find(ctrs{1}>1/pixel_to_mic & ctrs{1}<10/pixel_to_mic);
[~,ind_intersection]=min(abs(optimal_delta*p_dist_marg_same(ind_for_thresh)-(1-optimal_delta)*p_dist_marg_diff(ind_for_thresh)));
thresh_dist_from_intersection=round(100*pixel_to_mic*ctrs{1}(ind_intersection+ind_for_thresh(1)-1))/100;
set(handles.distance_threshold,'string',num2str(thresh_dist_from_intersection))
data_struct.thresh_dist_from_intersection=thresh_dist_from_intersection;
MSE_dist_model=sum(abs(((normalized_distance-p_dist_both))*(pixel_to_mic*(ctrs{1}(2)-ctrs{1}(1))+pixel_to_mic*(ctrs{1}(end)-ctrs{1}(1)))/num_bins))/2;
p_same_dist_model=params_final(1);

p_value_dist=1-params_final(1).*p_dist_marg_same./(params_final(1).*p_dist_marg_same+(1-params_final(1)).*p_dist_marg_diff);
p_value_dist(1)=p_value_dist(2);
    
figure('units','normalized','outerposition',[0.15 0.04 0.7 0.96])
if strcmp(imaging_technique,'one_photon');
    subplot(2,2,2)
    p_value_corr=1-paramEsts(1).*p_corr_marg_same./(paramEsts(1).*p_corr_marg_same+(1-paramEsts(1)).*p_corr_marg_diff);
    [n1,~]=hist(NN_correlations_temp,ctrs{2});
    [n2,~]=hist(NNN_correlations_temp,ctrs{2});
    bar(ctrs{2}+0.25/num_bins,n1,'g','EdgeColor','none','barwidth',0.5);
    hold on
    bar(ctrs{2}-0.25/num_bins,n2,'r','EdgeColor','none','barwidth',0.5);
    hold on
    plot([0.6 0.6],[0 max(n1)],'--','linewidth',2,'color','k')
    xlim([0 1])
    legend('Nearest neighbors','Other neighbors','location','northwest')
    legend('boxoff')
    x_label=linspace(0,1,6);
    x=linspace(0,1,6);
    set(gca,'fontsize',16)
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',16)
    set(gcf,'PaperPositionMode','auto')
    xlabel('Spatial correlation','FontWeight','Bold','fontsize',16)
    ylabel('Number of cell-pairs','FontWeight','Bold','fontsize',16)
    text(0.5,0.9*max(n1),[num2str(round(100*(1-more_than_06_fraction))) '%'],'fontsize',14,'fontweight','bold','HorizontalAlignment','Center','color','g')
    text(0.7,0.9*max(n1),[num2str(round(100*more_than_06_fraction)) '%'],'fontsize',14,'fontweight','bold','HorizontalAlignment','Center','color','g')
    text(0.5,0.8*max(n1),[num2str(round(100*less_than_06_fraction)) '%'],'fontsize',16,'fontweight','bold','HorizontalAlignment','Center','color','r')
    text(0.7,0.8*max(n1),[num2str(round(100*(1-less_than_06_fraction))) '%'],'fontsize',16,'fontweight','bold','HorizontalAlignment','Center','color','r')
    set(gca,'fontsize',16)
end
subplot(2,2,1)
[n1,~]=hist(NN_distances_temp,ctrs{1});
[n2,~]=hist(NNN_distances_temp,ctrs{1});
bar(pixel_to_mic*ctrs{1}+0.25*pixel_to_mic*temp_dist_thresh/num_bins,n1,'g','EdgeColor','none','barwidth',0.5);
hold on
bar(pixel_to_mic*ctrs{1}-0.25*pixel_to_mic*temp_dist_thresh/num_bins,n2,'r','EdgeColor','none','barwidth',0.5);
hold on
plot([7 7],[0 max(n1)],'--','linewidth',2,'color','k')
xlim([0 pixel_to_mic*temp_dist_thresh])
x_label=0:3:pixel_to_mic*temp_dist_thresh;
x=0:3:pixel_to_mic*temp_dist_thresh;
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',16)
xlabel('Centroids distance (\mum)','FontWeight','Bold','fontsize',16)
ylabel('Number of cell-pairs','FontWeight','Bold','fontsize',16)
text(6.2,0.9*max(n1),[num2str(round(100*less_than_7_fraction)) '%'],'fontsize',16,'fontweight','bold','HorizontalAlignment','Center','color','g')
text(8,0.9*max(n1),[num2str(round(100*(1-less_than_7_fraction))) '%'],'fontsize',16,'fontweight','bold','HorizontalAlignment','Center','color','g')
text(6,0.8*max(n1),[num2str(round(100*(1-more_than_7_fraction))) '%'],'fontsize',16,'fontweight','bold','HorizontalAlignment','Center','color','r')
text(8,0.8*max(n1),[num2str(round(100*more_than_7_fraction)) '%'],'fontsize',16,'fontweight','bold','HorizontalAlignment','Center','color','r')
set(gca,'fontsize',16)
if strcmp(imaging_technique,'one_photon');
    subplot(2,2,4)
    bar(ctrs{2},n_corr,'FaceColor','b','EdgeColor','none','barwidth',1);
    hold on
    plot(ctrs{2},paramEsts(1)*p_corr_marg_same,'--','linewidth',3,'color','g');
    hold on
    plot(ctrs{2},(1-paramEsts(1))*p_corr_marg_diff,'--','linewidth',3,'color','r');
    hold on
    plot(ctrs{2},p_corr_marg_both,'linewidth',3,'color','k');
    hold on
    plot(ctrs{2},paramEsts(1)*p_corr_marg_same,'--','linewidth',3,'color','g');
    hold on
    plot(ctrs{2},(1-paramEsts(1))*p_corr_marg_diff,'--','linewidth',3,'color','r');
    hold on
    plot([thresh_corr_from_intersection thresh_corr_from_intersection],[0 max(n_corr)],'--','linewidth',2,'color','k')
    xlim([0 1])
    xlabel('Spatial correlation','fontsize',18,'fontweight','bold')
    ylabel('Probability density','fontsize',18,'fontweight','bold')
    hold on
    x_label=linspace(0,1,6);
    x=linspace(0,1,6);
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',16)
    legend('Observed data','Same cell model','Different cells model','Overall model','location','northwest')
    legend('boxoff')
    xlim([0 1])
    set(gca,'fontsize',16)
    normalized_same=p_corr_marg_same./sum(p_corr_marg_same);
    normalized_diff=p_corr_marg_diff./sum(p_corr_marg_diff);
    same_more_than_thresh=sum(normalized_same(ctrs{2}>thresh_corr_from_intersection));
    diff_more_than_thresh=sum(normalized_diff(ctrs{2}>thresh_corr_from_intersection));
    text(thresh_corr_from_intersection+0.1,0.9*max(n_corr),[num2str(round(100*same_more_than_thresh)) '%'],'fontsize',16,'fontweight','bold','HorizontalAlignment','Center','color','g')
    text(thresh_corr_from_intersection-0.1,0.9*max(n_corr),[num2str(round(100*(1-same_more_than_thresh))) '%'],'fontsize',16,'fontweight','bold','HorizontalAlignment','Center','color','g')
    text(thresh_corr_from_intersection+0.1,0.8*max(n_corr),[num2str(round(100*(diff_more_than_thresh))) '%'],'fontsize',16,'fontweight','bold','HorizontalAlignment','Center','color','r')
    text(thresh_corr_from_intersection-0.1,0.8*max(n_corr),[num2str(round(100*(1-diff_more_than_thresh))) '%'],'fontsize',16,'fontweight','bold','HorizontalAlignment','Center','color','r')
end
subplot(2,2,3)
bar(pixel_to_mic*ctrs{1},normalized_distance,'FaceColor','b','EdgeColor','none','barwidth',1);
hold on
plot(pixel_to_mic*ctrs{1},params_final(1)*p_dist_marg_same,'--','color','g','linewidth',3)
hold on
plot(pixel_to_mic*ctrs{1},(1-params_final(1))*p_dist_marg_diff,'--','color','r','linewidth',3)
hold on
plot(pixel_to_mic*ctrs{1},p_dist_both,'color','r','linewidth',3,'color','k')
hold on
plot(pixel_to_mic*ctrs{1},params_final(1)*p_dist_marg_same,'--','color','g','linewidth',3)
hold on
plot(pixel_to_mic*ctrs{1},(1-params_final(1))*p_dist_marg_diff,'--','color','r','linewidth',3)
hold on
plot([thresh_dist_from_intersection thresh_dist_from_intersection],[0 max(normalized_distance)],'--','linewidth',2,'color','k')
xlim([0 pixel_to_mic*temp_dist_thresh])
xlabel('Centroids distance (\mum)','FontWeight','Bold','fontsize',18)
ylabel('Probability density','FontWeight','Bold','fontsize',18)
x_label=0:3:pixel_to_mic*temp_dist_thresh;
x=0:3:pixel_to_mic*temp_dist_thresh;
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',16)
set(gca,'fontsize',16)
normalized_same=p_dist_marg_same./sum(p_dist_marg_same);
normalized_diff=p_dist_marg_diff./sum(p_dist_marg_diff);
same_more_than_thresh=sum(normalized_same(ctrs{1}*pixel_to_mic>thresh_dist_from_intersection));
diff_more_than_thresh=sum(normalized_diff(ctrs{1}*pixel_to_mic>thresh_dist_from_intersection));
text(thresh_dist_from_intersection+1,0.9*max(normalized_distance),[num2str(round(100*same_more_than_thresh)) '%'],'fontsize',16,'fontweight','bold','HorizontalAlignment','Center','color','g')
text(thresh_dist_from_intersection-1,0.9*max(normalized_distance),[num2str(round(100*(1-same_more_than_thresh))) '%'],'fontsize',16,'fontweight','bold','HorizontalAlignment','Center','color','g')
text(thresh_dist_from_intersection+1,0.8*max(normalized_distance),[num2str(round(100*(diff_more_than_thresh))) '%'],'fontsize',16,'fontweight','bold','HorizontalAlignment','Center','color','r')
text(thresh_dist_from_intersection-1,0.8*max(normalized_distance),[num2str(round(100*(1-diff_more_than_thresh))) '%'],'fontsize',16,'fontweight','bold','HorizontalAlignment','Center','color','r')
set(gcf,'PaperPositionMode','auto')
savefig('Stage 3 - model')
saveas(gcf,'Stage 3 - model','tif')

% Computing an estimated ROC curve and fraction of uncertain cell-pairs in 1D:
num_p_values=1000;
p_thresh=0.05;
p_val_vec_temp=linspace(0,1,2*num_p_values+1);
p_val_vec=p_val_vec_temp(2:2:end);
if strcmp(imaging_technique,'one_photon');
    true_merge_corr=zeros(1,length(p_val_vec));
    false_merge_corr=zeros(1,length(p_val_vec));
    num_pairs_corr=zeros(1,length(p_val_vec));
end
true_merge_dist=zeros(1,length(p_val_vec));
false_merge_dist=zeros(1,length(p_val_vec));
num_pairs_dist=zeros(1,length(p_val_vec));

step=p_val_vec(2)-p_val_vec(1);
if strcmp(imaging_technique,'one_photon');
    [n1_corr,~]=hist(all_neighbor_correlations_temp,ctrs{2});
end
[n1_dist,~]=hist(all_neighbor_distances_temp,ctrs{1});
for n=1:length(p_val_vec)
    if strcmp(imaging_technique,'one_photon');        
        p_corr_marg_same_for_ROC=p_corr_marg_same./sum(p_corr_marg_same);
        p_corr_marg_diff_for_ROC=p_corr_marg_diff./sum(p_corr_marg_diff);
        true_merge_corr(n)=sum(p_corr_marg_same_for_ROC(p_value_corr<=p_val_vec(n)+step/2 & p_value_corr>=p_val_vec(n)-step/2));
        false_merge_corr(n)=sum(p_corr_marg_diff_for_ROC(p_value_corr<=p_val_vec(n)+step/2 & p_value_corr>=p_val_vec(n)-step/2));
        num_pairs_corr(n)=sum(n1_corr(1-p_value_corr<=p_val_vec(n)+step/2 & 1-p_value_corr>=p_val_vec(n)-step/2));
    end
    p_dist_marg_same_for_ROC=p_dist_marg_same./sum(p_dist_marg_same);
    p_dist_marg_diff_for_ROC=p_dist_marg_diff./sum(p_dist_marg_diff);
    true_merge_dist(n)=sum(p_dist_marg_same_for_ROC(p_value_dist<=p_val_vec(n)+step/2 & p_value_dist>=p_val_vec(n)-step/2));
    false_merge_dist(n)=sum(p_dist_marg_diff_for_ROC(p_value_dist<=p_val_vec(n)+step/2 & p_value_dist>=p_val_vec(n)-step/2));
    num_pairs_dist(n)=sum(n1_dist(1-p_value_dist<=p_val_vec(n)+step/2 & 1-p_value_dist>=p_val_vec(n)-step/2));
end

if strcmp(imaging_technique,'one_photon');
    uncertain_corr=(sum(num_pairs_corr(1-p_val_vec>p_thresh & 1-p_val_vec<1-p_thresh)))/sum(num_pairs_corr);
    [~,best_ind]=min(cumsum(false_merge_corr).^2+(1-cumsum(true_merge_corr)).^2);
    cumsum_false_merge_corr=cumsum(false_merge_corr);
    cumsum_false_split_corr=1-cumsum(true_merge_corr);
    false_positive_corr=cumsum_false_merge_corr(best_ind);
    false_negative_corr=cumsum_false_split_corr(best_ind);
end
uncertain_dist=(sum(num_pairs_dist(1-p_val_vec>p_thresh & 1-p_val_vec<1-p_thresh)))/sum(num_pairs_dist);
[~,best_ind]=min(cumsum(false_merge_dist).^2+(1-cumsum(true_merge_dist)).^2);
cumsum_false_merge_dist=cumsum(false_merge_dist);
cumsum_false_split_dist=1-cumsum(true_merge_dist);
false_positive_dist=cumsum_false_merge_dist(best_ind);
false_negative_dist=cumsum_false_split_dist(best_ind);

% Estimating the correlation distributions of same cell and different cells given a distance:
joint_flag=get(handles.use_joint_model,'Value');
if joint_flag==1 & strcmp(imaging_technique,'one_photon')==1;
end

figure('units','normalized','outerposition',[0.2 0.04 0.6 0.96]);
size_x=0.78;
size_y=0.78;
if strcmp(imaging_technique,'one_photon');
axes('position',[0.08 0.58 size_x/2.3 size_y/2.85/2.3])
start_x=0;
end_x=1;
y_vec=repmat(n_corr,[2 1]);
y_vec=y_vec(:);
x_vec=(ctrs{2}(2:end)+ctrs{2}(1:end-1))/2;
x_vec=repmat(x_vec,[2 1]);
x_vec=[start_x; x_vec(:); end_x];
for run_bins=1:length(x_vec)/2
    current_color=p_value_corr(run_bins)*[1 1 1];
    patch(x_vec([1 1 2 2]+2*(run_bins-1)),[0 [1 1]*n_corr(run_bins) 0],current_color,'EdgeColor',current_color)
    hold on
end
plot(x_vec,y_vec,'k-','linewidth',2)
xlabel('Spatial correlation','fontsize',16,'fontweight','bold')
x_label=linspace(0,1,6);
x=linspace(0,1,6);
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',16,'fontweight','bold')
xlim([0 1])
[~,ind_005]=min(abs(0.05-(1-p_value_corr)));
p_005=ctrs{2}(ind_005);
[~,ind_05]=min(abs(0.5-(1-p_value_corr)));
p_05=ctrs{2}(ind_05);
[~,ind_095]=min(abs(0.95-(1-p_value_corr)));
p_095=ctrs{2}(ind_095);
hold on
plot([p_005 p_005],[0 max(n_corr)],'--','linewidth',3,'color','k')
hold on
plot([p_05 p_05],[0 max(n_corr)],'--','linewidth',3,'color','k')
hold on
plot([p_095 p_095],[0 max(n_corr)],'--','linewidth',3,'color','k')
text(p_005,1.1*max(n_corr),'0.05','fontsize',14,'fontweight','bold','HorizontalAlignment','Center')
text(p_095,1.1*max(n_corr),'0.95','fontsize',14,'fontweight','bold','HorizontalAlignment','Center')
text(p_05,1.1*max(n_corr),'0.5','fontsize',14,'fontweight','bold','HorizontalAlignment','Center')
end

axes('position',[0.08 0.8 size_x/2.3 size_y/2.85/2.3])
start_x=0;
end_x=temp_dist_thresh*pixel_to_mic;
y_vec=repmat(normalized_distance,[2 1]);
y_vec=y_vec(:);
x_vec=(pixel_to_mic*ctrs{1}(2:end)+pixel_to_mic*ctrs{1}(1:end-1))/2;
x_vec=repmat(x_vec,[2 1]);
x_vec=[start_x; x_vec(:); end_x];
for run_bins=1:length(x_vec)/2
    current_color=p_value_dist(run_bins)*[1 1 1];
    patch(x_vec([1 1 2 2]+2*(run_bins-1)),[0 [1 1]*normalized_distance(run_bins) 0],current_color,'EdgeColor',current_color)
    hold on
end
plot(x_vec,y_vec,'k-','linewidth',2)
xlim([0 pixel_to_mic*temp_dist_thresh])
xlabel('Centroids distance (\mum)','FontWeight','Bold','fontsize',16)
x_label=0:3:pixel_to_mic*temp_dist_thresh;
x=0:3:pixel_to_mic*temp_dist_thresh;
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',16,'fontweight','bold')
[~,ind_005]=min(abs(0.05-(1-p_value_dist)));
p_005=pixel_to_mic*ctrs{1}(ind_005);
[~,ind_05]=min(abs(0.5-(1-p_value_dist)));
p_05=pixel_to_mic*ctrs{1}(ind_05);
[~,ind_095]=min(abs(0.95-(1-p_value_dist)));
p_095=pixel_to_mic*ctrs{1}(ind_095);
hold on
plot([p_005 p_005],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
hold on
plot([p_05 p_05],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
hold on
plot([p_095 p_095],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
text(p_005,1.1*max(normalized_distance),'0.05','fontsize',14,'fontweight','bold','HorizontalAlignment','Center')
text(p_095,1.1*max(normalized_distance),'0.95','fontsize',14,'fontweight','bold','HorizontalAlignment','Center')
text(p_05,1.1*max(normalized_distance),'0.5','fontsize',14,'fontweight','bold','HorizontalAlignment','Center')
axes('position',[0.44 0.58 0.03/2.3 size_y/2.3])
color_vec=linspace(1,0,num_bins);
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim([0 1])
ylim([0 1])
for n=1:num_bins
    hold on
    p=patch([0 1 1 0],[n/num_bins n/num_bins (n-1)/num_bins (n-1)/num_bins],[color_vec(n) color_vec(n) color_vec(n)]);
    set(p,'FaceAlpha',1,'EdgeColor','none');
end
text(2.5,0.5,'P_s_a_m_e','fontsize',18,'fontweight','bold','rotation',90,'HorizontalAlignment','Center')
text(1.5,0,'0','fontsize',18,'fontweight','bold','HorizontalAlignment','Left')
text(1.5,1,'1','fontsize',18,'fontweight','bold','HorizontalAlignment','Left')
text(-31.5,0.5,'Probability density','fontsize',18,'fontweight','bold','rotation',90,'HorizontalAlignment','Center')
plot([0 1 1 0 0],[0 0 1 1 0],'color','k')

if joint_flag==1 & strcmp(imaging_technique,'one_photon')==1;
    axes('position',[0.59 0.58 size_x/2.3 size_y/2.3])
    norm_factor_grid=1;
    log_grid=log(1+grid/norm_factor_grid);
    log_grid=log_grid./max(max(log_grid));
    imagesc(log_grid)
    cmap_jet=colormap('jet');
    y=round(linspace(1,num_bins,5));
    y_label=round(10*linspace(pixel_to_mic*temp_dist_thresh,0,5))/10;
    x=round(linspace(1,num_bins,6));
    x_label=linspace(0,1,6);
    set(gca,'YTick',y)
    set(gca,'YTickLabel',y_label,'fontsize',16,'FontWeight','Bold')
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',16,'FontWeight','Bold')
    xlabel('Spatial correlation','FontWeight','Bold','fontsize',20)
    ylabel('Centroids distance (\mum)','FontWeight','Bold','fontsize',20)
    freezeColors
    hold on
    contour_values=0.05:0.05:0.95;
    contour(p_value,contour_values,'linewidth',1)
    cmap_gray=flipud(colormap('gray'));
    hold on
    contour(p_value,[0.05 0.5 0.95],'linewidth',3)
    set(gca,'fontsize',16)
    axes('position',[0.935 0.58 0.03/2.3 size_y/2.3])
    num_colors=size(cmap_gray,1);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    xlim([0 1])
    ylim([0 1])
    for n=1:num_colors
        hold on
        p=patch([0 1 1 0],[n/num_colors n/num_colors (n-1)/num_colors (n-1)/num_colors],cmap_gray(n,:));
        set(p,'FaceAlpha',1,'EdgeColor','none');
    end
    text(3,0.5,'P_s_a_m_e','fontsize',20,'fontweight','bold','rotation',90,'HorizontalAlignment','Center')
    text(1.5,0,'0','fontsize',18,'fontweight','bold','HorizontalAlignment','Left')
    text(1.5,1,'1','fontsize',18,'fontweight','bold','HorizontalAlignment','Left')
    axes('position',[0.59 0.93 size_x/2.3 0.02/1.5])
    num_colors=size(cmap_gray,1);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    xlim([0 1])
    ylim([0 1])
    for n=1:size(cmap_jet,1)
        hold on
        p=patch([n/num_colors n/num_colors (n-1)/num_colors (n-1)/num_colors],[0 1 1 0],cmap_jet(n,:));
        set(p,'FaceAlpha',1,'EdgeColor','none');
    end
    text(0.5,3.5,'Number of cell-pairs (log)','fontsize',20,'fontweight','bold','HorizontalAlignment','Center')
    text(0,2.5,'0','fontsize',18,'fontweight','bold','HorizontalAlignment','Center')
    text(1,2.5,'Max','fontsize',18,'fontweight','bold','HorizontalAlignment','Center')
end

axes('position',[0.08 0.08 size_x/2.3 size_y/2.3])
plot([0,cumsum(num_pairs_dist)/sum(num_pairs_dist),1],[0,p_val_vec,1],'linewidth',2,'color','r')
hold on
if strcmp(imaging_technique,'one_photon');
    plot([0,cumsum(num_pairs_corr)/sum(num_pairs_corr),1],[0,p_val_vec,1],'linewidth',2,'color','b')
    hold on
end
if joint_flag==1 & strcmp(imaging_technique,'one_photon')==1;  
    plot([0,cumsum(num_pairs)/sum(num_pairs),1],[0,p_val_vec,1],'linewidth',2,'color',[0.1 0.7 0.1])
    hold on
    plot([1-sum(num_pairs(1-p_val_vec<p_thresh))/sum(num_pairs),1-sum(num_pairs(1-p_val_vec<p_thresh))/sum(num_pairs)],[0,1],'--','linewidth',2,'color','k')
    hold on
    plot([sum(num_pairs(1-p_val_vec>1-p_thresh))/sum(num_pairs),sum(num_pairs(1-p_val_vec>1-p_thresh))/sum(num_pairs)],[0,1],'--','linewidth',2,'color','k')
elseif strcmp(imaging_technique,'one_photon');
    plot([1-sum(num_pairs_corr(1-p_val_vec<p_thresh))/sum(num_pairs_corr),1-sum(num_pairs_corr(1-p_val_vec<p_thresh))/sum(num_pairs_corr)],[0,1],'--','linewidth',2,'color','k')
    hold on
    plot([sum(num_pairs_corr(1-p_val_vec>1-p_thresh))/sum(num_pairs_corr),sum(num_pairs_corr(1-p_val_vec>1-p_thresh))/sum(num_pairs_corr)],[0,1],'--','linewidth',2,'color','k')
end
xlabel('Fraction of cell pairs ','fontsize',16,'fontweight','bold')
ylabel('P_s_a_m_e','fontsize',16,'fontweight','bold')
hold on
plot([0 1],[p_thresh p_thresh],'--','linewidth',2,'color','k')
hold on
plot([0 1],[1-p_thresh 1-p_thresh],'--','linewidth',2,'color','k')
if joint_flag==1 & strcmp(imaging_technique,'one_photon')==1;
    legend('Dist.','Corr.','Joint','Location','Northeast')
elseif strcmp(imaging_technique,'one_photon');
    legend('Dist.','Corr.','Location','Northeast')
end
legend('boxoff')
y=linspace(0,1,6);
y_label=linspace(0,1,6);
x=linspace(0,1,6);
x_label=linspace(0,1,6);
set(gca,'YTick',y)
set(gca,'YTickLabel',y_label,'fontsize',16,'fontweight','bold')
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',16,'fontweight','bold')
axes('position',[0.3 0.15 size_x/2.3/3 size_y/2.3/3.5])
labels=cell(1,3);
labels{1}='Dist.';
labels{2}='Corr.';
if joint_flag==1 & strcmp(imaging_technique,'one_photon')==1;
    labels{3}='Joint';
    bar([1 2 3],[uncertain_dist; uncertain_corr ; uncertain],0.8,'FaceColor','none','EdgeColor','k')
    xtick_vec=1:3;
elseif strcmp(imaging_technique,'one_photon');    
    bar([1 2],[uncertain_dist; uncertain_corr],0.8,'FaceColor','none','EdgeColor','k')
    xtick_vec=1:2;
else
    bar(1,uncertain_dist,0.8,'FaceColor','none','EdgeColor','k')
    xtick_vec=1;
end
box off
ylabel('Uncertain fraction ','FontSize',10,'fontweight','bold')
set(gca,'XTick',xtick_vec)
set(gca,'XTickLabel',labels,'FontSize',10,'fontweight','bold')

axes('position',[0.59 0.08 size_x/2.3 size_y/2.3])
if strcmp(imaging_technique,'one_photon');
    plot(cumsum(false_merge_corr),cumsum(true_merge_corr),'linewidth',2,'color','b')
    hold on
end
plot([0,cumsum(false_merge_dist)],[0,cumsum(true_merge_dist)],'linewidth',2,'color','r')
hold on
if joint_flag==1 & strcmp(imaging_technique,'one_photon')==1;
    plot(cumsum(false_merge),cumsum(true_merge),'linewidth',2,'color',[0.1 0.7 0.1])
    hold on
end
plot([0 0.265],[0.9 0.14],'--','color',[0.8 0.8 0.8],'linewidth',2)
hold on
plot([0.1 0.77],[1 0.63],'--','color',[0.8 0.8 0.8],'linewidth',2)
ylabel('True positive rate','fontsize',16,'fontweight','bold')
xlabel('False positive rate','fontsize',16,'fontweight','bold')
ylim([0 1])
xlim([0 1])
p=patch([0 0.1 0.1 0],[1 1 0.9 0.9],[0.8 0.8 0.8]);
set(p,'FaceAlpha',0.3,'EdgeColor','none');
y=linspace(0,1,6);
y_label=linspace(0,1,6);
x=linspace(0,1,6);
x_label=linspace(0,1,6);
set(gca,'YTick',y)
set(gca,'YTickLabel',y_label,'fontsize',16,'fontweight','bold')
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',16,'fontweight','bold')
axes('position',[0.682 0.125 size_x/2.3/2 size_y/2.3/2])
plot([0,cumsum(false_merge_dist)],[0,cumsum(true_merge_dist)],'linewidth',2,'color','r')
if strcmp(imaging_technique,'one_photon');
    hold on
    plot(cumsum(false_merge_corr),cumsum(true_merge_corr),'linewidth',2,'color','b')
end
if joint_flag==1 & strcmp(imaging_technique,'one_photon')==1;
    hold on
    plot(cumsum(false_merge),cumsum(true_merge),'linewidth',2,'color',[0.1 0.7 0.1])
end
hold on
cum_sum_temp=[0,cumsum(false_merge_dist)];
cum_sum_temp_2=[0,cumsum(true_merge_dist)];
plot(cum_sum_temp(500),cum_sum_temp_2(500),'*','markersize',8,'linewidth',2,'color','k')
if strcmp(imaging_technique,'one_photon');
    hold on
    cum_sum_temp=[0,cumsum(false_merge_corr)];
    cum_sum_temp_2=[0,cumsum(true_merge_corr)];
    plot(cum_sum_temp(500),cum_sum_temp_2(500),'*','markersize',8,'linewidth',2,'color','k')
end
if joint_flag==1 & strcmp(imaging_technique,'one_photon')==1;
    hold on
    cum_sum_temp=[0,cumsum(false_merge)];
    cum_sum_temp_2=[0,cumsum(true_merge)];
    plot(cum_sum_temp(500),cum_sum_temp_2(500),'*','markersize',8,'linewidth',2,'color','k')
end

xlim([0 0.1])
ylim([0.9 1])
x_label=0:0.1:0.1;
x=0:0.1:0.2;
y=0.9:0.1:1;
y_label=0.9:0.1:1;
set(gca,'YTick',y)
set(gca,'YTickLabel',y_label,'fontsize',16,'fontweight','bold')
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',16,'fontweight','bold')
p=patch([0 1 1 0],[1 1 0 0],[0.8 0.8 0.8]);
set(p,'FaceAlpha',0.3,'EdgeColor','none');
set(gcf,'PaperPositionMode','auto')
savefig('Stage 3 - registration certainty')
saveas(gcf,'Stage 3 - registration certainty','tif')

axes(handles.axes3)
start_x=0;
end_x=temp_dist_thresh*pixel_to_mic;
y_vec=repmat(normalized_distance,[2 1]);
y_vec=y_vec(:);
x_vec=(pixel_to_mic*ctrs{1}(2:end)+pixel_to_mic*ctrs{1}(1:end-1))/2;
x_vec=repmat(x_vec,[2 1]);
x_vec=[start_x; x_vec(:); end_x];
for run_bins=1:length(x_vec)/2
    current_color=p_value_dist(run_bins)*[1 1 1];
    patch(x_vec([1 1 2 2]+2*(run_bins-1)),[0 [1 1]*normalized_distance(run_bins) 0],current_color,'EdgeColor',current_color)
    hold on
end
plot(x_vec,y_vec,'k-','linewidth',2)
xlim([0 pixel_to_mic*temp_dist_thresh])
xlabel('Centroids distance (\mum)','FontWeight','Bold','fontsize',12)
x_label=0:3:pixel_to_mic*temp_dist_thresh;
x=0:3:pixel_to_mic*temp_dist_thresh;
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
[~,ind_005]=min(abs(0.05-(1-p_value_dist)));
p_005=pixel_to_mic*ctrs{1}(ind_005);
[~,ind_05]=min(abs(0.5-(1-p_value_dist)));
p_05=pixel_to_mic*ctrs{1}(ind_05);
[~,ind_095]=min(abs(0.95-(1-p_value_dist)));
p_095=pixel_to_mic*ctrs{1}(ind_095);
hold on
plot([p_005 p_005],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
hold on
plot([p_05 p_05],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
hold on
plot([p_095 p_095],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
text(p_005,1.1*max(normalized_distance),'0.05','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
text(p_095,1.1*max(normalized_distance),'0.95','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
text(p_05,1.1*max(normalized_distance),'0.5','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')

if strcmp(imaging_technique,'one_photon');
    axes(handles.axes4)
    start_x=0;
    end_x=1;
    y_vec=repmat(n_corr,[2 1]);
    y_vec=y_vec(:);
    x_vec=(ctrs{2}(2:end)+ctrs{2}(1:end-1))/2;
    x_vec=repmat(x_vec,[2 1]);
    x_vec=[start_x; x_vec(:); end_x];
    for run_bins=1:length(x_vec)/2
        current_color=p_value_corr(run_bins)*[1 1 1];
        patch(x_vec([1 1 2 2]+2*(run_bins-1)),[0 [1 1]*n_corr(run_bins) 0],current_color,'EdgeColor',current_color)
        hold on
    end
    plot(x_vec,y_vec,'k-','linewidth',2)
    xlabel('Spatial correlation','fontsize',12,'fontweight','bold')
    x_label=linspace(0,1,6);
    x=linspace(0,1,6);
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
    xlim([0 1])
    [~,ind_005]=min(abs(0.05-(1-p_value_corr)));
    p_005=ctrs{2}(ind_005);
    [~,ind_05]=min(abs(0.5-(1-p_value_corr)));
    p_05=ctrs{2}(ind_05);
    [~,ind_095]=min(abs(0.95-(1-p_value_corr)));
    p_095=ctrs{2}(ind_095);
    hold on
    plot([p_005 p_005],[0 max(n_corr)],'--','linewidth',3,'color','k')
    hold on
    plot([p_05 p_05],[0 max(n_corr)],'--','linewidth',3,'color','k')
    hold on
    plot([p_095 p_095],[0 max(n_corr)],'--','linewidth',3,'color','k')
    text(p_005,1.1*max(n_corr),'0.05','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
    text(p_095,1.1*max(n_corr),'0.95','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
    text(p_05,1.1*max(n_corr),'0.5','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
end

if joint_flag==1 & strcmp(imaging_technique,'one_photon')==1;
    axes(handles.axes5)
    imagesc(log_grid)
    colormap('jet')
    y=round(linspace(1,num_bins,5));
    y_label=round(10*linspace(pixel_to_mic*temp_dist_thresh,0,5))/10;
    x=round(linspace(1,num_bins,6));
    x_label=linspace(0,1,6);
    set(gca,'YTick',y)
    set(gca,'YTickLabel',y_label,'fontsize',14,'FontWeight','Bold')
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',14,'FontWeight','Bold')
    xlabel('Spatial correlation','FontWeight','Bold','fontsize',14)
    ylabel('Centroids distance (\mum)','FontWeight','Bold','fontsize',14)
    freezeColors
    hold on
    contour_values=0.05:0.05:0.95;
    contour(p_value,contour_values,'linewidth',1)
    colormap('gray')
    hold on
    contour(p_value,[0.05 0.5 0.95],'linewidth',3)
    set(gca,'fontsize',14)
    
    axes(handles.axes6)
    cla;
    plot([0,cumsum(num_pairs_dist)/sum(num_pairs_dist),1],[0,p_val_vec,1],'linewidth',2,'color','r')
    hold on
    plot([0,cumsum(num_pairs_corr)/sum(num_pairs_corr),1],[0,p_val_vec,1],'linewidth',2,'color','b')
    hold on
    plot([0,cumsum(num_pairs)/sum(num_pairs),1],[0,p_val_vec,1],'linewidth',2,'color',[0.1 0.7 0.1])
    hold on
    plot([1-sum(num_pairs(1-p_val_vec<p_thresh))/sum(num_pairs),1-sum(num_pairs(1-p_val_vec<p_thresh))/sum(num_pairs)],[0,1],'--','linewidth',2,'color','k')
    hold on
    plot([sum(num_pairs(1-p_val_vec>1-p_thresh))/sum(num_pairs),sum(num_pairs(1-p_val_vec>1-p_thresh))/sum(num_pairs)],[0,1],'--','linewidth',2,'color','k')
    xlabel('Fraction of cell pairs ','fontsize',14,'fontweight','bold')
    ylabel('P_s_a_m_e','fontsize',14,'fontweight','bold')
    hold on
    plot([0 1],[p_thresh p_thresh],'--','linewidth',2,'color','k')
    hold on
    plot([0 1],[1-p_thresh 1-p_thresh],'--','linewidth',2,'color','k')
    y=linspace(0,1,6);
    y_label=linspace(0,1,6);
    x=linspace(0,1,6);
    x_label=linspace(0,1,6);
    set(gca,'YTick',y)
    set(gca,'YTickLabel',y_label,'fontsize',14)
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',14)
    legend('Dist.','Corr.','Joint','Location',[0.91 0.32 0.01 0.01])
    legend('boxoff')
    axes('position',[0.91 0.2 0.04 0.05])
    labels=cell(1,3);
    labels{1}='Dist.';
    labels{2}='Corr.';
    labels{3}='Joint';
    bar([1 2 3],100*[uncertain_dist; uncertain_corr ; uncertain],0.8,'FaceColor','none','EdgeColor','k')
    xtick_vec=1:3;
    box off
    ylabel('Uncertain %','FontSize',10,'fontweight','bold')
    set(gca,'XTick',xtick_vec)
    set(gca,'XTickLabel',labels,'FontSize',6,'fontweight','bold')
else
    axes(handles.axes5)
    cla;
    plot([0,cumsum(num_pairs_dist)/sum(num_pairs_dist),1],[0,p_val_vec,1],'linewidth',2,'color','r')
    if strcmp(imaging_technique,'one_photon');
        hold on
        plot([0,cumsum(num_pairs_corr)/sum(num_pairs_corr),1],[0,p_val_vec,1],'linewidth',2,'color','b')
        hold on
        plot([1-sum(num_pairs_corr(1-p_val_vec<p_thresh))/sum(num_pairs_corr),1-sum(num_pairs_corr(1-p_val_vec<p_thresh))/sum(num_pairs_corr)],[0,1],'--','linewidth',2,'color','k')
        hold on
        plot([sum(num_pairs_corr(1-p_val_vec>1-p_thresh))/sum(num_pairs_corr),sum(num_pairs_corr(1-p_val_vec>1-p_thresh))/sum(num_pairs_corr)],[0,1],'--','linewidth',2,'color','k')
    end
    xlabel('Fraction of cell pairs ','fontsize',14,'fontweight','bold')
    ylabel('P_s_a_m_e','fontsize',14,'fontweight','bold')
    hold on
    plot([0 1],[p_thresh p_thresh],'--','linewidth',2,'color','k')
    hold on
    plot([0 1],[1-p_thresh 1-p_thresh],'--','linewidth',2,'color','k')
    y=linspace(0,1,6);
    y_label=linspace(0,1,6);
    x=linspace(0,1,6);
    x_label=linspace(0,1,6);
    set(gca,'YTick',y)
    set(gca,'YTickLabel',y_label,'fontsize',14)
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',14)
    if strcmp(imaging_technique,'one_photon');
        legend('Dist.','Corr.','Location',[0.72 0.32 0.01 0.01])
        legend('boxoff')
    end
    axes('position',[0.72 0.2 0.04 0.05])
    labels=cell(1,3);
    labels{1}='Dist.';
    labels{2}='Corr.';
    if strcmp(imaging_technique,'one_photon');
        bar([1 2],100*[uncertain_dist; uncertain_corr],0.8,'FaceColor','none','EdgeColor','k')
        xtick_vec=1:2;
    else
        bar(1,100*uncertain_dist,0.8,'FaceColor','none','EdgeColor','k')
        xtick_vec=1;
    end
    box off
    ylabel('Uncertain %','FontSize',10,'fontweight','bold')
    set(gca,'XTick',xtick_vec)
    set(gca,'XTickLabel',labels,'FontSize',6,'fontweight','bold')
    
    axes(handles.axes6)
    if strcmp(imaging_technique,'one_photon');
        plot(cumsum(false_merge_corr),cumsum(true_merge_corr),'linewidth',2,'color','b')
        hold on
    end
    plot([0,cumsum(false_merge_dist)],[0,cumsum(true_merge_dist)],'linewidth',2,'color','r')
    ylabel('True positive rate','fontsize',16,'fontweight','bold')
    xlabel('False positive rate','fontsize',16,'fontweight','bold')
    ylim([0 1])
    xlim([0 1])
    set(gca,'fontsize',14)
    set(gcf,'PaperPositionMode','auto')
end

if strcmp(imaging_technique,'one_photon');
    data_struct.false_merge_corr=false_merge_corr;
    data_struct.true_merge_corr=true_merge_corr;
    data_struct.num_pairs_corr=num_pairs_corr;
    data_struct.uncertain_corr=uncertain_corr;
    data_struct.false_positive_corr=false_positive_corr;
    data_struct.false_negative_corr=false_negative_corr;
    data_struct.all_neighbor_correlations_temp=all_neighbor_correlations_temp;
    data_struct.MSE_corr_model=MSE_corr_model;
    data_struct.p_same_corr_model=p_same_corr_model;
    data_struct.p_value_corr=p_value_corr;
    data_struct.more_than_06_fraction=more_than_06_fraction;
    data_struct.less_than_06_fraction=less_than_06_fraction;
    data_struct.n_corr=n_corr;
end

data_struct.false_merge_dist=false_merge_dist;
data_struct.true_merge_dist=true_merge_dist;
data_struct.num_pairs_dist=num_pairs_dist;
data_struct.uncertain_dist=uncertain_dist;
data_struct.false_positive_dist=false_positive_dist;
data_struct.false_negative_dist=false_negative_dist;
data_struct.p_value_dist=p_value_dist;
data_struct.all_neighbor_distances_temp=all_neighbor_distances_temp;
data_struct.MSE_dist_model=MSE_dist_model;
data_struct.p_same_dist_model=p_same_dist_model;
data_struct.less_than_7_fraction=less_than_7_fraction;
data_struct.more_than_7_fraction=more_than_7_fraction;
data_struct.normalized_distance=normalized_distance;

if joint_flag==1 & strcmp(imaging_technique,'one_photon')==1;
    data_struct.false_merge=false_merge;
    data_struct.true_merge=true_merge;
    data_struct.num_pairs=num_pairs;
    data_struct.uncertain=uncertain;
    data_struct.false_positive=false_positive;
    data_struct.false_negative=false_negative;
    data_struct.MSE_2d_model=MSE_2d_model;
    data_struct.p_value=p_value;
    data_struct.log_grid=log_grid;
end
data_struct.ctrs=ctrs;

% Saving the pairwise correlations, distances, and P_same:
if strcmp(imaging_technique,'one_photon');
    all_to_all_p_value_corr_multi=cell(1,num_sessions);
end
all_to_all_p_value_dist_multi=cell(1,num_sessions);
if joint_flag==1 & strcmp(imaging_technique,'one_photon')==1;
    all_to_all_p_value_multi=cell(1,num_sessions);
end

h = waitbar(0,'Computing P_s_a_m_e for all cell-pairs','Units', 'normalized', 'Position',[0.4 0.5 0.2 0.07]);
for n=1:num_sessions
    waitbar((n)/num_sessions,h,['Computing P_s_a_m_e - session number ' num2str(n) '/' num2str(num_sessions)])
    num_cells=size(all_to_all_distance_multi{n},1);
    if strcmp(imaging_technique,'one_photon');
        all_to_all_p_value_corr_multi{n}=cell(num_cells,num_sessions);
    end
    all_to_all_p_value_dist_multi{n}=cell(num_cells,num_sessions);
    if joint_flag==1 & strcmp(imaging_technique,'one_photon')==1;
        all_to_all_p_value_multi{n}=cell(num_cells,num_sessions);
    end
    sessions_to_compare=1:num_sessions;
    sessions_to_compare(n)=[];
    h2 = waitbar(0,'Computing P_s_a_m_e for all cell-pairs','Units', 'normalized', 'Position',[0.4 0.4 0.2 0.07]);
    for k=1:num_cells
        waitbar((k)/num_cells,h2,['Computing P_s_a_m_e - cell number ' num2str(k) '/' num2str(num_cells)]);
        for m=1:length(sessions_to_compare)
            this_session=sessions_to_compare(m);
            if strcmp(imaging_technique,'one_photon');
                temp_corr_vec=all_to_all_correlation_multi{n}{k,this_session};
            end
            temp_dist_vec=all_to_all_distance_multi{n}{k,this_session};
            length_vec=length(temp_dist_vec);
            if length_vec>0
                if strcmp(imaging_technique,'one_photon');                    
                    p_value_corr_vec=zeros(1,length(temp_corr_vec));
                end
                p_value_dist_vec=zeros(1,length(temp_dist_vec));
                if joint_flag==1
                    p_value_vec=zeros(1,length(temp_corr_vec));
                end
                for p=1:length_vec
                    if strcmp(imaging_technique,'one_photon');                        
                        [~,this_corr_ind]=min(abs(ctrs{2}-temp_corr_vec(p)));
                        this_p_value_corr=p_value_corr(this_corr_ind);
                        p_value_corr_vec(p)=this_p_value_corr;
                    end
                    [~,this_dist_ind]=min(abs(ctrs{1}-temp_dist_vec(p)));
                    this_p_value_dist=p_value_dist(this_dist_ind);
                    p_value_dist_vec(p)=this_p_value_dist;
                    
                    if joint_flag==1 & strcmp(imaging_technique,'one_photon')==1;
                        [~,this_dist_ind]=min(abs(ctrs{1}-temp_dist_vec(p)));
                        [~,this_corr_ind]=min(abs(ctrs{2}-temp_corr_vec(p)));
                        this_p_value=p_value(num_bins+1-this_dist_ind,this_corr_ind);
                        p_value_vec(p)=this_p_value;
                    end
                end
                if strcmp(imaging_technique,'one_photon');                    
                    all_to_all_p_value_corr_multi{n}{k,this_session}=1-p_value_corr_vec;
                end
                all_to_all_p_value_dist_multi{n}{k,this_session}=1-p_value_dist_vec;
                if joint_flag==1 & strcmp(imaging_technique,'one_photon')==1;
                    all_to_all_p_value_multi{n}{k,this_session}=1-p_value_vec;
                end
            end
        end
    end
    close(h2);
end
close(h);

if strcmp(imaging_technique,'one_photon');
    data_struct.all_to_all_p_value_corr_multi=all_to_all_p_value_corr_multi;
    data_struct.all_to_all_correlation_multi=all_to_all_correlation_multi;
end
data_struct.all_to_all_p_value_dist_multi=all_to_all_p_value_dist_multi;
data_struct.all_to_all_distance_multi=all_to_all_distance_multi;
if joint_flag==1 & strcmp(imaging_technique,'one_photon')==1;
    data_struct.all_to_all_p_value_multi=all_to_all_p_value_multi;
end
data_struct.all_to_all_matrix_multi=all_to_all_matrix_multi;

if joint_flag==1 & strcmp(imaging_technique,'one_photon')==1;
    if MSE_2d_model>0.25
        warndlg('There is large discrepancy between the joint model and the data')
    end
end
if strcmp(imaging_technique,'one_photon');
if MSE_corr_model>0.1
    warndlg('There is large discrepancy between the spatial correlations model and the data')
end
end
if MSE_dist_model>0.1
    warndlg('There is large discrepancy between the centroid distances model and the data')
end

handles.data_struct=data_struct;
guidata(hObject, handles)
results_dir=data_struct.results_dir;
cd(results_dir);

save('model_data_struct.mat','data_struct','-v7.3')
msgbox('Finished computing probabilistic model')


% --- Executes on button press in register_cells_initial.
function register_cells_initial_Callback(hObject,~,handles)
% hObject    handle to register_cells_initial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Stage 4: Obtaining an initial cell registration according to an optimized
% registration threshold.

% This callback performs initial cell registration according to either
% spatial correlations or centroid distances:
data_struct=handles.data_struct;
  
pixel_to_mic=str2double(get(handles.pixel_to_mic,'string'));
correlation_thresh=str2double(get(handles.correlation_threshold,'string'));
distance_thresh=str2double(get(handles.distance_threshold,'string'))/pixel_to_mic;
maximal_distance_thresh=str2double(get(handles.model_maximal_distance,'string'))/pixel_to_mic;
all_filters_corrected=data_struct.all_filters_corrected;
all_centroids_corrected=data_struct.all_centroids_corrected;
num_sessions=data_struct.num_sessions;
temp_dist_thresh=str2double(get(handles.model_maximal_distance,'string'))/pixel_to_mic;

if isfield(data_struct,'number_of_bins')
    num_bins=data_struct.number_of_bins;
else
    if num_sessions>7
        num_bins=100;
    elseif num_sessions>5
        num_bins=80;
    elseif num_sessions>3
        num_bins=60;
    else
        num_bins=50;
    end
    data_struct.number_of_bins=num_bins;
end

results_dir=data_struct.results_dir;

mkdir('Figures');
figures_dir=fullfile(data_struct.results_dir,'Figures');
data_struct.figures_dir=figures_dir;

if get(handles.spatial_correlations,'Value')==1 % if spatial correlations are used the function "initial_clustering_corr" is called
    if ~isfield(data_struct,'all_to_all_p_value_corr_multi')
        if get(handles.one_photon,'value')==1
            errordlg('Please compute the spatial correlations probability model before performing final cell registration')
            error_variable=non_existing_variable;
        else
            errordlg('The spatial correlations model is only applied for 1-photon imaging data')
            error_variable=non_existing_variable;
        end
    else
    [cell_to_index_map,~,all_neighbor_correlations,all_neighbor_distances,all_assigned_correlations,non_assigned_correlations,~,~]=initial_clustering_corr(maximal_distance_thresh,correlation_thresh,all_filters_corrected,all_centroids_corrected,num_sessions);        
    figure
    xout=linspace(0,1,num_bins);
    [n1,~]=hist(all_assigned_correlations,xout);
    [n2,~]=hist(non_assigned_correlations,xout);
    bar(xout+0.25/num_bins,n1,'g','EdgeColor','none','barwidth',0.5);
    hold on
    bar(xout-0.25/num_bins,n2,'r','EdgeColor','none','barwidth',0.5);
    xlim([0 1])
    x_label=linspace(0,1,6);
    x=linspace(0,1,6);
    set(gca,'fontsize',16)
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',16)
    xlabel('Spatial correlation','FontWeight','Bold','fontsize',16)
    ylabel('Number of cell-pairs','FontWeight','Bold','fontsize',16)
    set(gca,'fontsize',16)
    legend('Same Cell','Different Cells','location','northwest')
    legend('boxoff')   
    end
else % if centroid distances are used it is performed here:
    initial_cell_num=size(all_filters_corrected{1},1);
    cell_to_index_map=zeros(initial_cell_num,num_sessions);
    cell_to_index_map(:,1)=1:initial_cell_num;
    distance_map=zeros(initial_cell_num,num_sessions);
    all_sessions_centroids=all_centroids_corrected{1};
    distance_thresh_high=maximal_distance_thresh;
    distance_thresh_low=distance_thresh;
    count=0;
    duplicate_match_count=0;
    all_neighbor_distances=zeros(1,num_sessions^2*5000);
    all_assigned_distances=zeros(1,num_sessions^2*5000);
    non_assigned_distances=zeros(1,num_sessions^2*5000);
    
    neighbor_count=0;
    assigned_count=0;
    non_assigned_count=0;
    num_candidates=0;
    
    % Registering cells according to the closest neighbors in the order of
    % the sessions:
    for n=2:num_sessions;
        new_centroids=all_centroids_corrected{n};
        for k=1:size(new_centroids,1)
            is_assigned=0;
            centroid=repmat(new_centroids(k,:),size(all_sessions_centroids,1),1);
            all_distances=sqrt(sum((centroid-all_sessions_centroids).^2,2));
            filters_to_check=find(all_distances<distance_thresh_high);
            if ~isempty(filters_to_check)
                num_candidates=num_candidates+length(filters_to_check);
                distance_vec=zeros(1,length(filters_to_check));
                for m=1:length(filters_to_check)
                    neighbor_count=neighbor_count+1;
                    distance_vec(m)=all_distances(filters_to_check(m));
                    all_neighbor_distances(neighbor_count)=all_distances(filters_to_check(m));
                end
                [lowest_dist,lowest_dist_ind]=min(distance_vec);
                if lowest_dist>distance_thresh_low
                    count=count+1;
                    all_sessions_centroids(initial_cell_num+count,:)=new_centroids(k,:);
                    cell_to_index_map(initial_cell_num+count,:)=zeros(1,num_sessions);
                    cell_to_index_map(initial_cell_num+count,n)=k;
                    distance_map(initial_cell_num+count,:)=zeros(1,num_sessions);
                    distance_map(initial_cell_num+count,n)=lowest_dist;
                else
                    index=filters_to_check(lowest_dist_ind);
                    if cell_to_index_map(index,n)==0
                        cell_to_index_map(index,n)=k;
                        distance_map(index,n)=lowest_dist;
                        assigned_count=assigned_count+1;
                        is_assigned=1;
                        all_assigned_distances(1,assigned_count)=lowest_dist;
                    else
                        duplicate_match_count=duplicate_match_count+1;
                        if lowest_dist<distance_map(index,n) % switch between cells
                            count=count+1;
                            switch_cell=cell_to_index_map(index,n);
                            switch_centroid=(all_centroids_corrected{n}(switch_cell,:));
                            all_sessions_centroids(initial_cell_num+count,:)=switch_centroid;
                            cell_to_index_map(initial_cell_num+count,n)=switch_cell;
                            distance_map(initial_cell_num+count,:)=zeros(1,num_sessions);
                            cell_to_index_map(index,n)=k;
                            distance_map(index,n)=lowest_dist;
                            assigned_count=assigned_count+1;
                            is_assigned=1;
                            all_assigned_distances(1,assigned_count)=lowest_dist;
                        else
                            count=count+1;
                            all_sessions_centroids(initial_cell_num+count,:)=new_centroids(k,:);
                            cell_to_index_map(initial_cell_num+count,:)=zeros(1,num_sessions);
                            distance_map(initial_cell_num+count,:)=zeros(1,num_sessions);
                            cell_to_index_map(initial_cell_num+count,n)=k;
                            distance_map(initial_cell_num+count,:)=zeros(1,num_sessions);
                            distance_map(initial_cell_num+count,n)=lowest_dist;
                        end
                    end
                end
                if is_assigned==0
                    non_assigned_count=non_assigned_count+length(filters_to_check);
                    non_assigned_distances(non_assigned_count-length(filters_to_check)+1:non_assigned_count)=distance_vec;
                else
                    temp_distance_vec=distance_vec;
                    temp_distance_vec(lowest_dist_ind)=[];
                    non_assigned_count=non_assigned_count+length(filters_to_check)-1;
                    non_assigned_distances(non_assigned_count-length(temp_distance_vec)+1:non_assigned_count)=temp_distance_vec;
                end
            else
                count=count+1;
                all_sessions_centroids(initial_cell_num+count,:)=new_centroids(k,:);
                cell_to_index_map(initial_cell_num+count,:)=zeros(1,num_sessions);
                cell_to_index_map(initial_cell_num+count,n)=k;
                distance_map(initial_cell_num+count,:)=zeros(1,num_sessions);
            end
        end
    end
    
    all_neighbor_distances(neighbor_count+1:end)=[];
    all_assigned_distances(assigned_count+1:end)=[];
    non_assigned_distances(non_assigned_count+1:end)=[];  
    
    figure
    xout=linspace(0,distance_thresh_high,num_bins);
    [n1,~]=hist(all_assigned_distances,xout);
    [n2,~]=hist(non_assigned_distances,xout);
    bar(pixel_to_mic*xout+0.25*pixel_to_mic*temp_dist_thresh/num_bins,n1,'g','EdgeColor','none','barwidth',0.5);
    hold on
    bar(pixel_to_mic*xout-0.25*pixel_to_mic*temp_dist_thresh/num_bins,n2,'r','EdgeColor','none','barwidth',0.5);
    xlim([0 pixel_to_mic*temp_dist_thresh])
    x_label=0:3:pixel_to_mic*temp_dist_thresh;
    x=0:3:pixel_to_mic*temp_dist_thresh;
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',16)
    xlabel('Centroids distance (\mum)','FontWeight','Bold','fontsize',16)
    ylabel('Number of cell-pairs','FontWeight','Bold','fontsize',16)
    set(gca,'fontsize',16)
    legend('Same Cell','Different Cells','location','northwest')
    legend('boxoff')
end
cd(figures_dir)
set(gcf,'PaperPositionMode','auto')
savefig('Stage 4 - same versus different cells')
saveas(gcf,'Stage 4 - same versus different cells','tif')
close;
    
% Plotting the registration results with the cell maps from all sessions:
filter_thresh=0.2;
all_projections_partial=cell(1,num_sessions);
mutual_projections_partial=cell(1,num_sessions);
cells_in_all_days=find(sum(cell_to_index_map'>0)==num_sessions);
other_cells=cell(1,num_sessions);
for n=1:num_sessions
    logical_1=sum(cell_to_index_map'>0)<num_sessions;
    other_cells{n}=find(cell_to_index_map(:,n)'>0 & logical_1);
end

for n=1:num_sessions
    this_session_filters=all_filters_corrected{n};
    num_filters=size(this_session_filters,1);
    normalized_filters=zeros(size(this_session_filters));
    for k=1:num_filters
        this_filter=this_session_filters(k,:,:);
        this_filter(this_filter<filter_thresh*max(max(this_filter)))=0;
        if max(max(this_filter))>0
            normalized_filters(k,:,:)=this_filter/max(max(this_filter));
        end
    end
    all_projections_partial{n}=zeros(size(this_filter,2),size(this_filter,3),3);
    mutual_projections_partial{n}=zeros(size(this_filter,2),size(this_filter,3),3);
    all_projections_partial{n}(:,:,2)=squeeze(sum(normalized_filters(cell_to_index_map(cells_in_all_days,n),:,:),1));
    all_projections_partial{n}(:,:,1)=squeeze(sum(normalized_filters(cell_to_index_map(other_cells{n},n),:,:),1));
    all_projections_partial{n}(:,:,2)=squeeze(sum(normalized_filters(cell_to_index_map(other_cells{n},n),:,:),1))+squeeze(sum(normalized_filters(cell_to_index_map(cells_in_all_days,n),:,:),1));
    all_projections_partial{n}(:,:,3)=squeeze(sum(normalized_filters(cell_to_index_map(other_cells{n},n),:,:),1));
    mutual_projections_partial{n}(:,:,2)=squeeze(sum(normalized_filters(cell_to_index_map(cells_in_all_days,n),:,:),1));
    all_projections_partial{n}(all_projections_partial{n}>1)=1;
end

subx=4;
suby=ceil(num_sessions/subx);
if num_sessions>4
    figure('units','normalized','outerposition',[0.1 0.04 0.8 0.96])
    for n=1:num_sessions
        subplot(suby,subx,n)
        imagesc(all_projections_partial{n})
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        this_session_name=data_struct.sessions_list{n+1};
        title(strrep(this_session_name,'_','\_'),'fontsize',14,'fontweight','bold')
        if n==1
            text(0.01*size(all_projections_partial{n},1),0.02*size(all_projections_partial{n},2),'Detected in','fontsize',16,'color','g','fontweight','bold')
            text(0.01*size(all_projections_partial{n},1),0.06*size(all_projections_partial{n},2),'all sessions','fontsize',16,'color','g','fontweight','bold')
        end
    end
else
    figure('units','normalized','outerposition',[0.1 0.2 0.9 0.5])
    for n=1:num_sessions
        subplot(1,num_sessions,n)
        imagesc(all_projections_partial{n})
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        colormap('gray')
        this_session_name=data_struct.sessions_list{n+1};
        title(strrep(this_session_name,'_','\_'),'fontsize',14,'fontweight','bold')
        if n==1
            text(0.01*size(all_projections_partial{n},1),0.02*size(all_projections_partial{n},2),'Detected in','fontsize',16,'color','g','fontweight','bold')
            text(0.01*size(all_projections_partial{n},1),0.06*size(all_projections_partial{n},2),'all sessions','fontsize',16,'color','g','fontweight','bold')
        end
    end
end
set(gcf,'PaperPositionMode','auto')
savefig('Stage 4 - projcetions - initial registration')
saveas(gcf,'Stage 4 - projcetions - initial registration','tif')
close;
cd(results_dir)

data_struct.cell_to_index_map=cell_to_index_map;
handles.data_struct=data_struct;
guidata(hObject, handles)
msgbox(['Finished performing initial cell registration - ' num2str(size(cell_to_index_map,1)) ' were found'])


% --- Executes on button press in register_cells_final.
function register_cells_final_Callback(hObject,~, handles)
% hObject    handle to register_cells_final (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Stage 5: Obtaining the final cell registration based on a correlation clustering algorithm.

% This callback performs the final cell registration according to the
% probability model for same cells and different cells. P_same can be
% either according to centroid distances, spatial correlations or both:
data_struct=handles.data_struct;
results_dir=data_struct.results_dir;

% Different options for clustering criterion:
decision_type='Maximal correlation';
% decision_type='Minimal correlation'; 
% decision_type='Average correlation';
data_struct.decision_type=decision_type;

cd(results_dir);

mkdir('Figures');
figures_dir=fullfile(data_struct.results_dir,'Figures');
data_struct.figures_dir=figures_dir;

if ~isfield(data_struct,'cell_to_index_map')
    errordlg('Final registration cannot be performed before initial registration')
end
cell_to_index_map=data_struct.cell_to_index_map;

if get(handles.p_value,'Value')==1
    if isfield(data_struct,'all_to_all_p_value_multi')
        all_to_all_p_value_multi=data_struct.all_to_all_p_value_multi;
        all_to_all_p_value_corr_multi=data_struct.all_to_all_p_value_corr_multi;
        all_to_all_p_value_dist_multi=data_struct.all_to_all_p_value_dist_multi;
        all_to_all_matrix_multi=data_struct.all_to_all_matrix_multi;
        all_to_all_correlation_multi=data_struct.all_to_all_correlation_multi;
        all_to_all_distance_multi=data_struct.all_to_all_distance_multi;
    else
        errordlg('Please compute the joint probability model before performing final cell registration')
    end
elseif  get(handles.spatial_correlations_2,'Value')==1
    if isfield(data_struct,'all_to_all_p_value_corr_multi')
        all_to_all_p_value_corr_multi=data_struct.all_to_all_p_value_corr_multi;
        all_to_all_matrix_multi=data_struct.all_to_all_matrix_multi;
        all_to_all_correlation_multi=data_struct.all_to_all_correlation_multi;
    else
        if get(handles.one_photon,'value')==1
            errordlg('Please compute the spatial correlations probability model before performing final cell registration')
        else
            errordlg('The spatial correlations model is only applied for 1-photon imaging data')
        end
    end
elseif get(handles.centroid_distances_2,'Value')==1
    if isfield(data_struct,'all_to_all_p_value_dist_multi')
        all_to_all_p_value_dist_multi=data_struct.all_to_all_p_value_dist_multi;
        all_to_all_matrix_multi=data_struct.all_to_all_matrix_multi;
        all_to_all_distance_multi=data_struct.all_to_all_distance_multi;
    end
end

all_centroids_corrected=data_struct.all_centroids_corrected;
all_filters_corrected=data_struct.all_filters_corrected;
num_sessions=data_struct.num_sessions;
pixel_to_mic=str2double(get(handles.pixel_to_mic,'string'));
max_iterations=10;
dist_thresh=str2double(get(handles.model_maximal_distance,'string'))/pixel_to_mic;
cluster_distance_thresh=1.7*dist_thresh;
p_value_thresh=str2double(get(handles.decision_thresh,'string'));

is_figure=1;
cd(results_dir);
% Cell registration is performed by the "find_optimal_clustering" function
if get(handles.spatial_correlations_2,'Value')==1 && get(handles.use_model,'Value')==1;
    [optimal_cell_to_index_map,~,~,~,~,register_scores,~,~,~,~,~,~,cell_scores_positive,cell_scores_negative,cell_scores_exclusive,all_clusters_centroids]=find_optimal_clustering(cell_to_index_map,all_to_all_p_value_corr_multi,all_to_all_matrix_multi,max_iterations,cluster_distance_thresh,p_value_thresh,all_centroids_corrected,num_sessions,decision_type,pixel_to_mic,is_figure,results_dir,figures_dir);
elseif get(handles.p_value,'Value')==1 && get(handles.use_model,'Value')==1;
    errordlg('The joint model is not implemented in this version. Please choose a different model')
    [optimal_cell_to_index_map,~,~,~,~,register_scores,~,~,~,~,~,~,cell_scores_positive,cell_scores_negative,cell_scores_exclusive,all_clusters_centroids]=find_optimal_clustering(cell_to_index_map,all_to_all_p_value_multi,all_to_all_matrix_multi,max_iterations,cluster_distance_thresh,p_value_thresh,all_centroids_corrected,num_sessions,decision_type,pixel_to_mic,is_figure,results_dir,figures_dir);
elseif get(handles.centroid_distances_2,'Value')==1 && get(handles.use_model,'Value')==1;
    [optimal_cell_to_index_map,~,~,~,~,register_scores,~,~,~,~,~,~,cell_scores_positive,cell_scores_negative,cell_scores_exclusive,all_clusters_centroids]=find_optimal_clustering(cell_to_index_map,all_to_all_p_value_dist_multi,all_to_all_matrix_multi,max_iterations,cluster_distance_thresh,p_value_thresh,all_centroids_corrected,num_sessions,decision_type,pixel_to_mic,is_figure,results_dir,figures_dir);
elseif get(handles.spatial_correlations_2,'Value')==1 && get(handles.use_simple_threshold,'Value')==1;
    is_figure=0;
    correlation_thresh=str2double(get(handles.simple_threshold,'string'));
    [optimal_cell_to_index_map,~,~,~,~,~,~,~,~,~,~,~,cell_scores_positive,cell_scores_negative,cell_scores_exclusive,all_clusters_centroids]=find_optimal_clustering(cell_to_index_map,all_to_all_correlation_multi,all_to_all_matrix_multi,max_iterations,cluster_distance_thresh,correlation_thresh,all_centroids_corrected,num_sessions,decision_type,pixel_to_mic,is_figure,results_dir,figures_dir);
else
    errordlg('Registration with a simple threshold is only enabled with spatial correlations')
end

cell_registered_struct=struct;
cell_registered_struct.optimal_cell_to_index_map=optimal_cell_to_index_map;

if get(handles.use_model,'Value')==1;
    cell_registered_struct.cell_scores=register_scores';
    cell_registered_struct.cell_true_positive_scores=cell_scores_positive';
    cell_registered_struct.cell_true_negative_scores=cell_scores_negative';
    cell_registered_struct.cell_exclusivity_scores=cell_scores_exclusive';
end

num_cells=size(optimal_cell_to_index_map,1);
if isfield(data_struct,'overlapping_matrix')
    overlapping_matrix=data_struct.overlapping_matrix;
    is_in_overlapping_FOV=zeros(1,num_cells);
    for n=1:num_cells
        if round(all_clusters_centroids(2,n))>0 && round(all_clusters_centroids(2,n))<=size(overlapping_matrix,1) && round(all_clusters_centroids(1,n))>0 && round(all_clusters_centroids(1,n))<=size(overlapping_matrix,2);
            if overlapping_matrix(round(all_clusters_centroids(2,n)),round(all_clusters_centroids(1,n)))>0
                is_in_overlapping_FOV(n)=1;
            end
        end
    end
    cell_registered_struct.is_cell_in_overlapping_FOV=is_in_overlapping_FOV';
end
cell_registered_struct.registered_cells_centroids=all_clusters_centroids';
cell_registered_struct.all_centroids_corrected=all_centroids_corrected';
cell_registered_struct.all_filters_corrected=all_filters_corrected';

save(['cellRegistered_Final_' datestr(clock,'yyyymmdd_HHMMss') '.mat'],'cell_registered_struct','-v7.3')

reference_session=str2num(get(handles.reference_session,'string'));
if isempty(reference_session) || reference_session<1 || reference_session>num_sessions
    msgbox('Please insert a valid number of sessions/folders required for pre-processing and press enter ');
    waitfor(handles.reference_session,'value',1);
    reference_session=str2double(get(handles.reference_session,'string'));
end
data_struct.reference_session=reference_session;

% Plotting the registration results with the cell maps from all sessions:
filter_thresh=0.2;
all_projections_partial=cell(1,num_sessions);
mutual_projections_partial=cell(1,num_sessions);
cells_in_all_days=find(sum(optimal_cell_to_index_map'>0)==num_sessions);
other_cells=cell(1,num_sessions);
for n=1:num_sessions
    logical_1=sum(optimal_cell_to_index_map'>0)<num_sessions;
    other_cells{n}=find(optimal_cell_to_index_map(:,n)'>0 & logical_1);
end

for n=1:num_sessions
    this_session_filters=all_filters_corrected{n};
    num_filters=size(this_session_filters,1);
    normalized_filters=zeros(size(this_session_filters));
    for k=1:num_filters
        this_filter=this_session_filters(k,:,:);
        this_filter(this_filter<filter_thresh*max(max(this_filter)))=0;
        if max(max(this_filter))>0
            normalized_filters(k,:,:)=this_filter/max(max(this_filter));
        end
    end
    all_projections_partial{n}=zeros(size(this_filter,2),size(this_filter,3),3);
    mutual_projections_partial{n}=zeros(size(this_filter,2),size(this_filter,3),3);
    all_projections_partial{n}(:,:,2)=squeeze(sum(normalized_filters(optimal_cell_to_index_map(cells_in_all_days,n),:,:),1));
    all_projections_partial{n}(:,:,1)=squeeze(sum(normalized_filters(optimal_cell_to_index_map(other_cells{n},n),:,:),1));
    all_projections_partial{n}(:,:,2)=squeeze(sum(normalized_filters(optimal_cell_to_index_map(other_cells{n},n),:,:),1))+squeeze(sum(normalized_filters(optimal_cell_to_index_map(cells_in_all_days,n),:,:),1));
    all_projections_partial{n}(:,:,3)=squeeze(sum(normalized_filters(optimal_cell_to_index_map(other_cells{n},n),:,:),1));
    mutual_projections_partial{n}(:,:,2)=squeeze(sum(normalized_filters(optimal_cell_to_index_map(cells_in_all_days,n),:,:),1));
    all_projections_partial{n}(all_projections_partial{n}>1)=1;
end

subx=4;
suby=ceil(num_sessions/subx);
if num_sessions>4
    figure('units','normalized','outerposition',[0.1 0.04 0.8 0.96])
    for n=1:num_sessions
        subplot(suby,subx,n)
        imagesc(all_projections_partial{n})
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        colormap('gray')
        this_session_name=data_struct.sessions_list{n+1};
        title(strrep(this_session_name,'_','\_'),'fontsize',14,'fontweight','bold')
        if n==1
            text(0.01*size(all_projections_partial{n},1),0.02*size(all_projections_partial{n},2),'Detected in','fontsize',16,'color','g','fontweight','bold')
            text(0.01*size(all_projections_partial{n},1),0.06*size(all_projections_partial{n},2),'all sessions','fontsize',16,'color','g','fontweight','bold')
        end
    end
else
    figure('units','normalized','outerposition',[0.1 0.2 0.9 0.5])
    for n=1:num_sessions
        subplot(1,num_sessions,n)
        imagesc(all_projections_partial{n})
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        colormap('gray')
        this_session_name=data_struct.sessions_list{n+1};
        title(strrep(this_session_name,'_','\_'),'fontsize',14,'fontweight','bold')
        if n==1
            text(0.01*size(all_projections_partial{n},1),0.02*size(all_projections_partial{n},2),'Detected in','fontsize',16,'color','g','fontweight','bold')
            text(0.01*size(all_projections_partial{n},1),0.06*size(all_projections_partial{n},2),'all sessions','fontsize',16,'color','g','fontweight','bold')
        end
    end
end

cd(figures_dir);    
set(gcf,'PaperPositionMode','auto')
savefig('Stage 5 - projcetions - final registration')
saveas(gcf,'Stage 5 - projcetions - final registration','tif')
cd(results_dir);

data_struct.optimal_cell_to_index_map=optimal_cell_to_index_map;

% Saving a log file with all the chosen parameters:
cd(results_dir)
logFile=fopen(['logFile - final registration_' datestr(clock,'yyyymmdd_HHMMss') '.txt'], 'wt' );

% General data parameters:
registered_sessions=data_struct.sessions_list;
for n=1:num_sessions+1
    fprintf(logFile,[registered_sessions{n} '\n']);
    if n==1
        fprintf(logFile,'-----------------\n');
    end
end
fprintf(logFile,'\n\nGeneral data parameters:\n------------------------\n');
if strcmp(data_struct.imaging_technique,'one_photon')
    fprintf(logFile,'Imaging technique - 1-photon\n');
else
    fprintf(logFile,'Imaging technique - 2-photon\n');
end
fprintf(logFile,['Number of sessions - ' num2str(num_sessions) '\n']);
FOV_size=[data_struct.x_scale data_struct.y_scale];
fprintf(logFile,['FOV size - ' num2str(FOV_size(1)) 'X' num2str(FOV_size(2)) ' [mic] ; (x,y)\n']);
image_size=[data_struct.max_x data_struct.max_y];
fprintf(logFile,['Image size - ' num2str(image_size(1)) 'X' num2str(image_size(2)) ' [pixels] ; (x,y)\n\n']);

% Session registration parameters:
fprintf(logFile,'\nSession transformation parameters:\n----------------------------------\n');
reference_session=data_struct.reference_session;
fprintf(logFile,['Reference session - ' num2str(reference_session) '\n']);
if get(handles.translations,'Value')==1
    registration_type='Translation';
else
    registration_type='Translation and Rotation';
end
fprintf(logFile,['Registration Type - ' registration_type '\n\n']);

% Probabilistic model parameters:
fprintf(logFile,'\nCell registration parameters:\n-----------------------------\n');
model_maximal_distance=str2double(get(handles.model_maximal_distance,'String'));
fprintf(logFile,['Model maximal distance - ' num2str(model_maximal_distance) ' [mic]\n']);
num_of_bins=data_struct.number_of_bins;
fprintf(logFile,['Number of bins - ' num2str(num_of_bins) '\n']);

% Initial alignment parameters:
correlation_thresh=str2double(get(handles.correlation_threshold,'string'));
distance_thresh=str2double(get(handles.distance_threshold,'string'));
if get(handles.spatial_correlations,'Value')==1
    initial_type='Spatial correlations';
    fprintf(logFile,['initial registration type - ' initial_type '\n']);
    fprintf(logFile,['Correlation threshold - ' num2str(correlation_thresh) '\n']);
else
    initial_type='Centroid distances';
    fprintf(logFile,['initial registration type - ' initial_type '\n']);
    fprintf(logFile,['Distance threshold - ' num2str(distance_thresh) ' [mic]\n']);
end

% Final alignment parameters:
num_bins_p_same=length(data_struct.true_merge_dist);
if get(handles.use_model,'Value')==1;
    decision_thresh=str2double(get(handles.decision_thresh,'string'));
    fprintf(logFile,'Registration approach - Probabilistic modeling\n');
    if get(handles.spatial_correlations_2,'Value')==1
        final_type='Spatial correlations';
        uncertain_pairs_fraction=round(100*data_struct.uncertain_corr)/100;
        cumsum_false_merge_corr=cumsum(data_struct.false_merge_corr);
        cumsum_false_split_corr=1-cumsum(data_struct.true_merge_corr);
        false_positives=round(100*cumsum_false_merge_corr(round(num_bins_p_same*(1-decision_thresh))))/100;
        false_negatives=round(100*cumsum_false_split_corr(round(num_bins_p_same*(1-decision_thresh))))/100;
    elseif get(handles.p_value,'Value')==1
        final_type='Joint';
        uncertain_pairs_fraction=round(100*data_struct.uncertain)/100;
        cumsum_false_merge=cumsum(data_struct.false_merge);
        cumsum_false_split=1-cumsum(data_struct.true_merge);
        false_positives=round(100*cumsum_false_merge(round(num_bins_p_same*(1-decision_thresh))))/100;
        false_negatives=round(100*cumsum_false_split(round(num_bins_p_same*(1-decision_thresh))))/100;
    elseif get(handles.centroid_distances_2,'Value')==1
        final_type='Centroids distances';
        uncertain_pairs_fraction=round(100*data_struct.uncertain_dist)/100;
        cumsum_false_merge_dist=cumsum(data_struct.false_merge_dist);
        cumsum_false_split_dist=1-cumsum(data_struct.true_merge_dist);
        false_positives=round(100*cumsum_false_merge_dist(round(num_bins_p_same*(1-decision_thresh))))/100;
        false_negatives=round(100*cumsum_false_split_dist(round(num_bins_p_same*(1-decision_thresh))))/100;
    end
    if  strcmp(final_type,'2 dimensions');
        fprintf(logFile,['Final registration type - Joint' '\n']);
    else
        fprintf(logFile,['Final registration type - ' final_type '\n']);
    end    
    fprintf(logFile,['P_same threshold - ' num2str(decision_thresh) ' \n\n']);
else
    fprintf(logFile,'Registration approach - Simple correlation threshold\n');
    simple_thresh=str2double(get(handles.simple_threshold,'string'));
    fprintf(logFile,['Correlation threshold - ' num2str(simple_thresh) ' \n\n']);
    uncertain_pairs_fraction=round(100*data_struct.uncertain_corr)/100;
    cumsum_false_merge_corr=cumsum(data_struct.false_merge_corr);
    cumsum_false_split_corr=1-cumsum(data_struct.true_merge_corr);
    [~,ind_corresponding_p]=min(abs(simple_thresh-data_struct.ctrs{2}));
    p_value_corresponding_to_corr_thresh=1-data_struct.p_value_corr(ind_corresponding_p);
    false_positives=round(100*cumsum_false_merge_corr(round(num_bins_p_same*(1-p_value_corresponding_to_corr_thresh))))/100;
    false_negatives=round(100*cumsum_false_split_corr(round(num_bins_p_same*(1-p_value_corresponding_to_corr_thresh))))/100;
end

% General results:
fprintf(logFile,'\nGeneral results:\n----------------\n');
fprintf(logFile,['Final number of cells - ' num2str(num_cells) '\n']);
average_fraction_active_cells=round(100*sum(sum(optimal_cell_to_index_map>0))/num_sessions/num_cells)/100;
fprintf(logFile,['Average fraction of active cells - ' num2str(average_fraction_active_cells) '\n\n']);

fprintf(logFile,['Fraction of false positives - ' num2str(false_positives) '\n']);
fprintf(logFile,['Fraction of false negatives - ' num2str(false_negatives) '\n']);
fprintf(logFile,['Fraction of uncertain cell-pairs - ' num2str(uncertain_pairs_fraction) '\n']);

if get(handles.use_model,'Value')==1;
    average_score=round(100*mean(register_scores))/100;
    fprintf(logFile,['\nAverage register score - ' num2str(average_score) '\n']);
    average_score_positive=round(100*mean(cell_scores_positive(cell_scores_positive>=0)))/100;
    fprintf(logFile,['Average true positive score - ' num2str(average_score_positive) '\n']);
    average_score_negative=round(100*mean(cell_scores_negative(cell_scores_negative>=0)))/100;
    fprintf(logFile,['Average true negative score - ' num2str(average_score_negative) '\n']);
    average_score_exclusive=round(100*mean(cell_scores_exclusive(cell_scores_exclusive>=0)))/100;
    fprintf(logFile,['Average exclusivity score - ' num2str(average_score_exclusive) '\n']);
    if get(handles.spatial_correlations_2,'Value')==1
        fprintf(logFile,['\nDiscrepancy of the model - ' num2str(round(100*data_struct.MSE_corr_model)/100) '\n']);
    elseif get(handles.centroid_distances_2,'Value')==1
        fprintf(logFile,['Discrepancy of the model - ' num2str(round(100*data_struct.MSE_dist_model)/100) '\n']);
    else
        fprintf(logFile,['Discrepancy of the model - ' num2str(round(100*data_struct.MSE_2d_model)/100) '\n\n']);
    end
end

comments=get(handles.comments,'string');
fprintf(logFile,'\n\nComments:\n---------\n');
num_rows=size(comments,1);
for n=1:num_rows
    fprintf(logFile,comments(n,:));
    fprintf(logFile,'\n');
end
fclose(logFile);

handles.data_struct=data_struct;
guidata(hObject, handles)
msgbox(['Finished performing final cell registration - ' num2str(size(optimal_cell_to_index_map,1)) ' were found'])


% --- Executes on button press in reset.
function reset_Callback(hObject,~, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);
data_struct=handles.data_struct;
data_struct=struct;
data_struct.sessions_list=[];

% Defining the main path of the cell registration scripts:
scriptName=mfilename('fullpath');
[currentpath, ~, ~]=fileparts(scriptName);
cd(currentpath);
addpath(currentpath);

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

sessions_list=cell(1,1);
sessions_list{1}='List of sessions:';
set(handles.list_of_sessions,'value',1)
set(handles.list_of_sessions,'string',sessions_list);

set(handles.red_session,'string','1')
set(handles.green_session,'string','2')
set(handles.blue_session,'string','3')
set(handles.decision_thresh,'string','0.5')
set(handles.initial_p_same_slider,'value',0.5);
set(handles.initial_p_same_threshold,'string','0.5');
set(handles.final_p_same_slider,'value',0.5);
set(handles.model_maximal_distance,'string','12')
set(handles.distance_threshold,'string','5')
set(handles.correlation_threshold,'string','0.65')
set(handles.one_photon,'Value',1);
set(handles.translations_rotations,'Value',1);
set(handles.spatial_correlations_2,'Value',1);
set(handles.use_model,'Value',1);
set(handles.use_joint_model,'Value',0);
set(handles.spatial_correlations,'Value',1);
set(handles.frame_rate,'string',[])
set(handles.frame_rate,'value',0);
set(handles.x_scale,'value',0)
set(handles.y_scale,'value',0)
set(handles.x_scale,'string',[])
set(handles.y_scale,'string',[])
set(handles.pixel_to_mic,'value',0)
set(handles.pixel_to_mic,'string',[])
set(handles.reference_session,'string','1')
set(handles.maximal_rotation','string','30')
set(handles.maximal_rotation','enable','on')
set(handles.distance_threshold,'enable','off')
set(handles.correlation_threshold,'enable','on')
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


function x_scale_Callback(~,~, handles)
% hObject    handle to x_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_scale as text
%        str2double(get(hObject,'String')) returns contents of x_scale as a double
x_scale=get(handles.x_scale,'string');
if ~isempty(x_scale)
    set(handles.x_scale,'value',1)
end

% --- Executes during object creation, after setting all properties.
function x_scale_CreateFcn(hObject,~,~)
% hObject    handle to x_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function y_scale_Callback(~,~, handles)
% hObject    handle to y_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_scale as text
%        str2double(get(hObject,'String')) returns contents of y_scale as a double
y_scale=get(handles.y_scale,'string');
if ~isempty(y_scale)
    set(handles.y_scale,'value',1)
end

% --- Executes during object creation, after setting all properties.
function y_scale_CreateFcn(hObject,~,~)
% hObject    handle to y_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pixel_to_mic_Callback(~,~, handles)
% hObject    handle to pixel_to_mic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pixel_to_mic as text
%        str2double(get(hObject,'String')) returns contents of pixel_to_mic as a double
pixel_to_mic=get(handles.pixel_to_mic,'string');
if ~isempty(pixel_to_mic)
    set(handles.pixel_to_mic,'value',1)
end


% --- Executes during object creation, after setting all properties.
function pixel_to_mic_CreateFcn(hObject,~,~)
% hObject    handle to pixel_to_mic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function reference_session_Callback(~,~, handles)
% hObject    handle to reference_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of reference_session as text
%        str2double(get(hObject,'String')) returns contents of reference_session as a double

reference_session=get(handles.reference_session,'string');
if ~isempty(reference_session)
    set(handles.reference_session,'value',1)
end

% --- Executes during object creation, after setting all properties.
function reference_session_CreateFcn(hObject, ~, ~)
% hObject    handle to reference_session (see GCBO)
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
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in list_of_sessions.
function list_of_sessions_Callback(hObject, eventdata, handles)
% hObject    handle to list_of_sessions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list_of_sessions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_of_sessions


function model_maximal_distance_Callback(hObject, eventdata, handles)
% hObject    handle to model_maximal_distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of model_maximal_distance as text
%        str2double(get(hObject,'String')) returns contents of model_maximal_distance as a double


function number_of_bins_Callback(hObject, eventdata, handles)
% hObject    handle to number_of_bins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_of_bins as text
%        str2double(get(hObject,'String')) returns contents of number_of_bins as a double


% --- Executes when selected object is changed in transformation_type.
function transformation_type_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in transformation_type
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

if get(handles.translations,'Value')==1;
    set(handles.maximal_rotation','enable','off')
else
    set(handles.maximal_rotation','enable','on')
    set(handles.maximal_rotation','string','30')
end


function maximal_rotation_Callback(hObject, eventdata, handles)
% hObject    handle to maximal_rotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maximal_rotation as text
%        str2double(get(hObject,'String')) returns contents of maximal_rotation as a double


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


function frame_rate_Callback(hObject, eventdata, handles)
% hObject    handle to frame_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_rate as text
%        str2double(get(hObject,'String')) returns contents of frame_rate as a double
frame_rate=get(handles.frame_rate,'string');
if ~isempty(frame_rate)
    set(handles.frame_rate,'value',1)
end


% --- Executes during object creation, after setting all properties.
function frame_rate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_rate (see GCBO)
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
% handles    structure with handles and user data (see GUIDATA)

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
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of comments as text
%        str2double(get(hObject,'String')) returns contents of comments as a double


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
% handles    structure with handles and user data (see GUIDATA)

data_struct=handles.data_struct;

if ~isfield(data_struct,'all_projections_corrected')
    errordlg('RGB overlay cannot be displayed before transformation is performed')
end
all_projections_corrected=data_struct.all_projections_corrected;
num_sessions=data_struct.num_sessions;
red_session=str2num(get(handles.red_session,'string'));
green_session=str2num(get(handles.green_session,'string'));
blue_session=str2num(get(handles.blue_session,'string'));
max_y=data_struct.max_y;
max_x=data_struct.max_x;

all_projections_corrected_rgb=zeros(max_y,max_x,3);
all_projections_corrected_rgb(:,:,1)=all_projections_corrected{red_session};
all_projections_corrected_rgb(:,:,2)=all_projections_corrected{green_session};
if num_sessions>2
    all_projections_corrected_rgb(:,:,3)=all_projections_corrected{blue_session};
end

if isfield(data_struct,'overlapping_matrix')
    overlapping_matrix=data_struct.overlapping_matrix;
    all_projections_corrected_rgb=all_projections_corrected_rgb+0.25*repmat(overlapping_matrix,1,1,3);
end
all_projections_corrected_rgb(all_projections_corrected_rgb>1)=1;

axes(handles.axes1);
imagesc(all_projections_corrected_rgb)
title('RGB overlay','FontWeight','Bold','fontsize',16)
set(gca,'xtick',[])
set(gca,'ytick',[])

figure
imshow(all_projections_corrected_rgb)
title_string=['RGB overlay - sessions ' num2str(red_session) ', ' num2str(green_session) ', and ' num2str(blue_session)];
if num_sessions>3
    title(title_string,'FontWeight','Bold','fontsize',18)
else
    title('RGB overlay: Post-transformation','FontWeight','Bold','fontsize',18)
end
handles.data_struct=data_struct;
guidata(hObject, handles)


function red_session_Callback(hObject, eventdata, handles)
% hObject    handle to red_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of red_session as text
%        str2double(get(hObject,'String')) returns contents of red_session as a double


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
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of green_session as text
%        str2double(get(hObject,'String')) returns contents of green_session as a double


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
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of blue_session as text
%        str2double(get(hObject,'String')) returns contents of blue_session as a double


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
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

data_struct=handles.data_struct;
initial_p_same_slider_value=round(100*get(handles.initial_p_same_slider,'value'))/100;
set(handles.initial_p_same_threshold,'string',num2str(initial_p_same_slider_value));
ctrs=data_struct.ctrs;
pixel_to_mic=data_struct.pixel_to_mic;
n_corr=data_struct.n_corr;
normalized_distance=data_struct.normalized_distance;
temp_dist_thresh=str2double(get(handles.model_maximal_distance,'string'))/pixel_to_mic;

if get(handles.spatial_correlations,'Value')==1 
    p_value_corr=data_struct.p_value_corr;
    [~,p_same_ind]=min(abs(initial_p_same_slider_value-(1-p_value_corr)));
else
    p_value_dist=data_struct.p_value_dist;
    [~,p_same_ind]=min(abs(initial_p_same_slider_value-(1-p_value_dist)));
end

if get(handles.spatial_correlations,'Value')==1 
    corr_thresh=round(100*ctrs{2}(p_same_ind))/100;
    set(handles.correlation_threshold,'string',num2str(corr_thresh));
else
    dist_thresh=round(pixel_to_mic*100*ctrs{1}(p_same_ind))/100;
    set(handles.distance_threshold,'string',num2str(dist_thresh));
end

if get(handles.spatial_correlations,'Value')==1
    axes(handles.axes4)
    plot([],[])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    cla(handles.axes4,'reset')
    start_x=0;
    end_x=1;
    y_vec=repmat(n_corr,[2 1]);
    y_vec=y_vec(:);
    x_vec=(ctrs{2}(2:end)+ctrs{2}(1:end-1))/2;
    x_vec=repmat(x_vec,[2 1]);
    x_vec=[start_x; x_vec(:); end_x];
    for run_bins=1:length(x_vec)/2
        current_color=p_value_corr(run_bins)*[1 1 1];
        patch(x_vec([1 1 2 2]+2*(run_bins-1)),[0 [1 1]*n_corr(run_bins) 0],current_color,'EdgeColor',current_color)
        hold on
    end
    plot(x_vec,y_vec,'k-','linewidth',2)
    xlabel('Spatial correlation','fontsize',12,'fontweight','bold')
    x_label=linspace(0,1,6);
    x=linspace(0,1,6);
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
    xlim([0 1])
    [~,ind_005]=min(abs(0.05-(1-p_value_corr)));
    p_005=ctrs{2}(ind_005);
    [~,ind_05]=min(abs(0.5-(1-p_value_corr)));
    p_05=ctrs{2}(ind_05);
    [~,ind_095]=min(abs(0.95-(1-p_value_corr)));
    p_095=ctrs{2}(ind_095);
    hold on
    plot([p_005 p_005],[0 max(n_corr)],'--','linewidth',3,'color','k')
    hold on
    plot([p_05 p_05],[0 max(n_corr)],'--','linewidth',3,'color','k')
    hold on
    plot([p_095 p_095],[0 max(n_corr)],'--','linewidth',3,'color','k')
    hold on
    plot([corr_thresh corr_thresh],[0 max(n_corr)],'linewidth',2,'color','r')
    text(p_005,1.1*max(n_corr),'0.05','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
    text(p_095,1.1*max(n_corr),'0.95','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
    text(p_05,1.1*max(n_corr),'0.5','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')    
else
    axes(handles.axes3)
    plot([],[])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    cla(handles.axes3,'reset')
    start_x=0;
    end_x=temp_dist_thresh*pixel_to_mic;
    y_vec=repmat(normalized_distance,[2 1]);
    y_vec=y_vec(:);
    x_vec=(pixel_to_mic*ctrs{1}(2:end)+pixel_to_mic*ctrs{1}(1:end-1))/2;
    x_vec=repmat(x_vec,[2 1]);
    x_vec=[start_x; x_vec(:); end_x];
    for run_bins=1:length(x_vec)/2
        current_color=p_value_dist(run_bins)*[1 1 1];
        patch(x_vec([1 1 2 2]+2*(run_bins-1)),[0 [1 1]*normalized_distance(run_bins) 0],current_color,'EdgeColor',current_color)
        hold on
    end
    plot(x_vec,y_vec,'k-','linewidth',2)
    xlim([0 pixel_to_mic*temp_dist_thresh])
    xlabel('Centroids distance (\mum)','FontWeight','Bold','fontsize',12)
    x_label=0:3:pixel_to_mic*temp_dist_thresh;
    x=0:3:pixel_to_mic*temp_dist_thresh;
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
    [~,ind_005]=min(abs(0.05-(1-p_value_dist)));
    p_005=pixel_to_mic*ctrs{1}(ind_005);
    [~,ind_05]=min(abs(0.5-(1-p_value_dist)));
    p_05=pixel_to_mic*ctrs{1}(ind_05);
    [~,ind_095]=min(abs(0.95-(1-p_value_dist)));
    p_095=pixel_to_mic*ctrs{1}(ind_095);
    hold on
    plot([p_005 p_005],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
    hold on
    plot([p_05 p_05],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
    hold on
    plot([p_095 p_095],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
    hold on
    plot([dist_thresh dist_thresh],[0 max(normalized_distance)],'linewidth',2,'color','r')
    text(p_005,1.1*max(normalized_distance),'0.05','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
    text(p_095,1.1*max(normalized_distance),'0.95','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
    text(p_05,1.1*max(normalized_distance),'0.5','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
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
% hObject    handle to initial_p_same_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of initial_p_same_threshold as text
%        str2double(get(hObject,'String')) returns contents of initial_p_same_threshold as a double

data_struct=handles.data_struct;
initial_p_same_threshold_value=str2num(get(handles.initial_p_same_threshold,'string'));
set(handles.initial_p_same_slider,'value',initial_p_same_threshold_value);

ctrs=data_struct.ctrs;
pixel_to_mic=data_struct.pixel_to_mic;
n_corr=data_struct.n_corr;
normalized_distance=data_struct.normalized_distance;
temp_dist_thresh=str2double(get(handles.model_maximal_distance,'string'))/pixel_to_mic;

if get(handles.spatial_correlations,'Value')==1 
    p_value_corr=data_struct.p_value_corr;
    [~,p_same_ind]=min(abs(initial_p_same_slider_value-(1-p_value_corr)));
else
    p_value_dist=data_struct.p_value_dist;
    [~,p_same_ind]=min(abs(initial_p_same_slider_value-(1-p_value_dist)));
end

if get(handles.spatial_correlations,'Value')==1 
    corr_thresh=round(100*ctrs{2}(p_same_ind))/100;
    set(handles.correlation_threshold,'string',num2str(corr_thresh));
else
    dist_thresh=round(pixel_to_mic*100*ctrs{1}(p_same_ind))/100;
    set(handles.distance_threshold,'string',num2str(dist_thresh));
end

if get(handles.spatial_correlations,'Value')==1
    axes(handles.axes4)
    plot([],[])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    cla(handles.axes4,'reset')
    start_x=0;
    end_x=1;
    y_vec=repmat(n_corr,[2 1]);
    y_vec=y_vec(:);
    x_vec=(ctrs{2}(2:end)+ctrs{2}(1:end-1))/2;
    x_vec=repmat(x_vec,[2 1]);
    x_vec=[start_x; x_vec(:); end_x];
    for run_bins=1:length(x_vec)/2
        current_color=p_value_corr(run_bins)*[1 1 1];
        patch(x_vec([1 1 2 2]+2*(run_bins-1)),[0 [1 1]*n_corr(run_bins) 0],current_color,'EdgeColor',current_color)
        hold on
    end
    plot(x_vec,y_vec,'k-','linewidth',2)
    xlabel('Spatial correlation','fontsize',12,'fontweight','bold')
    x_label=linspace(0,1,6);
    x=linspace(0,1,6);
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
    xlim([0 1])
    [~,ind_005]=min(abs(0.05-(1-p_value_corr)));
    p_005=ctrs{2}(ind_005);
    [~,ind_05]=min(abs(0.5-(1-p_value_corr)));
    p_05=ctrs{2}(ind_05);
    [~,ind_095]=min(abs(0.95-(1-p_value_corr)));
    p_095=ctrs{2}(ind_095);
    hold on
    plot([p_005 p_005],[0 max(n_corr)],'--','linewidth',3,'color','k')
    hold on
    plot([p_05 p_05],[0 max(n_corr)],'--','linewidth',3,'color','k')
    hold on
    plot([p_095 p_095],[0 max(n_corr)],'--','linewidth',3,'color','k')
    hold on
    plot([corr_thresh corr_thresh],[0 max(n_corr)],'linewidth',2,'color','r')
    text(p_005,1.1*max(n_corr),'0.05','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
    text(p_095,1.1*max(n_corr),'0.95','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
    text(p_05,1.1*max(n_corr),'0.5','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')    
else
    axes(handles.axes3)
    plot([],[])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    cla(handles.axes3,'reset')
    start_x=0;
    end_x=temp_dist_thresh*pixel_to_mic;
    y_vec=repmat(normalized_distance,[2 1]);
    y_vec=y_vec(:);
    x_vec=(pixel_to_mic*ctrs{1}(2:end)+pixel_to_mic*ctrs{1}(1:end-1))/2;
    x_vec=repmat(x_vec,[2 1]);
    x_vec=[start_x; x_vec(:); end_x];
    for run_bins=1:length(x_vec)/2
        current_color=p_value_dist(run_bins)*[1 1 1];
        patch(x_vec([1 1 2 2]+2*(run_bins-1)),[0 [1 1]*normalized_distance(run_bins) 0],current_color,'EdgeColor',current_color)
        hold on
    end
    plot(x_vec,y_vec,'k-','linewidth',2)
    xlim([0 pixel_to_mic*temp_dist_thresh])
    xlabel('Centroids distance (\mum)','FontWeight','Bold','fontsize',12)
    x_label=0:3:pixel_to_mic*temp_dist_thresh;
    x=0:3:pixel_to_mic*temp_dist_thresh;
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
    [~,ind_005]=min(abs(0.05-(1-p_value_dist)));
    p_005=pixel_to_mic*ctrs{1}(ind_005);
    [~,ind_05]=min(abs(0.5-(1-p_value_dist)));
    p_05=pixel_to_mic*ctrs{1}(ind_05);
    [~,ind_095]=min(abs(0.95-(1-p_value_dist)));
    p_095=pixel_to_mic*ctrs{1}(ind_095);
    hold on
    plot([p_005 p_005],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
    hold on
    plot([p_05 p_05],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
    hold on
    plot([p_095 p_095],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
    hold on
    plot([dist_thresh dist_thresh],[0 max(normalized_distance)],'linewidth',2,'color','r')
    text(p_005,1.1*max(normalized_distance),'0.05','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
    text(p_095,1.1*max(normalized_distance),'0.95','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
    text(p_05,1.1*max(normalized_distance),'0.5','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
end

handles.data_struct=data_struct;
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function initial_p_same_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initial_p_same_threshold (see GCBO)
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
% handles    structure with handles and user data (see GUIDATA)

initial_p_same_slider_value=0.5;
data_struct=handles.data_struct;
set(handles.initial_p_same_slider,'value',initial_p_same_slider_value);
set(handles.initial_p_same_threshold,'string',num2str(initial_p_same_slider_value));

ctrs=data_struct.ctrs;
pixel_to_mic=data_struct.pixel_to_mic;
n_corr=data_struct.n_corr;
normalized_distance=data_struct.normalized_distance;
temp_dist_thresh=str2double(get(handles.model_maximal_distance,'string'))/pixel_to_mic;

p_value_corr=data_struct.p_value_corr;
[~,p_same_ind_corr]=min(abs(initial_p_same_slider_value-(1-p_value_corr)));
p_value_dist=data_struct.p_value_dist;
[~,p_same_ind_dist]=min(abs(initial_p_same_slider_value-(1-p_value_dist)));

corr_thresh=round(100*ctrs{2}(p_same_ind_corr))/100;
set(handles.correlation_threshold,'string',num2str(corr_thresh));
dist_thresh=round(pixel_to_mic*100*ctrs{1}(p_same_ind_dist))/100;
set(handles.distance_threshold,'string',num2str(dist_thresh));

axes(handles.axes4)
plot([],[])
set(gca,'xtick',[])
set(gca,'ytick',[])
cla(handles.axes4,'reset')
start_x=0;
end_x=1;
y_vec=repmat(n_corr,[2 1]);
y_vec=y_vec(:);
x_vec=(ctrs{2}(2:end)+ctrs{2}(1:end-1))/2;
x_vec=repmat(x_vec,[2 1]);
x_vec=[start_x; x_vec(:); end_x];
for run_bins=1:length(x_vec)/2
    current_color=p_value_corr(run_bins)*[1 1 1];
    patch(x_vec([1 1 2 2]+2*(run_bins-1)),[0 [1 1]*n_corr(run_bins) 0],current_color,'EdgeColor',current_color)
    hold on
end
plot(x_vec,y_vec,'k-','linewidth',2)
xlabel('Spatial correlation','fontsize',12,'fontweight','bold')
x_label=linspace(0,1,6);
x=linspace(0,1,6);
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
xlim([0 1])
[~,ind_005]=min(abs(0.05-(1-p_value_corr)));
p_005=ctrs{2}(ind_005);
[~,ind_05]=min(abs(0.5-(1-p_value_corr)));
p_05=ctrs{2}(ind_05);
[~,ind_095]=min(abs(0.95-(1-p_value_corr)));
p_095=ctrs{2}(ind_095);
hold on
plot([p_005 p_005],[0 max(n_corr)],'--','linewidth',3,'color','k')
hold on
plot([p_05 p_05],[0 max(n_corr)],'--','linewidth',3,'color','k')
hold on
plot([p_095 p_095],[0 max(n_corr)],'--','linewidth',3,'color','k')
hold on
plot([corr_thresh corr_thresh],[0 max(n_corr)],'linewidth',2,'color','r')
text(p_005,1.1*max(n_corr),'0.05','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
text(p_095,1.1*max(n_corr),'0.95','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
text(p_05,1.1*max(n_corr),'0.5','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
axes(handles.axes3)
plot([],[])
set(gca,'xtick',[])
set(gca,'ytick',[])
cla(handles.axes3,'reset')
start_x=0;
end_x=temp_dist_thresh*pixel_to_mic;
y_vec=repmat(normalized_distance,[2 1]);
y_vec=y_vec(:);
x_vec=(pixel_to_mic*ctrs{1}(2:end)+pixel_to_mic*ctrs{1}(1:end-1))/2;
x_vec=repmat(x_vec,[2 1]);
x_vec=[start_x; x_vec(:); end_x];
for run_bins=1:length(x_vec)/2
    current_color=p_value_dist(run_bins)*[1 1 1];
    patch(x_vec([1 1 2 2]+2*(run_bins-1)),[0 [1 1]*normalized_distance(run_bins) 0],current_color,'EdgeColor',current_color)
    hold on
end
plot(x_vec,y_vec,'k-','linewidth',2)
xlim([0 pixel_to_mic*temp_dist_thresh])
xlabel('Centroids distance (\mum)','FontWeight','Bold','fontsize',12)
x_label=0:3:pixel_to_mic*temp_dist_thresh;
x=0:3:pixel_to_mic*temp_dist_thresh;
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
[~,ind_005]=min(abs(0.05-(1-p_value_dist)));
p_005=pixel_to_mic*ctrs{1}(ind_005);
[~,ind_05]=min(abs(0.5-(1-p_value_dist)));
p_05=pixel_to_mic*ctrs{1}(ind_05);
[~,ind_095]=min(abs(0.95-(1-p_value_dist)));
p_095=pixel_to_mic*ctrs{1}(ind_095);
hold on
plot([p_005 p_005],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
hold on
plot([p_05 p_05],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
hold on
plot([p_095 p_095],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
hold on
plot([dist_thresh dist_thresh],[0 max(normalized_distance)],'linewidth',2,'color','r')
text(p_005,1.1*max(normalized_distance),'0.05','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
text(p_095,1.1*max(normalized_distance),'0.95','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
text(p_05,1.1*max(normalized_distance),'0.5','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')

handles.data_struct=data_struct;
guidata(hObject, handles)


% --- Executes on slider movement.
function final_p_same_slider_Callback(hObject, eventdata, handles)
% hObject    handle to final_p_same_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

data_struct=handles.data_struct;
final_p_same_slider_value=round(100*get(handles.final_p_same_slider,'value'))/100;
set(handles.decision_thresh,'string',num2str(final_p_same_slider_value));

ctrs=data_struct.ctrs;
pixel_to_mic=data_struct.pixel_to_mic;
n_corr=data_struct.n_corr;
normalized_distance=data_struct.normalized_distance;
temp_dist_thresh=str2double(get(handles.model_maximal_distance,'string'))/pixel_to_mic;

if get(handles.spatial_correlations_2,'Value')==1
    p_value_corr=data_struct.p_value_corr;
    [~,p_same_ind]=min(abs(final_p_same_slider_value-(1-p_value_corr)));
elseif get(handles.centroid_distances_2,'Value')==1
    p_value_dist=data_struct.p_value_dist;
    [~,p_same_ind]=min(abs(final_p_same_slider_value-(1-p_value_dist)));
end

if get(handles.spatial_correlations_2,'Value')==1 
    corr_thresh=round(100*ctrs{2}(p_same_ind))/100;
    set(handles.correlation_threshold,'string',num2str(corr_thresh));
elseif get(handles.centroid_distances_2,'Value')==1
    dist_thresh=round(pixel_to_mic*100*ctrs{1}(p_same_ind))/100;
    set(handles.distance_threshold,'string',num2str(dist_thresh));
end

if get(handles.spatial_correlations_2,'Value')==1
    axes(handles.axes4)
    plot([],[])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    cla(handles.axes4,'reset')
    start_x=0;
    end_x=1;
    y_vec=repmat(n_corr,[2 1]);
    y_vec=y_vec(:);
    x_vec=(ctrs{2}(2:end)+ctrs{2}(1:end-1))/2;
    x_vec=repmat(x_vec,[2 1]);
    x_vec=[start_x; x_vec(:); end_x];
    for run_bins=1:length(x_vec)/2
        current_color=p_value_corr(run_bins)*[1 1 1];
        patch(x_vec([1 1 2 2]+2*(run_bins-1)),[0 [1 1]*n_corr(run_bins) 0],current_color,'EdgeColor',current_color)
        hold on
    end
    plot(x_vec,y_vec,'k-','linewidth',2)
    xlabel('Spatial correlation','fontsize',12,'fontweight','bold')
    x_label=linspace(0,1,6);
    x=linspace(0,1,6);
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
    xlim([0 1])
    [~,ind_005]=min(abs(0.05-(1-p_value_corr)));
    p_005=ctrs{2}(ind_005);
    [~,ind_05]=min(abs(0.5-(1-p_value_corr)));
    p_05=ctrs{2}(ind_05);
    [~,ind_095]=min(abs(0.95-(1-p_value_corr)));
    p_095=ctrs{2}(ind_095);
    hold on
    plot([p_005 p_005],[0 max(n_corr)],'--','linewidth',3,'color','k')
    hold on
    plot([p_05 p_05],[0 max(n_corr)],'--','linewidth',3,'color','k')
    hold on
    plot([p_095 p_095],[0 max(n_corr)],'--','linewidth',3,'color','k')
    hold on
    plot([corr_thresh corr_thresh],[0 max(n_corr)],'linewidth',2,'color','r')
    text(p_005,1.1*max(n_corr),'0.05','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
    text(p_095,1.1*max(n_corr),'0.95','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
    text(p_05,1.1*max(n_corr),'0.5','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')    
elseif get(handles.centroid_distances_2,'Value')==1
    axes(handles.axes3)
    plot([],[])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    cla(handles.axes3,'reset')
    start_x=0;
    end_x=temp_dist_thresh*pixel_to_mic;
    y_vec=repmat(normalized_distance,[2 1]);
    y_vec=y_vec(:);
    x_vec=(pixel_to_mic*ctrs{1}(2:end)+pixel_to_mic*ctrs{1}(1:end-1))/2;
    x_vec=repmat(x_vec,[2 1]);
    x_vec=[start_x; x_vec(:); end_x];
    for run_bins=1:length(x_vec)/2
        current_color=p_value_dist(run_bins)*[1 1 1];
        patch(x_vec([1 1 2 2]+2*(run_bins-1)),[0 [1 1]*normalized_distance(run_bins) 0],current_color,'EdgeColor',current_color)
        hold on
    end
    plot(x_vec,y_vec,'k-','linewidth',2)
    xlim([0 pixel_to_mic*temp_dist_thresh])
    xlabel('Centroids distance (\mum)','FontWeight','Bold','fontsize',12)
    x_label=0:3:pixel_to_mic*temp_dist_thresh;
    x=0:3:pixel_to_mic*temp_dist_thresh;
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
    [~,ind_005]=min(abs(0.05-(1-p_value_dist)));
    p_005=pixel_to_mic*ctrs{1}(ind_005);
    [~,ind_05]=min(abs(0.5-(1-p_value_dist)));
    p_05=pixel_to_mic*ctrs{1}(ind_05);
    [~,ind_095]=min(abs(0.95-(1-p_value_dist)));
    p_095=pixel_to_mic*ctrs{1}(ind_095);
    hold on
    plot([p_005 p_005],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
    hold on
    plot([p_05 p_05],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
    hold on
    plot([p_095 p_095],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
    hold on
    plot([dist_thresh dist_thresh],[0 max(normalized_distance)],'linewidth',2,'color','r')
    text(p_005,1.1*max(normalized_distance),'0.05','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
    text(p_095,1.1*max(normalized_distance),'0.95','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
    text(p_05,1.1*max(normalized_distance),'0.5','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
else
    if isfield(data_struct,'log_grid')
        log_grid=data_struct.log_grid;
        p_value=data_struct.p_value;
        num_bins=data_struct.number_of_bins;
        axes(handles.axes5)
        imagesc(log_grid)
        colormap('jet')
        y=round(linspace(1,num_bins,5));
        y_label=round(10*linspace(pixel_to_mic*temp_dist_thresh,0,5))/10;
        x=round(linspace(1,num_bins,6));
        x_label=linspace(0,1,6);
        set(gca,'YTick',y)
        set(gca,'YTickLabel',y_label,'fontsize',14,'FontWeight','Bold')
        set(gca,'XTick',x)
        set(gca,'XTickLabel',x_label,'fontsize',14,'FontWeight','Bold')
        xlabel('Spatial correlation','FontWeight','Bold','fontsize',14)
        ylabel('Centroids distance (\mum)','FontWeight','Bold','fontsize',14)
        freezeColors
        hold on
        contour_values=0.05:0.05:0.95;
        contour(p_value,contour_values,'linewidth',1)
        colormap('gray')
        hold on
        contour(p_value,[0.05 0.5 0.95],'linewidth',3)
        hold on
        contour(p_value,[1-(final_p_same_slider_value-0.005) 1-(final_p_same_slider_value+0.005)],'linewidth',2,'color','r')
        set(gca,'fontsize',14)
    end
end

handles.data_struct=data_struct;
guidata(hObject, handles)


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
% handles    structure with handles and user data (see GUIDATA)

final_p_same_slider_value=0.5;
data_struct=handles.data_struct;
set(handles.final_p_same_slider,'value',final_p_same_slider_value);
set(handles.decision_thresh,'string',num2str(final_p_same_slider_value));

optimal_threshold=data_struct.thresh_corr_from_intersection;
set(handles.simple_threshold,'string',optimal_threshold);

ctrs=data_struct.ctrs;
pixel_to_mic=data_struct.pixel_to_mic;
n_corr=data_struct.n_corr;
normalized_distance=data_struct.normalized_distance;
log_grid=data_struct.log_grid;
temp_dist_thresh=str2double(get(handles.model_maximal_distance,'string'))/pixel_to_mic;

p_value_corr=data_struct.p_value_corr;
[~,p_same_ind_corr]=min(abs(final_p_same_slider_value-(1-p_value_corr)));
p_value_dist=data_struct.p_value_dist;
[~,p_same_ind_dist]=min(abs(final_p_same_slider_value-(1-p_value_dist)));

corr_thresh=round(100*ctrs{2}(p_same_ind_corr))/100;
set(handles.correlation_threshold,'string',num2str(corr_thresh));
dist_thresh=round(pixel_to_mic*100*ctrs{1}(p_same_ind_dist))/100;
set(handles.distance_threshold,'string',num2str(dist_thresh));

axes(handles.axes4)
plot([],[])
set(gca,'xtick',[])
set(gca,'ytick',[])
cla(handles.axes4,'reset')
start_x=0;
end_x=1;
y_vec=repmat(n_corr,[2 1]);
y_vec=y_vec(:);
x_vec=(ctrs{2}(2:end)+ctrs{2}(1:end-1))/2;
x_vec=repmat(x_vec,[2 1]);
x_vec=[start_x; x_vec(:); end_x];
for run_bins=1:length(x_vec)/2
    current_color=p_value_corr(run_bins)*[1 1 1];
    patch(x_vec([1 1 2 2]+2*(run_bins-1)),[0 [1 1]*n_corr(run_bins) 0],current_color,'EdgeColor',current_color)
    hold on
end
plot(x_vec,y_vec,'k-','linewidth',2)
xlabel('Spatial correlation','fontsize',12,'fontweight','bold')
x_label=linspace(0,1,6);
x=linspace(0,1,6);
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
xlim([0 1])
[~,ind_005]=min(abs(0.05-(1-p_value_corr)));
p_005=ctrs{2}(ind_005);
[~,ind_05]=min(abs(0.5-(1-p_value_corr)));
p_05=ctrs{2}(ind_05);
[~,ind_095]=min(abs(0.95-(1-p_value_corr)));
p_095=ctrs{2}(ind_095);
hold on
plot([p_005 p_005],[0 max(n_corr)],'--','linewidth',3,'color','k')
hold on
plot([p_05 p_05],[0 max(n_corr)],'--','linewidth',3,'color','k')
hold on
plot([p_095 p_095],[0 max(n_corr)],'--','linewidth',3,'color','k')
hold on
plot([corr_thresh corr_thresh],[0 max(n_corr)],'linewidth',2,'color','r')
text(p_005,1.1*max(n_corr),'0.05','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
text(p_095,1.1*max(n_corr),'0.95','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
text(p_05,1.1*max(n_corr),'0.5','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
axes(handles.axes3)
plot([],[])
set(gca,'xtick',[])
set(gca,'ytick',[])
cla(handles.axes3,'reset')
start_x=0;
end_x=temp_dist_thresh*pixel_to_mic;
y_vec=repmat(normalized_distance,[2 1]);
y_vec=y_vec(:);
x_vec=(pixel_to_mic*ctrs{1}(2:end)+pixel_to_mic*ctrs{1}(1:end-1))/2;
x_vec=repmat(x_vec,[2 1]);
x_vec=[start_x; x_vec(:); end_x];
for run_bins=1:length(x_vec)/2
    current_color=p_value_dist(run_bins)*[1 1 1];
    patch(x_vec([1 1 2 2]+2*(run_bins-1)),[0 [1 1]*normalized_distance(run_bins) 0],current_color,'EdgeColor',current_color)
    hold on
end
plot(x_vec,y_vec,'k-','linewidth',2)
xlim([0 pixel_to_mic*temp_dist_thresh])
xlabel('Centroids distance (\mum)','FontWeight','Bold','fontsize',12)
x_label=0:3:pixel_to_mic*temp_dist_thresh;
x=0:3:pixel_to_mic*temp_dist_thresh;
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
[~,ind_005]=min(abs(0.05-(1-p_value_dist)));
p_005=pixel_to_mic*ctrs{1}(ind_005);
[~,ind_05]=min(abs(0.5-(1-p_value_dist)));
p_05=pixel_to_mic*ctrs{1}(ind_05);
[~,ind_095]=min(abs(0.95-(1-p_value_dist)));
p_095=pixel_to_mic*ctrs{1}(ind_095);
hold on
plot([p_005 p_005],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
hold on
plot([p_05 p_05],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
hold on
plot([p_095 p_095],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
hold on
plot([dist_thresh dist_thresh],[0 max(normalized_distance)],'linewidth',2,'color','r')
text(p_005,1.1*max(normalized_distance),'0.05','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
text(p_095,1.1*max(normalized_distance),'0.95','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
text(p_05,1.1*max(normalized_distance),'0.5','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')

if isfield(data_struct,'log_grid')
    log_grid=data_struct.log_grid;
    p_value=data_struct.p_value;
    num_bins=data_struct.number_of_bins;
    axes(handles.axes5)
    imagesc(log_grid)
    colormap('jet')
    y=round(linspace(1,num_bins,5));
    y_label=round(10*linspace(pixel_to_mic*temp_dist_thresh,0,5))/10;
    x=round(linspace(1,num_bins,6));
    x_label=linspace(0,1,6);
    set(gca,'YTick',y)
    set(gca,'YTickLabel',y_label,'fontsize',14,'FontWeight','Bold')
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',14,'FontWeight','Bold')
    xlabel('Spatial correlation','FontWeight','Bold','fontsize',14)
    ylabel('Centroids distance (\mum)','FontWeight','Bold','fontsize',14)
    freezeColors
    hold on
    contour_values=0.05:0.05:0.95;
    contour(p_value,contour_values,'linewidth',1)
    colormap('gray')
    hold on
    contour(p_value,[0.05 0.5 0.95],'linewidth',3)
    hold on
    contour(p_value,[0.495 0.505],'linewidth',2,'color','r')
    set(gca,'fontsize',14)
end

handles.data_struct=data_struct;
guidata(hObject, handles)

function decision_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to decision_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of decision_thresh as text
%        str2double(get(hObject,'String')) returns contents of decision_thresh as a double

data_struct=handles.data_struct;
final_p_same_slider_value=str2num(get(handles.decision_thresh,'string'));
set(handles.final_p_same_slider,'value',final_p_same_slider_value);

ctrs=data_struct.ctrs;
pixel_to_mic=data_struct.pixel_to_mic;
n_corr=data_struct.n_corr;
normalized_distance=data_struct.normalized_distance;
temp_dist_thresh=str2double(get(handles.model_maximal_distance,'string'))/pixel_to_mic;

if get(handles.spatial_correlations_2,'Value')==1
    p_value_corr=data_struct.p_value_corr;
    [~,p_same_ind]=min(abs(final_p_same_slider_value-(1-p_value_corr)));
elseif get(handles.centroid_distances_2,'Value')==1
    p_value_dist=data_struct.p_value_dist;
    [~,p_same_ind]=min(abs(final_p_same_slider_value-(1-p_value_dist)));
end

if get(handles.spatial_correlations_2,'Value')==1 
    corr_thresh=round(100*ctrs{2}(p_same_ind))/100;
    set(handles.correlation_threshold,'string',num2str(corr_thresh));
elseif get(handles.centroid_distances_2,'Value')==1
    dist_thresh=round(pixel_to_mic*100*ctrs{1}(p_same_ind))/100;
    set(handles.distance_threshold,'string',num2str(dist_thresh));
end

if get(handles.spatial_correlations_2,'Value')==1
    axes(handles.axes4)
    plot([],[])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    cla(handles.axes4,'reset')
    start_x=0;
    end_x=1;
    y_vec=repmat(n_corr,[2 1]);
    y_vec=y_vec(:);
    x_vec=(ctrs{2}(2:end)+ctrs{2}(1:end-1))/2;
    x_vec=repmat(x_vec,[2 1]);
    x_vec=[start_x; x_vec(:); end_x];
    for run_bins=1:length(x_vec)/2
        current_color=p_value_corr(run_bins)*[1 1 1];
        patch(x_vec([1 1 2 2]+2*(run_bins-1)),[0 [1 1]*n_corr(run_bins) 0],current_color,'EdgeColor',current_color)
        hold on
    end
    plot(x_vec,y_vec,'k-','linewidth',2)
    xlabel('Spatial correlation','fontsize',12,'fontweight','bold')
    x_label=linspace(0,1,6);
    x=linspace(0,1,6);
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
    xlim([0 1])
    [~,ind_005]=min(abs(0.05-(1-p_value_corr)));
    p_005=ctrs{2}(ind_005);
    [~,ind_05]=min(abs(0.5-(1-p_value_corr)));
    p_05=ctrs{2}(ind_05);
    [~,ind_095]=min(abs(0.95-(1-p_value_corr)));
    p_095=ctrs{2}(ind_095);
    hold on
    plot([p_005 p_005],[0 max(n_corr)],'--','linewidth',3,'color','k')
    hold on
    plot([p_05 p_05],[0 max(n_corr)],'--','linewidth',3,'color','k')
    hold on
    plot([p_095 p_095],[0 max(n_corr)],'--','linewidth',3,'color','k')
    hold on
    plot([corr_thresh corr_thresh],[0 max(n_corr)],'linewidth',2,'color','r')
    text(p_005,1.1*max(n_corr),'0.05','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
    text(p_095,1.1*max(n_corr),'0.95','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
    text(p_05,1.1*max(n_corr),'0.5','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')    
elseif get(handles.centroid_distances_2,'Value')==1
    axes(handles.axes3)
    plot([],[])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    cla(handles.axes3,'reset')
    start_x=0;
    end_x=temp_dist_thresh*pixel_to_mic;
    y_vec=repmat(normalized_distance,[2 1]);
    y_vec=y_vec(:);
    x_vec=(pixel_to_mic*ctrs{1}(2:end)+pixel_to_mic*ctrs{1}(1:end-1))/2;
    x_vec=repmat(x_vec,[2 1]);
    x_vec=[start_x; x_vec(:); end_x];
    for run_bins=1:length(x_vec)/2
        current_color=p_value_dist(run_bins)*[1 1 1];
        patch(x_vec([1 1 2 2]+2*(run_bins-1)),[0 [1 1]*normalized_distance(run_bins) 0],current_color,'EdgeColor',current_color)
        hold on
    end
    plot(x_vec,y_vec,'k-','linewidth',2)
    xlim([0 pixel_to_mic*temp_dist_thresh])
    xlabel('Centroids distance (\mum)','FontWeight','Bold','fontsize',12)
    x_label=0:3:pixel_to_mic*temp_dist_thresh;
    x=0:3:pixel_to_mic*temp_dist_thresh;
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
    [~,ind_005]=min(abs(0.05-(1-p_value_dist)));
    p_005=pixel_to_mic*ctrs{1}(ind_005);
    [~,ind_05]=min(abs(0.5-(1-p_value_dist)));
    p_05=pixel_to_mic*ctrs{1}(ind_05);
    [~,ind_095]=min(abs(0.95-(1-p_value_dist)));
    p_095=pixel_to_mic*ctrs{1}(ind_095);
    hold on
    plot([p_005 p_005],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
    hold on
    plot([p_05 p_05],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
    hold on
    plot([p_095 p_095],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
    hold on
    plot([dist_thresh dist_thresh],[0 max(normalized_distance)],'linewidth',2,'color','r')
    text(p_005,1.1*max(normalized_distance),'0.05','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
    text(p_095,1.1*max(normalized_distance),'0.95','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
    text(p_05,1.1*max(normalized_distance),'0.5','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
else
    if isfield(data_struct,'log_grid')
        log_grid=data_struct.log_grid;
        p_value=data_struct.p_value;
        num_bins=data_struct.number_of_bins;
        axes(handles.axes5)
        imagesc(log_grid)
        colormap('jet')
        y=round(linspace(1,num_bins,5));
        y_label=round(10*linspace(pixel_to_mic*temp_dist_thresh,0,5))/10;
        x=round(linspace(1,num_bins,6));
        x_label=linspace(0,1,6);
        set(gca,'YTick',y)
        set(gca,'YTickLabel',y_label,'fontsize',14,'FontWeight','Bold')
        set(gca,'XTick',x)
        set(gca,'XTickLabel',x_label,'fontsize',14,'FontWeight','Bold')
        xlabel('Spatial correlation','FontWeight','Bold','fontsize',14)
        ylabel('Centroids distance (\mum)','FontWeight','Bold','fontsize',14)
        freezeColors
        hold on
        contour_values=0.05:0.05:0.95;
        contour(p_value,contour_values,'linewidth',1)
        colormap('gray')
        hold on
        contour(p_value,[0.05 0.5 0.95],'linewidth',3)
        hold on
        contour(p_value,[1-(final_p_same_slider_value-0.005) 1-(final_p_same_slider_value+0.005)],'linewidth',2,'color','r')
        set(gca,'fontsize',14)
    end
end

handles.data_struct=data_struct;
guidata(hObject, handles)


% --- Executes on button press in use_simple_threshold.
function use_simple_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to use_simple_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_simple_threshold


% --- Executes on button press in use_model.
function use_model_Callback(hObject, eventdata, handles)
% hObject    handle to use_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_model


% --- Executes during object creation, after setting all properties.
function Registration_approach_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Registration_approach (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes when selected object is changed in Registration_approach.
function Registration_approach_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in transformation_type
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

if get(handles.use_model,'Value')==1;
    set(handles.use_simple_threshold,'value',0)
    set(handles.decision_thresh,'enable','on')
    set(handles.simple_threshold,'enable','off')
elseif get(handles.use_simple_threshold,'Value')==1;
    set(handles.use_model,'value',0)
    set(handles.decision_thresh,'enable','off')
    set(handles.simple_threshold,'enable','on')
end


function simple_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to simple_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of simple_threshold as text
%        str2double(get(hObject,'String')) returns contents of simple_threshold as a double


% --- Executes during object creation, after setting all properties.
function simple_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to simple_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in use_joint_model.
function use_joint_model_Callback(hObject, eventdata, handles)
% hObject    handle to use_joint_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_joint_model

if get(handles.use_joint_model,'value')==1
    errordlg('The joint model is not implemented in this version. Please uncheck this radio button')
end

% --- Executes when selected object is changed in imaging_technique_select.
function imaging_technique_select_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in imaging_technique_select 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.one_photon,'Value')==1;
    set(handles.spatial_correlations,'Value',1);
    set(handles.spatial_correlations_2,'Value',1);
else
    set(handles.use_joint_model,'Value',0);
    set(handles.centroid_distances,'Value',1);
    set(handles.centroid_distances_2,'Value',1);
end


% --- Executes when selected object is changed in probability_model_select.
function probability_model_select_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in probability_model_select 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

