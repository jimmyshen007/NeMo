function varargout = import_colormap(varargin)
% IMPORT_COLORMAP MATLAB code for import_colormap.fig
%      IMPORT_COLORMAP, by itself, creates a new IMPORT_COLORMAP or raises the existing
%      singleton*.
%
%      H = IMPORT_COLORMAP returns the handle to a new IMPORT_COLORMAP or the handle to
%      the existing singleton*.
%
%      IMPORT_COLORMAP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMPORT_COLORMAP.M with the given input arguments.
%
%      IMPORT_COLORMAP('Property','Value',...) creates a new IMPORT_COLORMAP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before import_colormap_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to import_colormap_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help import_colormap

% Last Modified by GUIDE v2.5 08-May-2013 00:22:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @import_colormap_OpeningFcn, ...
                   'gui_OutputFcn',  @import_colormap_OutputFcn, ...
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


% --- Executes just before import_colormap is made visible.
function import_colormap_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to import_colormap (see VARARGIN)

% Choose default command line output for import_colormap
handles.output = hObject;
handles.mainHandle = getappdata(0,'mainHandle');
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes import_colormap wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = import_colormap_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.mat'});

if isequal(filename,0)
    disp('User selected Cancel');
elseif size(filename,1) > 1
    disp('Please choose only one file at a time'); 
else
    CMvar = whos('-file',[pathname filesep filename]); %Need error handling for non-symmetric and dim > 2
    matIn = load([pathname filesep filename]);
    eval(['cMap = matIn.' CMvar.name ';']);
    if size(cMap,2) ~= 3
        disp('Incorrect dimensions for RGB colormap');
    else
        setappdata(handles.output,'newColorMap',cMap);
    end
end
    

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
importCM = inputdlg({'Please enter the variable to import:'});
if ~isempty(importCM)
    try
        CM = evalin('base', importCM{1});
    catch
        CM=[];
        disp(['Variable name ''' importCM{1} ''' not found.']);
    end
    if ~isempty(CM) && size(CM,2) == 3
        setappdata(handles.output,'newColorMap',CM);
    elseif size(CM,2) ~= 3 
        disp('Imported variable incorrect dimensions for RGB colormap');
    end
end

% --- Executes on button press in pushbutton3.  Cancel button
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.output);

% --- Executes on button press in pushbutton4.  OK button
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(getappdata(handles.output,'newColorMap'))
    B = guidata(handles.mainHandle);
    B(B(end).currentVol).custom_colormap = getappdata(handles.output,'newColorMap');
    guidata(handles.mainHandle,B);
end
close(handles.output);
