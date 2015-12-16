function varargout = CellAnnotatorGUI(varargin)
% CELLANNOTATORGUI MATLAB code for CellAnnotatorGUI.fig
%      CELLANNOTATORGUI, by itself, creates a new CELLANNOTATORGUI or raises the existing
%      singleton*.
%
%      H = CELLANNOTATORGUI returns the handle to a new CELLANNOTATORGUI or the handle to
%      the existing singleton*.
%
%      CELLANNOTATORGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CELLANNOTATORGUI.M with the given input arguments.
%
%      CELLANNOTATORGUI('Property','Value',...) creates a new CELLANNOTATORGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CellAnnotatorGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CellAnnotatorGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CellAnnotatorGUI

% Last Modified by GUIDE v2.5 10-Dec-2015 17:20:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CellAnnotatorGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @CellAnnotatorGUI_OutputFcn, ...
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

% --- Executes just before CellAnnotatorGUI is made visible.
function CellAnnotatorGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CellAnnotatorGUI (see VARARGIN)

% Choose default command line output for CellAnnotatorGUI
set(hObject,'Toolbar','figure');
handles.output = hObject;

handles.ImageSequence.FileName=[];
% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using CellAnnotatorGUI.
% % if strcmp(get(hObject,'Visible'),'off')
% %     plot(rand(5));
% % end


% UIWAIT makes CellAnnotatorGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CellAnnotatorGUI_OutputFcn(hObject, eventdata, handles)
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
axes(handles.axes1);
cla;

popup_sel_index = get(handles.popupmenu1, 'Value');
switch popup_sel_index
    case 1
        plot(rand(5));
    case 2
        plot(sin(1:0.01:25.99));
    case 3
        bar(1:.5:10);
    case 4
        plot(membrane);
    case 5
        surf(peaks);
end


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('*.tif',[],'X:\raid0\expk\expk0256\0001_threshold_focused\focused_one_test\img_0001.tif');
if ~isequal(FileName, 0)
    handles.ImageSequence.FileName=FileName;
    handles.ImageSequence.PathName=PathName;
    
    handles=ReadImageSequenceParameters(handles);
    
    imagesc(handles.ImageSequence.ImageData{handles.CurrentImageIndex});
    colormap('gray');
    axis equal off;
    
    set(handles.slider1,'Value',handles.CurrentImageIndex,'Min',1,'Max',handles.ImageSequence.NumImages,'SliderStep',[1/(handles.ImageSequence.NumImages-1) 1/10]);
    set(get(handles.axes1,'Children'),'ButtonDownFcn',@position_and_button);
%     set(handles.figure1,'WindowButtonD',[])    
    
    guidata(hObject, handles);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});


%% Model functions
function handles=ReadImageSequenceParameters(handles)

files=dir([handles.ImageSequence.PathName handles.ImageSequence.FileName(1) '*']);

% find files
[~, Name,handles.ImageSequence.FileExtension]=fileparts(files(1).name);

% find correct file name
NotDigits=regexp(Name,'\D');
handles.ImageSequence.NumDigits=length(Name)-NotDigits(end);
handles.ImageSequence.FileName=handles.ImageSequence.FileName(1:NotDigits(end));

% sort files into number order

FileTable=(1:length(files))';
FileTable(:,2)=0;

for FileIndex=1:size(FileTable,1)
    FileTable(FileIndex,2)=str2num(files(FileIndex).name((NotDigits(end)+1):(NotDigits(end)+handles.ImageSequence.NumDigits)));
end
FileTable=sortrows(FileTable,2);

for FileIndex=1:size(FileTable,1)
    handles.ImageSequence.ImageFileNames(FileIndex)=files(FileTable(FileIndex,1));
end

% other image sequence parameters
handles.ImageSequence.NumImages=size(FileTable,1);
handles.CurrentImageIndex=1;

ImageInfo=imfinfo(ImageSequenceGetFilePath(handles,handles.CurrentImageIndex));

handles.ImageSequence.Width=ImageInfo.Width;
handles.ImageSequence.Height=ImageInfo.Height;

handles.ImageSequence.ImageData=cell(handles.ImageSequence.NumImages,1);

for FileIndex=1:handles.ImageSequence.NumImages
    % tell user that images are loaded
    disp(FileIndex);
    handles.ImageSequence.ImageData{FileIndex}=uint8(ImageSequenceGetImage(handles,FileIndex));
end

handles.ImageSequence.ImageDataLoaded=ones(handles.ImageSequence.NumImages,1);




function FileName=ImageSequenceGetFilePath(handles,FileIndex)
if (FileIndex<1) || (FileIndex>length(handles.ImageSequence.ImageFileNames))
    FileName=[];
    return 
else    
    FileName=[handles.ImageSequence.PathName handles.ImageSequence.ImageFileNames(FileIndex).name];
end

function ImageData=ImageSequenceGetImage(handles,FileIndex)
if (FileIndex<1) || (FileIndex>length(handles.ImageSequence.ImageFileNames))
    ImageData=[];
    return 
else    
    success=false;
    while(success==false)
        try
            ImageData=uint8(imread(ImageSequenceGetFilePath(handles,FileIndex)));
            success=true;
        catch
            pause(0.1);
        end
    end
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.CurrentImageIndex=round(get(handles.slider1,'value'));
% check if image data is loaded
if handles.ImageSequence.ImageDataLoaded(handles.CurrentImageIndex)==0
    handles.ImageSequence.ImageData{handles.CurrentImageIndex}=uint8(ImageSequenceGetImage(handles,handles.CurrentImageIndex));
    handles.ImageSequence.ImageDataLoaded(handles.CurrentImageIndex)=1;
end

set(get(handles.axes1,'Children'),'CData',handles.ImageSequence.ImageData{handles.CurrentImageIndex});






% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axis equal

function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

% keyboard


% % switch eventdata.Key
% % % zoom in
% %     case 'equal'
% %         CameraViewAngle=get(handles.axes1,'CameraViewAngle');
% %         CameraViewAngle=CameraViewAngle/2;
% %         set(handles.axes1,'CameraViewAngle',CameraViewAngle);
% % % zoom out
% %     case 'hyphen'
% %         CameraViewAngle=get(handles.axes1,'CameraViewAngle');
% %         CameraViewAngle=CameraViewAngle*2;
% %         set(handles.axes1,'CameraViewAngle',CameraViewAngle);
% % % move to next image
% %     case 'rightarrow'
% %         if (get(handles.slider1,'value')+1)<=get(handles.slider1,'max')
% %             set(handles.slider1,'value',get(handles.slider1,'value')+1);
% %             slider1_Callback(handles.slider1,[],handles);
% %         end
% % % move to previous image
% %     case 'leftarrow'
% %         if (get(handles.slider1,'value')-1)>=get(handles.slider1,'min')
% %             set(handles.slider1,'value',get(handles.slider1,'value')-1);
% %             slider1_Callback(handles.slider1,[],handles);
% %         end
% % end


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


switch eventdata.Key
% zoom in    
    case 'equal'                
        CameraViewAngle=get(handles.axes1,'CameraViewAngle');        
        CameraViewAngle=CameraViewAngle/2;        
        set(handles.axes1,'CameraViewAngle',CameraViewAngle);
% zoom out        
    case 'hyphen'
        CameraViewAngle=get(handles.axes1,'CameraViewAngle');        
        CameraViewAngle=CameraViewAngle*2;        
        set(handles.axes1,'CameraViewAngle',CameraViewAngle);
% move to next image
    case 'rightarrow'
        if (get(handles.slider1,'value')+1)<=get(handles.slider1,'max')
            set(handles.slider1,'value',get(handles.slider1,'value')+1);
            slider1_Callback(handles.slider1,[],handles);
        end
% move to previous image
    case 'leftarrow'
        if (get(handles.slider1,'value')-1)>=get(handles.slider1,'min')
            set(handles.slider1,'value',get(handles.slider1,'value')-1);
            slider1_Callback(handles.slider1,[],handles);
        end
end


function position_and_button(hObject,eventdata)

Position = get( ancestor(hObject,'axes'), 'CurrentPoint' );
Button = get( ancestor(hObject,'figure'), 'SelectionType' );

handles=guidata(gcf);
% left click
if strcmp(Button,'normal')
	disp('left click, select track');
    disp(Position(1,1:2));    
% right click    
elseif strcmp(Button,'alt')
    disp('right click')
    disp('open menu')
    
    RightClickMenu = uicontextmenu;
    
%     hcb1=['disp(''Add new track'')'];
    item1=uimenu(RightClickMenu,'Label','Add new track','Callback',@AddNewTrack);
    set(get(gca,'children'),'uicontextmenu',RightClickMenu)
    
%     ImageHandle=get(handles.axes1,'Children');    
%     set(handles.figure1,'UIContextMenu',hcmenu );
%     m1 = uimenu(handles.figure1,'Label','Add new track','Callback',@AddNewTrack);
%     m1 = uimenu(handles.axes1,'Label','Add new track');
    
end

function AddNewTrack(source,callbackdata)
disp('Add track callback');


    
