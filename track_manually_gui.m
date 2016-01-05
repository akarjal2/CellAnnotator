function varargout = track_manually_gui(varargin)
% TRACK_MANUALLY_GUI MATLAB code for track_manually_gui.fig
%      TRACK_MANUALLY_GUI, by itself, creates a new TRACK_MANUALLY_GUI or raises the existing
%      singleton*.
%
%      H = TRACK_MANUALLY_GUI returns the handle to a new TRACK_MANUALLY_GUI or the handle to
%      the existing singleton*.
%
%      TRACK_MANUALLY_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRACK_MANUALLY_GUI.M with the given input arguments.
%
%      TRACK_MANUALLY_GUI('Property','Value',...) creates a new TRACK_MANUALLY_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before track_manually_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to track_manually_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help track_manually_gui

% Last Modified by GUIDE v2.5 04-Jan-2016 15:24:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',        mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @track_manually_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @track_manually_gui_OutputFcn, ...
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
end


% --- Executes just before track_manually_gui is made visible.
function track_manually_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to track_manually_gui (see VARARGIN)

% Choose default command line output for track_manually_gui
handles.output = hObject;

set(handles.figure1,'WindowButtonD',[])

guidata(hObject, handles);
end

% UIWAIT makes track_manually_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = track_manually_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global param

% if get(handles.checkbox1,'value')==1
%     img=antti_merge_RGBG(handles.s_imgs(:,:,round(get(handles.slider2,'Value'))-handles.time_interval(1)+1),[],[],handles.o_imgs(:,:,round(get(handles.slider2,'Value'))-handles.time_interval(1)+1));
%     set(handles.cur_img,'CData',img);
% else
%     set(handles.cur_img,'CData',handles.o_imgs(:,:,round(get(handles.slider2,'Value'))-handles.time_interval(1)+1));    
% end

draw_image(handles);

set(handles.text1,'string',num2str(round(get(hObject,'Value'))));

if ~isempty(handles.selection)
    time=round(get(handles.slider2,'Value'))-handles.time_interval(1)+1;

    for i_ind=1:size(handles.selection,1);
        t_inds=param.tracks(handles.selection(i_ind,3)).t==time;
        if sum(t_inds)>0
            handles.selection(i_ind,[1 2 4])=[param.tracks(handles.selection(i_ind,3)).cent(t_inds,1:2) param.tracks(handles.selection(i_ind,3)).t(t_inds)];
        else
            handles.selection(i_ind,[1 2 4])=0;
        end
    end

    handles.centroid_scatter=draw_selection(handles);
    handles.centroid_tracks_lines=draw_tracks(handles);
    handles.centroid_neighbour_lines=draw_neighbours(handles);
    handles.centroid_polygons=draw_polygons(handles);
        
    guidata(hObject, handles);
end

if (round(get(handles.slider2,'Value'))-handles.time_interval(1)+1)==handles.time & get(handles.checkbox7,'value')==1
    set(handles.projection_scatter,'visible','on');
else
    set(handles.projection_scatter,'visible','off');
end

if get(handles.checkbox9,'value')==1
%     keyboard
    movegui(handles.figure1);
%     keyboard
    aa=getframe(handles.axis1);
%     keyboard
%     aa=aa.cdata;
%     aa=aa.cdata(238:512,1:1101,:);
%     aa=aa.cdata(1:749,30:1015,:);
    aa=aa.cdata;
    imwrite(aa,['.' filesep 'saving' filesep 'img_' num2str(round(get(handles.slider2,'value')),'%.03d') '.tif']);
    
%     print(handles.figure1,['.' filesep 'saving' filesep 'img_' num2str(round(get(handles.slider2,'value')),'%.03d') '.tif'],'-dtiff');
end
end


% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

function centroid_scatter=draw_selection(handles)
global param
if ishandle(handles.centroid_scatter)        
    delete(handles.centroid_scatter);
end

if ~isempty(handles.selection)
    s_inds=handles.selection(:,1)>0;

    if get(handles.checkbox8,'value')==1
        centroid_scatter=scatter(handles.axis1,handles.selection(s_inds,1),handles.selection(s_inds,2),[],param.track_c_map(param.c_map_inds(handles.selection(s_inds,3)),:));
        set(centroid_scatter,'hittest','off');
    else
        centroid_scatter=[];
    end
else
   centroid_scatter=[];
end
end

function track_lines=draw_tracks(handles)
global param

try
    delete(handles.centroid_tracks_lines);
catch
end
track_lines=[];

if ~isempty(handles.selection)
    s_inds=handles.selection(:,1)>0;

    if get(handles.checkbox6,'value')==1
        
        cell_tracks_x=zeros(0);
        cell_tracks_y=zeros(0);
        t_ind=1;
        for c_ind=handles.selection(s_inds,3)'
            temp=param.tracks(c_ind).cent;
            cell_tracks_x(1:size(temp,1),t_ind)=temp(:,1);
            cell_tracks_y(1:size(temp,1),t_ind)=temp(:,2);
            t_ind=t_ind+1;
        end
        cell_tracks_x(cell_tracks_x(:)==0)=nan;
        cell_tracks_y(cell_tracks_y(:)==0)=nan;
        
        track_lines=line(cell_tracks_x,cell_tracks_y);
%         param.track_c_map(param.c_map_inds(handles.selection(s_inds,3)),:)        
        colours_to_use=param.track_c_map(param.c_map_inds(handles.selection(s_inds,3)),:);
        for t_ind=1:size(track_lines,1)
            set(track_lines(t_ind),'color',colours_to_use(t_ind,:));
        end
        
        set(track_lines,'hittest','off');
    end
end
end

function neighbour_lines=draw_neighbours(handles)
global param
try 
    delete(handles.centroid_neighbour_lines);
catch
end
neighbour_lines=[];

if ~isempty(handles.selection)
    s_inds=handles.selection(:,1)>0;

    if get(handles.checkbox10,'value')==1
        neighbours_x=zeros(2,100);
        neighbours_y=zeros(2,100);
        n_ind=1;        
        
        all_cell_inds=handles.selection(s_inds,3)';
        for c_ind1_index=1:length(all_cell_inds)
            c_ind1=all_cell_inds(c_ind1_index);
            cell_time_1=get(handles.slider2,'value')-handles.time_interval(1)+1-param.tracks(c_ind1).t(1)+1;
            for c_ind2_index=(c_ind1_index+1):length(all_cell_inds)
                c_ind2=all_cell_inds(c_ind2_index);              
                cell_time_2=get(handles.slider2,'value')-handles.time_interval(1)+1-param.tracks(c_ind2).t(1)+1;
% neighbours found
                if sum(param.tracks(c_ind1).neighs{cell_time_1}==c_ind2)>0
                    cent_1=double(param.tracks(c_ind1).cent(cell_time_1,:));
                    cent_2=double(param.tracks(c_ind2).cent(cell_time_2,:));                    
                    neighbours_x(:,n_ind)=[cent_1(1);cent_2(1)];
                    neighbours_y(:,n_ind)=[cent_1(2);cent_2(2)];
                    n_ind=n_ind+1;
                end
            end            
        end
        neighbours_x(:,n_ind:end)=[];
        neighbours_y(:,n_ind:end)=[];
        
        neighbour_lines=line(neighbours_x,neighbours_y);        
        set(neighbour_lines,'color',[1 1 0.9999]);
        
        set(neighbour_lines,'hittest','off');
    end
end
end

function polygon_coordinates=draw_polygons(handles)
global param
try 
    delete(handles.centroid_polygons);
catch
end
polygon_coordinates=[];

if ~isempty(handles.selection)
    s_inds=handles.selection(:,1)>0;

    if get(handles.checkbox11,'value')==1
        polygons_x=zeros(10,100);
        polygons_y=zeros(10,100);
        n_ind=1;        
        
        all_cell_inds=handles.selection(s_inds,3)';
        for c_ind1_index=1:length(all_cell_inds)
            c_ind1=all_cell_inds(c_ind1_index);
            cell_time_1=get(handles.slider2,'value')-handles.time_interval(1)+1-param.tracks(c_ind1).t(1)+1;
        
            temp_points=zeros(100,2);
            for b_ind=1:length(param.tracks(c_ind1).bounds{cell_time_1})                
                num_points_end=param.tracks(c_ind1).bounds{cell_time_1}{b_ind}(1,1);
                b_points_end=double(param.tracks(c_ind1).bounds{cell_time_1}{b_ind}(2:(num_points_end+1),:));
                if num_points_end==1
                    b_points_end(2,:)=b_points_end(1,:);
                    num_points_end=2;
                end
                
                distances=(repmat(b_points_end(:,1),[1 num_points_end])-repmat(b_points_end(:,1)',[num_points_end 1])).^2+...
                    (repmat(b_points_end(:,2),[1 num_points_end])-repmat(b_points_end(:,2)',[num_points_end 1])).^2;
                [vals, max_inds]=max(distances);
                [~,m_ind_j]=max(vals);
                m_ind_i=max_inds(m_ind_j);
                
                temp_points(b_ind*2-1,:)=b_points_end(m_ind_i,:);
                temp_points(b_ind*2,:)=b_points_end(m_ind_j,:);
            end
            temp_points((b_ind*2+1):end,:)=[];
            
            for b_ind=3:2:size(temp_points,1)
                find_single_points=(temp_points(b_ind:2:end,1)==temp_points((b_ind+1):2:end,1) & temp_points(b_ind:2:end,2)==temp_points((b_ind+1):2:end,2))';
                find_single_points(2,:)=0;
                find_single_points=find_single_points(:);                
                [~,m_ind]=min((temp_points(b_ind:end,1)-temp_points(b_ind-1,1)).^2+  (temp_points(b_ind:end,2)-temp_points(b_ind-1,2)).^2-find_single_points*0.01);
                m_ind=m_ind+b_ind-1;
                if mod(m_ind,2)==0
                    m_ind=m_ind-1;
                    temp_points(m_ind:(m_ind+1),:)=flipud(temp_points(m_ind:(m_ind+1),:));                    
                end
                if m_ind~=b_ind
                    temp=temp_points(b_ind:(b_ind+1),:);
                    temp_points(b_ind:(b_ind+1),:)=temp_points(m_ind:(m_ind+1),:);
                    temp_points(m_ind:(m_ind+1),:)=temp;
                end                    
            end

            polygons_x(1:size(temp_points,1)/2,n_ind)=temp_points(1:2:end,1);
            polygons_y(1:size(temp_points,1)/2,n_ind)=temp_points(1:2:end,2);
                       
            n_ind=n_ind+1;
        end
        polygons_x(:,n_ind:end)=[];
        polygons_y(:,n_ind:end)=[];
        
        for n_ind=1:size(polygons_x,2)
            polygons_x(polygons_x(:,n_ind)==0,n_ind)=polygons_x(find(polygons_x(:,n_ind),1,'last'),n_ind);
            polygons_y(polygons_y(:,n_ind)==0,n_ind)=polygons_y(find(polygons_y(:,n_ind),1,'last'),n_ind);
        end
                
        polygon_coordinates=fill(polygons_x,polygons_y,[1 0 1]);
        
        colours_to_use=param.track_c_map(param.c_map_inds(handles.selection(s_inds,3)),:);        
        for f_ind=1:length(polygon_coordinates)
            set(polygon_coordinates(f_ind),'Facecolor',colours_to_use(f_ind,:))
        end
        
        set(polygon_coordinates,'hittest','off');
    end
end
end

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);    
end
end

% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

% keyboard

if isempty(eventdata.Modifier)
    switch eventdata.Key
        case 'rightarrow'        
            if (get(handles.slider2,'value')+1)>get(handles.slider2,'max')
                set(handles.slider2,'value',get(handles.slider2,'max'));
            else
                set(handles.slider2,'value',get(handles.slider2,'value')+1);            
            end
            slider2_Callback(handles.slider2,[],handles);
        case 'leftarrow'        
            if (get(handles.slider2,'value')-1)<get(handles.slider2,'min')
                set(handles.slider2,'value',get(handles.slider2,'min'));
            else
                set(handles.slider2,'value',get(handles.slider2,'value')-1);            
            end
            slider2_Callback(handles.slider2,[],handles);        
    end
elseif strcmp(eventdata.Modifier{1},'alt')
    switch eventdata.Key
        case 'rightarrow'        
            if (get(handles.slider2,'value')+10)>get(handles.slider2,'max')
                set(handles.slider2,'value',get(handles.slider2,'max'));
            else
                set(handles.slider2,'value',get(handles.slider2,'value')+10);            
            end
            slider2_Callback(handles.slider2,[],handles);
        case 'leftarrow'        
            if (get(handles.slider2,'value')-10)<get(handles.slider2,'min')
                set(handles.slider2,'value',get(handles.slider2,'min'));
            else
                set(handles.slider2,'value',get(handles.slider2,'value')-10);            
            end
            slider2_Callback(handles.slider2,[],handles);        
    end    
end
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global param;


% keyboard

if isempty(handles.time)
    % parameter list
    param=struct('w_s',[100 100],...    
        'o_filename','',...
        'f_filename','',...
        's_filename','',...
        'filter_script','',...
        'image_interval',[1 200],...
        'image_interval_start_in_original_seq',15,...
        'img_s',[],...
        'velocity_field',[],...
        'filter_piv_result','true',...
        'moving_average_window_size',1,...
        'ingression_area_size_threshold',15,...
        'tracks',[],...
        'computing_time',[],...
        'seq_name','squareD');

    param.vel_field_scale_factor=1;
    param.img_s=[size(handles.o_imgs,2) size(handles.o_imgs,1)];
    
    param.c_map_max=16;    
    param.track_c_map=jet(param.c_map_max);

    handles.time=1;

    % initial segmentation
    f_img=handles.f_imgs(:,:,handles.time);
    f_img_m=imhmin(f_img,10);
    segment_temp=watershed(f_img_m);
    handles.s_imgs(:,:,handles.time)=uint8(~logical(segment_temp))*255;
    handles.s_imgs_independent(:,:,handles.time)=uint8(~logical(segment_temp))*255;
    handles.labeled_independent_image=segment_temp;
    info_temp=regionprops(segment_temp,'centroid','area','perimeter','orientation','majoraxislength','minoraxislength');    
    collect_initial_segmentation_bounaries(info_temp,segment_temp);
    param.c_map_inds=mod(1:length(param.tracks),param.c_map_max)'+1;
    
    handles.cents{handles.time}=zeros(length(param.tracks),1); % [x y id t_ind]
    for i_ind=1:length(param.tracks)
        handles.cents{handles.time}(i_ind,1:4)=[param.tracks(i_ind).cent(1,:) i_ind 1];
    end
    
%% make intial selection    

% %     handles.selection=handles.cents{handles.time}(...
% %         (handles.cents{handles.time}(:,1)>=(param.img_s(1)/2-diff(handles.min_max_mean(1,1:2))/2)) & ...
% %         (handles.cents{handles.time}(:,1)<=(param.img_s(1)/2+diff(handles.min_max_mean(1,1:2))/2)) & ...
% %         (handles.cents{handles.time}(:,2)>=(param.img_s(2)/2-diff(handles.min_max_mean(1,3:4))/2)) & ...
% %         (handles.cents{handles.time}(:,2)<=(param.img_s(2)/2+diff(handles.min_max_mean(1,3:4))/2)),:);
        
%%        
    
    set(handles.text3,'string',num2str(handles.time+handles.time_interval(1)-1));
    
    set(handles.slider2,'value',handles.time+handles.time_interval(1)-1);
    slider2_Callback(handles.slider2,[],handles);    
    
    if isempty(handles.selection)
        guidata(hObject, handles);
    end
    
    set(handles.pushbutton4,'enable','off');
    set(handles.pushbutton12,'enable','on');
    
elseif (handles.time+handles.time_interval(1)-1)==round(get(handles.slider2,'value'))
%       keyboard
% tracking segmentation    
    [handles.s_imgs(:,:,handles.time)]=antti_track_next_time_point(handles.f_imgs(:,:,handles.time),handles.time,handles.projected_points);

    handles.cents{handles.time}=zeros(size(handles.projected_points,1),4);
    for i_ind=1:size(handles.projected_points,1)
        handles.cents{handles.time}(i_ind,1:4)=[param.tracks(handles.projected_points(i_ind,1)).cent(end,:) handles.projected_points(i_ind,1) param.tracks(handles.projected_points(i_ind,1)).t(end)];
    end    
    
    draw_image(handles);
    handles.centroid_scatter=draw_selection(handles);
    
    guidata(hObject, handles);
    
    set(handles.pushbutton12,'enable','on');
%     set(handles.pushbutton4,'enable','off');        
end
end
    
% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
% keyboard

% % if get(hObject,'value')==1
% %     img=antti_merge_RGBG(handles.s_imgs(:,:,round(get(handles.slider2,'Value'))-handles.time_interval(1)+1),[],[],handles.o_imgs(:,:,round(get(handles.slider2,'Value'))-handles.time_interval(1)+1));
% %     set(handles.cur_img,'CData',img);
% % else
% % %     keyboard
% %     set(handles.cur_img,'CData',handles.o_imgs(:,:,round(get(handles.slider2,'Value'))-handles.time_interval(1)+1));    
% % end
draw_image(handles);
end

% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)

% keyboard

if eventdata.VerticalScrollCount>0
    if (get(handles.slider2,'value')+eventdata.VerticalScrollCount)>get(handles.slider2,'max')
        set(handles.slider2,'value',get(handles.slider2,'max'));
    else
        set(handles.slider2,'value',get(handles.slider2,'value')+eventdata.VerticalScrollCount);
    end
    slider2_Callback(handles.slider2,[],handles);
else
    if (get(handles.slider2,'value')+eventdata.VerticalScrollCount)<get(handles.slider2,'min')
        set(handles.slider2,'value',get(handles.slider2,'min'));
    else
        set(handles.slider2,'value',get(handles.slider2,'value')+eventdata.VerticalScrollCount);
    end
    slider2_Callback(handles.slider2,[],handles);
end
end
   
function position_and_button(hObject,eventdata)

global param;

Position = get( ancestor(hObject,'axes'), 'CurrentPoint' );
Button = get( ancestor(hObject,'figure'), 'SelectionType' );
%# do stuff with Position and Button   
% Position
% Button

handles=guidata(gcf);

%     keyboard
% select and unselect a track
if ~isempty(handles.time) & handles.time>=round((get(handles.slider2,'value')-handles.time_interval(1)+1)) & get(handles.checkbox7,'value')==0
    if strcmp(Button,'normal')
%         keyboard
        [~,m_ind]=min(sum((double(handles.cents{round(get(handles.slider2,'value'))-handles.time_interval(1)+1}(:,1:2))-repmat(Position(1,1:2),[length(handles.cents{round(get(handles.slider2,'Value'))-handles.time_interval(1)+1}) 1])).^2,2));
        candidate=handles.cents{round(get(handles.slider2,'Value'))-handles.time_interval(1)+1}(m_ind,:);

        if ~isempty(handles.selection) & sum(handles.selection(:,3)==candidate(1,3))>0
            handles.selection(handles.selection(:,3)==candidate(1,3),:)=[];
        else
            handles.selection(end+1,:)=candidate;
        end

        handles.centroid_scatter=draw_selection(handles);
        handles.centroid_tracks_lines=draw_tracks(handles);
        handles.centroid_neighbour_lines=draw_neighbours(handles);
        handles.centroid_polygons=draw_polygons(handles);
    end
elseif ~isempty(handles.time) & handles.time==round((get(handles.slider2,'value')-handles.time_interval(1)+1)) & get(handles.checkbox7,'value')==1 & handles.time~=1
    
    [~,m_ind]=min(sum((double(handles.projected_points(:,2:3))-repmat(Position(1,1:2),[size(handles.projected_points,1) 1])).^2,2));

    cell_i=handles.projected_points(m_ind,1);
    
%     handles.projected_points(:,2:3)=handles.projected_points(:,2:3)+round([interp2(param.vel_field_x_coords,param.vel_field_y_coords,param.velocity_field(:,:,1,handles.time),handles.projected_points(:,2),handles.projected_points(:,3),'linear') interp2(param.vel_field_x_coords,param.vel_field_y_coords,param.velocity_field(:,:,2,handles.time),handles.projected_points(:,2),handles.projected_points(:,3),'linear')]);
% 
%     handles.projected_points(handles.projected_points(:,2)<1,2)=1;
%     handles.projected_points(handles.projected_points(:,2)>param.img_s(1),2)=1;
%     handles.projected_points(handles.projected_points(:,3)<1,3)=1;
%     handles.projected_points(handles.projected_points(:,3)>param.img_s(1),3)=1;
% 
%     delete(handles.projection_scatter);
%     handles.projection_scatter=scatter(handles.axis1,handles.projected_points(:,2),handles.projected_points(:,3),[],param.track_c_map(handles.projected_points(:,4),:),'.');
%     set(handles.projection_scatter,'hittest','off');
%     if get(handles.checkbox7,'value')==0
%         set(handles.projection_scatter,'visible','off');
%     end
    
    
    if strcmp(Button,'normal') % just left click, select, add new location, move, deselect
        [x,y,key]=ginput(1);
        if key==1 % left click, move cell
            handles.projected_points(m_ind,2:3)=round([x y]);
        elseif key==97 % a pressed, add new cell        
            param.tracks(end+1).t=[];
            param.tracks(end).birth=-1;
            param.tracks(end).daughters=[0 0];
            param.tracks(end).death=0;
%             handles.selection(end+1,1:4)=[0 0 length(param.tracks) 0];
            param.c_map_inds(end+1)=mod(length(param.tracks),param.c_map_max)+1;
            handles.projected_points(end+1,:)=[length(param.tracks) round(Position(1,1:2)) param.c_map_inds(end)];        
        end
        
    elseif strcmp(Button,'extend') % divisions, select mother, add daughter cell locations, deselect
        param.tracks(cell_i).death=-3;
        if param.tracks(cell_i).t(end)==handles.time
            param.tracks(cell_i).t(end)=[];
            param.tracks(cell_i).A(end)=[];
            param.tracks(cell_i).cent(end,:)=[];
            param.tracks(cell_i).ellipse(end,:)=[];
            param.tracks(cell_i).perim(end)=[];
            param.tracks(cell_i).neighs(end)=[];
            param.tracks(cell_i).bounds(end)=[];            
        end     
        
        param.tracks(cell_i).daughters=[length(param.tracks)+1 length(param.tracks)+2];
        
        param.tracks(end+1).t=[];              
        param.tracks(end).birth=cell_i;
        param.tracks(end).daughters=[0 0];
        param.tracks(end).death=0;        
        
        param.tracks(end+1).t=[];
        param.tracks(end).birth=cell_i;
        param.tracks(end).daughters=[0 0];
        param.tracks(end).death=0;
        
        [x,y]=ginput(2);
        
        param.c_map_inds((end+1):(end+2))=param.c_map_inds(cell_i);
        handles.projected_points((end+1):(end+2),:)=[[length(param.tracks)-1 round([x(1) y(1)]) param.c_map_inds(end-1)];[length(param.tracks) round([x(2) y(2)]) param.c_map_inds(end)]];        
        handles.projected_points(m_ind,:)=[];
        
        if ~isempty(handles.selection) & isempty(param.tracks(cell_i).t)
            handles.selection(handles.selection(:,3)==cell_i,:)=[];
        end

% select daughter only if mother is selected        
        if ~isempty(handles.selection) & sum(handles.selection(:,3)==cell_i)>0
            handles.selection(end+1,1:4)=[0 0 length(param.tracks)-1 0];
            handles.selection(end+1,1:4)=[0 0 length(param.tracks) 0];
        end
        
    elseif strcmp(Button,'alt') % ingression: remove marker
        param.tracks(cell_i).death=-1;
        if ~isempty(param.tracks(cell_i).t) & param.tracks(cell_i).t(end)==handles.time
            param.tracks(cell_i).t(end)=[];
            param.tracks(cell_i).A(end)=[];
            param.tracks(cell_i).cent(end,:)=[];
            param.tracks(cell_i).ellipse(end,:)=[];
            param.tracks(cell_i).perim(end)=[];
            param.tracks(cell_i).neighs(end)=[];
            param.tracks(cell_i).bounds(end)=[];            
        end        
        handles.projected_points(m_ind,:)=[];
        if ~isempty(handles.selection) & isempty(param.tracks(cell_i).t)
            handles.selection(handles.selection(:,3)==cell_i,:)=[];
        end
% %     elseif strcmp(Button,'open') % appearing cell, introduce new marker, do not select                
% %         param.tracks(end+1).t=[];
% %         param.tracks(end).birth=-1;
% %         handles.selection(end+1,1:4)=[0 0 length(param.tracks) 0];
% %         param.c_map_inds(end+1)=mod(length(param.tracks),param.c_map_max)+1;
% %         handles.projected_points(end+1,:)=[length(param.tracks) round(Position(1,1:2)) param.c_map_inds(end)];
    end
    
%     redraw markers

    handles.projection_scatter=draw_projection(handles);
end

guidata(handles.figure1, handles);
end   
    

% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
draw_image(handles);
end


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4
draw_image(handles);
end


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5
end

% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6
handles.centroid_tracks_lines=draw_tracks(handles);
% guidata(handles.figure1, handles);
handles.centroid_neighbour_lines=draw_neighbours(handles);
handles.centroid_polygons=draw_polygons(handles);
guidata(handles.figure1, handles);
end

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles1)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global param;
handles1.slider2_value=get(handles1.slider2,'value');


[filename, pathname, index]=uiputfile('*.mat');

if index~=0
    save([pathname filename],'handles1','param');    
end
end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles1)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

choice = questdlg('Close tracker?','','Yes','No','Yes');
% Handle response
switch choice
    case 'Yes'
        choice1 = questdlg('Save data?','', 'Yes','No','Yes');
        switch choice1
            case 'Yes'
                global param;
                handles1.slider2_value=get(handles1.slider2,'value');
                
                [filename, pathname, index]=uiputfile('*.mat');
                
                if index~=0
                    save([pathname filename],'handles1','param');
                end
        end
        delete(hObject);
end
end
% Hint: delete(hObject) closes the figure


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global param
global param1
param1=param;

keyboard
end

% guidata(handles.figure1, handles);


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles) % project points
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global param;

if handles.time>=(handles.time_interval(2)-handles.time_interval(1)+1)
    return
end
% keyboard
[param.vel_field_x_coords, param.vel_field_y_coords, u_filtered, v_filtered]=antti_PIVlab_vel_field(imresize(handles.f_imgs(:,:,handles.time),param.vel_field_scale_factor),imresize(handles.f_imgs(:,:,handles.time+1),param.vel_field_scale_factor),0);

param.vel_field_x_coords=param.vel_field_x_coords/param.vel_field_scale_factor;
param.vel_field_y_coords=param.vel_field_y_coords/param.vel_field_scale_factor;
param.velocity_field(:,:,1,handles.time)=u_filtered/param.vel_field_scale_factor;
param.velocity_field(:,:,2,handles.time)=v_filtered/param.vel_field_scale_factor;

handles.projected_points=[];
for cell_i=1:length(param.tracks)
    if param.tracks(cell_i).death==0
        handles.projected_points(end+1,1:4)=[cell_i double(param.tracks(cell_i).cent(end,:)) param.c_map_inds(cell_i)];
    end
end

handles.projected_points(:,2:3)=handles.projected_points(:,2:3)+round([interp2(param.vel_field_x_coords,param.vel_field_y_coords,param.velocity_field(:,:,1,handles.time),handles.projected_points(:,2),handles.projected_points(:,3),'linear') interp2(param.vel_field_x_coords,param.vel_field_y_coords,param.velocity_field(:,:,2,handles.time),handles.projected_points(:,2),handles.projected_points(:,3),'linear')]);

% ingress points flowing out

inds_out=(handles.projected_points(:,2)<1) | (handles.projected_points(:,2)>param.img_s(1)) | (handles.projected_points(:,3)<1) | (handles.projected_points(:,3)>param.img_s(2));

for out_ind=find(inds_out)'
    param.tracks(handles.projected_points(out_ind,1)).death=-2;
    if ~isempty(handles.selection) & isempty(param.tracks(handles.projected_points(out_ind,1)).t)
        handles.selection(handles.selection(:,3)==handles.projected_points(out_ind,1),:)=[];
    end
end

handles.projected_points(inds_out,:)=[];

% check that distance between points is sufficient

mask=zeros(fliplr(param.img_s)+[2 2]);
for i_ind=1:size(handles.projected_points,1)
    mask(handles.projected_points(i_ind,3)+1,handles.projected_points(i_ind,2)+1)=mask(handles.projected_points(i_ind,3)+1,handles.projected_points(i_ind,2)+1)+1;
end

for i_ind=1:(size(handles.projected_points,1)-1)
    while sum(sum(mask((handles.projected_points(i_ind,3)+1-1):(handles.projected_points(i_ind,3)+1+1),(handles.projected_points(i_ind,2)+1-1):(handles.projected_points(i_ind,2)+1+1))))>1
        mask(handles.projected_points(i_ind,3)+1,handles.projected_points(i_ind,2)+1)=mask(handles.projected_points(i_ind,3)+1,handles.projected_points(i_ind,2)+1)-1;

        handles.projected_points(i_ind,3)=handles.projected_points(i_ind,3)+sign(rand(1)-0.5)*2;
        if handles.projected_points(i_ind,3)<2, handles.projected_points(i_ind,3)=2; end
        if handles.projected_points(i_ind,3)>(param.img_s(2)+1), handles.projected_points(i_ind,3)=param.img_s(2)+1; end
    
        handles.projected_points(i_ind,2)=handles.projected_points(i_ind,2)+sign(rand(1)-0.5)*2;
        if handles.projected_points(i_ind,2)<2, handles.projected_points(i_ind,2)=2; end
        if handles.projected_points(i_ind,2)>(param.img_s(1)+1), handles.projected_points(i_ind,2)=param.img_s(1)+1; end        
        
        mask(handles.projected_points(i_ind,3)+1,handles.projected_points(i_ind,2)+1)=mask(handles.projected_points(i_ind,3)+1,handles.projected_points(i_ind,2)+1)+1;
    end        
end

handles.projection_scatter=draw_projection(handles);

f_img=handles.f_imgs(:,:,handles.time+1);
f_img_m=imhmin(f_img,10);
handles.s_imgs_independent(:,:,handles.time+1)=uint8(~logical(watershed(f_img_m)))*255;

handles.time=handles.time+1;
set(handles.text3,'string',num2str(handles.time+handles.time_interval(1)-1));

set(handles.slider2,'value',get(handles.slider2,'value')+1);
slider2_Callback(handles.slider2,[],handles);

set(handles.pushbutton12,'enable','off');
set(handles.pushbutton4,'enable','on');
guidata(handles.figure1, handles);

% if auto tracking is checked segment automatically
if get(handles.checkbox5,'value')==1
% segment
    pushbutton4_Callback(hObject, eventdata, handles)
    handles=guidata(gcf);
% refresh selection    
    slider2_Callback(hObject, eventdata, handles);
    handles=guidata(gcf);    
end
handles.centroid_scatter=draw_selection(handles);
guidata(handles.figure1, handles);
end


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.checkbox7,'value')==0
    set(handles.projection_scatter,'visible','off');
elseif (round(get(handles.slider2,'Value'))-handles.time_interval(1)+1)==handles.time
    set(handles.projection_scatter,'visible','on');
end
end
% Hint: get(hObject,'Value') returns toggle state of checkbox7


% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% keyboard
% Hint: get(hObject,'Value') returns toggle state of checkbox8
if get(handles.checkbox8,'value')==0
    set(handles.centroid_scatter,'visible','off');
else        
    handles.centroid_scatter=draw_selection(handles);
    set(handles.centroid_scatter,'visible','on');
    guidata(handles.figure1, handles);
end
end

function draw_image(handles)
img=antti_merge_RGBG(handles.s_imgs(:,:,round(get(handles.slider2,'Value'))-handles.time_interval(1)+1)*uint8(get(handles.checkbox1,'value')==1),handles.s_imgs_independent(:,:,round(get(handles.slider2,'Value'))-handles.time_interval(1)+1)*uint8(get(handles.checkbox3,'value')==1),[],handles.o_imgs(:,:,round(get(handles.slider2,'Value'))-handles.time_interval(1)+1)*uint8(get(handles.checkbox4,'value')==1));
set(handles.cur_img,'CData',img);
end

function projection_scatter=draw_projection(handles)
global param

delete(handles.projection_scatter);
if ~isempty(handles.projected_points)
    projection_scatter=scatter(handles.axis1,handles.projected_points(:,2),handles.projected_points(:,3),[],param.track_c_map(handles.projected_points(:,4),:),'.');
    set(projection_scatter,'hittest','off');
else
    projection_scatter=[];
end
if get(handles.checkbox7,'value')==0
    set(projection_scatter,'visible','off');
end
end


% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9
end

function collect_initial_segmentation_bounaries(info_temp,segment_temp)

global param;

param.tracks=struct('t',uint16(1),'A',uint32(info_temp(1).Area),'cent',uint16(info_temp(1).Centroid),'ellipse',int16([info_temp(1).MajorAxisLength info_temp(1).MinorAxisLength info_temp(1).Orientation]),'perim',uint32(info_temp(1).Perimeter),'birth',int32(-1),'death',int16(0),'daughters',uint32([0 0]),'neighs',cell(1),'bounds',cell(1));
param.tracks(1).neighs=cell(1);param.tracks(1).bounds=cell(1);
for cell_i=2:size(info_temp,1)
    param.tracks(cell_i)=struct('t',uint16(1),'A',uint32(info_temp(cell_i).Area),'cent',uint16(info_temp(cell_i).Centroid),'ellipse',int16([info_temp(cell_i).MajorAxisLength info_temp(cell_i).MinorAxisLength info_temp(cell_i).Orientation]),'perim',uint32(info_temp(cell_i).Perimeter),'birth',int32(-1),'death',int16(0),'daughters',uint32([0 0]),'neighs',cell(1),'bounds',cell(1));
    param.tracks(cell_i).neighs=cell(1);param.tracks(cell_i).bounds=cell(1);
end

% membrane edges
edge_ind=max(segment_temp(:))+1;
segment_temp1=zeros(size(segment_temp)+[4 4]);
segment_temp1(3:end-2,3:end-2)=segment_temp;
segment_temp1([1 end],:)=edge_ind;
segment_temp1(:,[1 end])=edge_ind;

m_inds=find(segment_temp1==0);
m_points=[segment_temp1(m_inds-param.img_s(2)-4-1) segment_temp1(m_inds-param.img_s(2)-4) segment_temp1(m_inds-param.img_s(2)-4+1)...
    segment_temp1(m_inds-1)  segment_temp1(m_inds) segment_temp1(m_inds+1)...
    segment_temp1(m_inds+param.img_s(2)+4-1) segment_temp1(m_inds+param.img_s(2)+4) segment_temp1(m_inds+param.img_s(2)+4+1)];
m_sort=sort(m_points,2);
m_diff=diff(sort(m_points,2),1,2)>0;

% search double membrane points
ind2=find(sum(m_diff,2)==2);
[x_ind y_ind]=find(m_diff(ind2,:)');
n_pairs=reshape(m_sort(sub2ind([size(m_diff,1) 9],ind2(y_ind),x_ind+1)),[2 size(ind2,1)])';
sort_pairs=sortrows([n_pairs [floor(m_inds(ind2)/(param.img_s(2)+4))+1-2 mod(m_inds(ind2)-1,(param.img_s(2)+4))+1-2]],[1 2]);
[uni_pairs, f_locs]=unique(sort_pairs(:,1:2),'rows','first');
[~, l_locs]=unique(sort_pairs(:,1:2),'rows','last');

% save edge pixels and neighbours to the tracklet struct
for e_ind=1:size(uni_pairs,1)
    % first half
    if uni_pairs(e_ind,1)~=edge_ind
        if isempty(param.tracks(uni_pairs(e_ind,1)).bounds{1})
            param.tracks(uni_pairs(e_ind,1)).bounds{1}=cell(0);
        end
        % mark neighbours (zero if neighbour to edge of the image)
        param.tracks(uni_pairs(e_ind,1)).neighs{1}(end+1)=uint32((uni_pairs(e_ind,2)~=edge_ind)*uni_pairs(e_ind,2));
        % save membrane edge
        param.tracks(uni_pairs(e_ind,1)).bounds{1}{end+1}=uint16([0 0;sort_pairs(f_locs(e_ind):l_locs(e_ind),3:4)]);
    end

    % second half
    if uni_pairs(e_ind,2)~=edge_ind
        if isempty(param.tracks(uni_pairs(e_ind,2)).bounds{1})
            param.tracks(uni_pairs(e_ind,2)).bounds{1}=cell(0);
        end
        % mark neighbours (zero if neighbour to edge of the image)
        param.tracks(uni_pairs(e_ind,2)).neighs{1}(end+1)=uint32((uni_pairs(e_ind,1)~=edge_ind)*uni_pairs(e_ind,1));
        % save membrane edge
        param.tracks(uni_pairs(e_ind,2)).bounds{1}{end+1}=uint16([0 0;sort_pairs(f_locs(e_ind):l_locs(e_ind),3:4)]);
    end
end

% search quadruple and triple membrane points (membrane vertices)

% quadruple
ind4=find(sum(m_diff,2)==4);
[x_ind y_ind]=find(m_diff(ind4,:)');
n_quadruple=[reshape(m_sort(sub2ind([size(m_diff,1) 9],ind4(y_ind),x_ind+1)),[4 size(ind4,1)])' [floor(m_inds(ind4)/(param.img_s(2)+4))+1-2 mod(m_inds(ind4)-1,(param.img_s(2)+4))+1-2]];

% triple
ind3=find(sum(m_diff,2)==3);
[x_ind y_ind]=find(m_diff(ind3,:)');
n_triples=[reshape(m_sort(sub2ind([size(m_diff,1) 9],ind3(y_ind),x_ind+1)),[3 size(ind3,1)])' [floor(m_inds(ind3)/(param.img_s(2)+4))+1-2 mod(m_inds(ind3)-1,(param.img_s(2)+4))+1-2]];

% combine quadruple points to the triple points
n_triples=[n_triples;n_quadruple(:,[1 2 3 5 6]);n_quadruple(:,[1 2 4 5 6]);n_quadruple(:,[1 3 4 5 6]);];

for row_ind=1:size(n_triples,1)
    for col_ind=2:3
        % first half
        if n_triples(row_ind,1)~=edge_ind
            n_ind=find(param.tracks(n_triples(row_ind,1)).neighs{end}==((n_triples(row_ind,col_ind)~=edge_ind)*n_triples(row_ind,col_ind)));
            if isempty(n_ind)
                param.tracks(n_triples(row_ind,1)).bounds{end}{end+1}=uint16([1 0;n_triples(row_ind,4:5)]);
                param.tracks(n_triples(row_ind,1)).neighs{end}(end+1)=uint32(n_triples(row_ind,col_ind));
            else
                param.tracks(n_triples(row_ind,1)).bounds{end}{n_ind}(1,1)=param.tracks(n_triples(row_ind,1)).bounds{end}{n_ind}(1,1)+1;
                ind2add=param.tracks(n_triples(row_ind,1)).bounds{end}{n_ind}(1,1);
                param.tracks(n_triples(row_ind,1)).bounds{end}{n_ind}=[param.tracks(n_triples(row_ind,1)).bounds{end}{n_ind}(1:ind2add,:) ; uint16(n_triples(row_ind,4:5)) ; param.tracks(n_triples(row_ind,1)).bounds{end}{n_ind}((ind2add+1):end,:)];
            end
        end
        % second half
        if n_triples(row_ind,col_ind)~=edge_ind
            n_ind=find(param.tracks(n_triples(row_ind,col_ind)).neighs{end}==((n_triples(row_ind,1)~=edge_ind)*n_triples(row_ind,1)));
            if isempty(n_ind)
                param.tracks(n_triples(row_ind,col_ind)).bounds{end}{end+1}=uint16([1 0;n_triples(row_ind,4:5)]);
                param.tracks(n_triples(row_ind,col_ind)).neighs{end}(end+1)=uint32(n_triples(row_ind,1));
            else
                param.tracks(n_triples(row_ind,col_ind)).bounds{end}{n_ind}(1,1)=param.tracks(n_triples(row_ind,col_ind)).bounds{end}{n_ind}(1,1)+1;
                ind2add=param.tracks(n_triples(row_ind,col_ind)).bounds{end}{n_ind}(1,1);
                param.tracks(n_triples(row_ind,col_ind)).bounds{end}{n_ind}=[param.tracks(n_triples(row_ind,col_ind)).bounds{end}{n_ind}(1:ind2add,:) ; uint16(n_triples(row_ind,4:5)) ; param.tracks(n_triples(row_ind,col_ind)).bounds{end}{n_ind}((ind2add+1):end,:)];
            end
        end
    end
end

for row_ind=1:size(n_triples,1)
    % first half
    if n_triples(row_ind,2)~=edge_ind
        n_ind=find(param.tracks(n_triples(row_ind,2)).neighs{end}==((n_triples(row_ind,3)~=edge_ind)*n_triples(row_ind,3)));
        if isempty(n_ind)
            param.tracks(n_triples(row_ind,2)).bounds{end}{end+1}=uint16([1 0;n_triples(row_ind,4:5)]);
            param.tracks(n_triples(row_ind,2)).neighs{end}(end+1)=uint32(n_triples(row_ind,3));
        else
            param.tracks(n_triples(row_ind,2)).bounds{end}{n_ind}(1,1)=param.tracks(n_triples(row_ind,2)).bounds{end}{n_ind}(1,1)+1;
            ind2add=param.tracks(n_triples(row_ind,2)).bounds{end}{n_ind}(1,1);
            param.tracks(n_triples(row_ind,2)).bounds{end}{n_ind}=[param.tracks(n_triples(row_ind,2)).bounds{end}{n_ind}(1:ind2add,:) ; uint16(n_triples(row_ind,4:5)) ; param.tracks(n_triples(row_ind,2)).bounds{end}{n_ind}((ind2add+1):end,:)];
        end
    end
    % second half
    if n_triples(row_ind,3)~=edge_ind
        n_ind=find(param.tracks(n_triples(row_ind,3)).neighs{end}==((n_triples(row_ind,2)~=edge_ind)*n_triples(row_ind,2)));
        if isempty(n_ind)
            param.tracks(n_triples(row_ind,3)).bounds{end}{end+1}=uint16([1 0 ; n_triples(row_ind,4:5)]);
            param.tracks(n_triples(row_ind,3)).neighs{end}(end+1)=uint32(n_triples(row_ind,2));
        else
            param.tracks(n_triples(row_ind,3)).bounds{end}{n_ind}(1,1)=param.tracks(n_triples(row_ind,3)).bounds{end}{n_ind}(1,1)+1;
            ind2add=param.tracks(n_triples(row_ind,3)).bounds{end}{n_ind}(1,1);
            param.tracks(n_triples(row_ind,3)).bounds{end}{n_ind}=[param.tracks(n_triples(row_ind,3)).bounds{end}{n_ind}(1:ind2add,:) ; uint16(n_triples(row_ind,4:5)) ; param.tracks(n_triples(row_ind,3)).bounds{end}{n_ind}((ind2add+1):end,:)];
        end
    end
end

% 2x2 pixel squares are problem to the endpoint detection, thus these
% points must be marked as separate end points
square_inds=find(conv2(single(~logical(segment_temp1)),[1 1;1 1],'same')==4);
square_inds=[square_inds;square_inds+1;square_inds+param.img_s(2)+4;square_inds+param.img_s(2)+4+1];

m_points=[segment_temp1(square_inds-param.img_s(2)-4-1) segment_temp1(square_inds-param.img_s(2)-4) segment_temp1(square_inds-param.img_s(2)-4+1)...
    segment_temp1(square_inds-1)  segment_temp1(square_inds) segment_temp1(square_inds+1)...
    segment_temp1(square_inds+param.img_s(2)+4-1) segment_temp1(square_inds+param.img_s(2)+4) segment_temp1(square_inds+param.img_s(2)+4+1)];
m_sort=sort(m_points,2);
m_diff=diff(sort(m_points,2),1,2)>0;

[x_ind y_ind]=find(m_diff');
n_pairs=[reshape(m_sort(sub2ind([size(m_diff,1) 9],y_ind,x_ind+1)),[2 size(m_sort,1)])' [floor(square_inds/(param.img_s(2)+4))+1-2 mod(square_inds-1,(param.img_s(2)+4))+1-2]];


% save edge pixels and neighbours to the tracklet struct
for e_ind=1:size(n_pairs,1)
    % first half
    n_ind=find(param.tracks(n_pairs(e_ind,1)).neighs{end}==n_pairs(e_ind,2));
    param.tracks(n_pairs(e_ind,1)).bounds{end}{n_ind}(1,1)=param.tracks(n_pairs(e_ind,1)).bounds{end}{n_ind}(1,1)+1;
    ind2add=param.tracks(n_pairs(e_ind,1)).bounds{end}{n_ind}(1,1);
    param.tracks(n_pairs(e_ind,1)).bounds{end}{n_ind}=[param.tracks(n_pairs(e_ind,1)).bounds{end}{n_ind}(1:ind2add,:) ; uint16(n_pairs(e_ind,3:4)) ; param.tracks(n_pairs(e_ind,1)).bounds{end}{n_ind}((ind2add+1):end,:)];

    % second half
    n_ind=find(param.tracks(n_pairs(e_ind,2)).neighs{end}==n_pairs(e_ind,1));
    param.tracks(n_pairs(e_ind,2)).bounds{end}{n_ind}(1,1)=param.tracks(n_pairs(e_ind,2)).bounds{end}{n_ind}(1,1)+1;
    ind2add=param.tracks(n_pairs(e_ind,2)).bounds{end}{n_ind}(1,1);
    param.tracks(n_pairs(e_ind,2)).bounds{end}{n_ind}=[param.tracks(n_pairs(e_ind,2)).bounds{end}{n_ind}(1:ind2add,:) ; uint16(n_pairs(e_ind,3:4)) ; param.tracks(n_pairs(e_ind,2)).bounds{end}{n_ind}((ind2add+1):end,:)];
end

end

%%

function [s_img]=antti_track_next_time_point(f_img,time,projected_points)

global param;
inds_alive=projected_points;
bw_mask=zeros(fliplr(param.img_s));
bw_mask(sub2ind(fliplr(param.img_s),double(projected_points(:,3)),double(projected_points(:,2))))=1;

%%
for c_ind=1:size(projected_points,1) 
    projected_points(c_ind,5)=max(1,floor(sqrt(double(param.tracks(c_ind).A(end)/pi))*0.4));
end
for d_ind=2:max(projected_points(:,5))
    projected_points(projected_points(:,5)<d_ind,:)=[];
    for c_ind=1:size(projected_points,1)
        if (projected_points(c_ind,2)+d_ind-1)>=size(f_img,2) | (projected_points(c_ind,3)+d_ind-1)>=size(f_img,1) | (projected_points(c_ind,2)-d_ind+1)<=1 | (projected_points(c_ind,3)-d_ind+1)<=1
            continue
        end
        % right
        x_point=projected_points(c_ind,2)+d_ind;
        y_point=projected_points(c_ind,3);        
        if sum(bw_mask((y_point-1):(y_point+1),x_point))==0
            bw_mask(y_point,x_point-1)=1;
        end
        
        % left
        x_point=projected_points(c_ind,2)-d_ind;
        y_point=projected_points(c_ind,3);        
        if sum(bw_mask((y_point-1):(y_point+1),x_point))==0
            bw_mask(y_point,x_point+1)=1;
        end  
        
        % top
        x_point=projected_points(c_ind,2);
        y_point=projected_points(c_ind,3)-d_ind;        
        if sum(bw_mask(y_point,(x_point-1):(x_point+1)))==0
            bw_mask(y_point+1,x_point)=1;
        end        
        
        % bottom
        x_point=projected_points(c_ind,2);
        y_point=projected_points(c_ind,3)+d_ind;        
        if sum(bw_mask(y_point,(x_point-1):(x_point+1)))==0
            bw_mask(y_point-1,x_point)=1;
        end                
    end
    
end

%%

segment_temp=(watershed(imimposemin(f_img,bw_mask)));

s_img=uint8(~logical(segment_temp))*255;

info_temp=regionprops(segment_temp,'centroid','area','perimeter','orientation','majoraxislength','minoraxislength');

%% assign new cell parameters to existing struct
for cell_i=1:size(inds_alive,1)
    info_ind=segment_temp(inds_alive(cell_i,3),inds_alive(cell_i,2));
    inds_alive(cell_i,4)=info_ind;
    t_ind=1;
    if ~isempty(param.tracks(inds_alive(cell_i,1)).t)
        t_ind=time-param.tracks(inds_alive(cell_i,1)).t(1)+1;   
    end    
    param.tracks(inds_alive(cell_i,1)).t(t_ind)=time;
    param.tracks(inds_alive(cell_i,1)).A(t_ind)=uint32(info_temp(info_ind).Area);
    param.tracks(inds_alive(cell_i,1)).cent(t_ind,1:2)=uint16(info_temp(info_ind).Centroid);
    param.tracks(inds_alive(cell_i,1)).ellipse(t_ind,1:3)=int16([info_temp(info_ind).MajorAxisLength info_temp(info_ind).MinorAxisLength info_temp(info_ind).Orientation]);
    param.tracks(inds_alive(cell_i,1)).perim(t_ind)=uint32(info_temp(info_ind).Perimeter);
    param.tracks(inds_alive(cell_i,1)).bounds{t_ind}=cell(0);
    param.tracks(inds_alive(cell_i,1)).neighs{t_ind}=zeros(0);
end

% membrane edges
edge_ind=max(segment_temp(:))+1;
segment_temp1=zeros(size(segment_temp)+[4 4]);
segment_temp1(3:end-2,3:end-2)=segment_temp;
segment_temp1([1 end],:)=edge_ind;
segment_temp1(:,[1 end])=edge_ind;

m_inds=find(segment_temp1==0);
m_points=[segment_temp1(m_inds-param.img_s(2)-4-1) segment_temp1(m_inds-param.img_s(2)-4) segment_temp1(m_inds-param.img_s(2)-4+1)...
    segment_temp1(m_inds-1)  segment_temp1(m_inds) segment_temp1(m_inds+1)...
    segment_temp1(m_inds+param.img_s(2)+4-1) segment_temp1(m_inds+param.img_s(2)+4) segment_temp1(m_inds+param.img_s(2)+4+1)];
m_sort=sort(m_points,2);
m_diff=diff(sort(m_points,2),1,2)>0;

% search double membrane points
ind2=find(sum(m_diff,2)==2);
[x_ind y_ind]=find(m_diff(ind2,:)');
n_pairs=reshape(m_sort(sub2ind([size(m_diff,1) 9],ind2(y_ind),x_ind+1)),[2 size(ind2,1)])';
sort_pairs=sortrows([n_pairs [floor(m_inds(ind2)/(param.img_s(2)+4))+1-2 mod(m_inds(ind2)-1,(param.img_s(2)+4))+1-2]],[1 2]);
[uni_pairs, f_locs]=unique(sort_pairs(:,1:2),'rows','first');
[~, l_locs]=unique(sort_pairs(:,1:2),'rows','last');

inds_alive=sortrows(inds_alive,4);
inds_alive(end+1,1)=0;

% save edge pixels and neighbours to the tracklet struct
for e_ind=1:size(uni_pairs,1)
    % first half
    if uni_pairs(e_ind,1)~=edge_ind
        % mark neighbours (zero if neighbour to edge of the image)
        param.tracks(inds_alive(uni_pairs(e_ind,1),1)).neighs{end}(end+1)=uint32(inds_alive(uni_pairs(e_ind,2),1));
        % save membrane edge
        param.tracks(inds_alive(uni_pairs(e_ind,1),1)).bounds{end}{end+1}=uint16([0 0;sort_pairs(f_locs(e_ind):l_locs(e_ind),3:4)]);
    end
    
    % second half
    if uni_pairs(e_ind,2)~=edge_ind
        % mark neighbours (zero if neighbour to edge of the image)
        param.tracks(inds_alive(uni_pairs(e_ind,2),1)).neighs{end}(end+1)=uint32(inds_alive(uni_pairs(e_ind,1),1));
        % save membrane edge
        param.tracks(inds_alive(uni_pairs(e_ind,2),1)).bounds{end}{end+1}=uint16([0 0;sort_pairs(f_locs(e_ind):l_locs(e_ind),3:4)]);
    end
end

% search quadruple and triple membrane points (membrane vertices)

%     keyboard

% quadruple
ind4=find(sum(m_diff,2)==4);
[x_ind y_ind]=find(m_diff(ind4,:)');
n_quadruple=[reshape(m_sort(sub2ind([size(m_diff,1) 9],ind4(y_ind),x_ind+1)),[4 size(ind4,1)])' [floor(m_inds(ind4)/(param.img_s(2)+4))+1-2 mod(m_inds(ind4)-1,(param.img_s(2)+4))+1-2]];

% triple
ind3=find(sum(m_diff,2)==3);
[x_ind y_ind]=find(m_diff(ind3,:)');
n_triples=[reshape(m_sort(sub2ind([size(m_diff,1) 9],ind3(y_ind),x_ind+1)),[3 size(ind3,1)])' [floor(m_inds(ind3)/(param.img_s(2)+4))+1-2 mod(m_inds(ind3)-1,(param.img_s(2)+4))+1-2]];

% combine quadruple points to the triple points
n_triples=[n_triples;n_quadruple(:,[1 2 3 5 6]);n_quadruple(:,[1 2 4 5 6]);n_quadruple(:,[1 3 4 5 6]);];

for row_ind=1:size(n_triples,1)
    for col_ind=2:3
        % first half
        if n_triples(row_ind,1)~=edge_ind
            n_ind=find(param.tracks(inds_alive(n_triples(row_ind,1),1)).neighs{end}==inds_alive(n_triples(row_ind,col_ind),1));
            if isempty(n_ind)
                param.tracks(inds_alive(n_triples(row_ind,1),1)).bounds{end}{end+1}=uint16([1 0;n_triples(row_ind,4:5)]);
                param.tracks(inds_alive(n_triples(row_ind,1),1)).neighs{end}(end+1)=uint32(inds_alive(n_triples(row_ind,col_ind),1));
            else
                param.tracks(inds_alive(n_triples(row_ind,1),1)).bounds{end}{n_ind}(1,1)=param.tracks(inds_alive(n_triples(row_ind,1),1)).bounds{end}{n_ind}(1,1)+1;
                ind2add=param.tracks(inds_alive(n_triples(row_ind,1),1)).bounds{end}{n_ind}(1,1);
                param.tracks(inds_alive(n_triples(row_ind,1),1)).bounds{end}{n_ind}=[param.tracks(inds_alive(n_triples(row_ind,1),1)).bounds{end}{n_ind}(1:ind2add,:) ; uint16(n_triples(row_ind,4:5)) ; param.tracks(inds_alive(n_triples(row_ind,1),1)).bounds{end}{n_ind}((ind2add+1):end,:)];
            end
        end
        % second half
        if n_triples(row_ind,col_ind)~=edge_ind
            n_ind=find(param.tracks(inds_alive(n_triples(row_ind,col_ind),1)).neighs{end}==inds_alive(n_triples(row_ind,1),1));
            if isempty(n_ind)
                param.tracks(inds_alive(n_triples(row_ind,col_ind),1)).bounds{end}{end+1}=uint16([1 0;n_triples(row_ind,4:5)]);
                param.tracks(inds_alive(n_triples(row_ind,col_ind),1)).neighs{end}(end+1)=uint32(inds_alive(n_triples(row_ind,1),1));
            else
                param.tracks(inds_alive(n_triples(row_ind,col_ind),1)).bounds{end}{n_ind}(1,1)=param.tracks(inds_alive(n_triples(row_ind,col_ind),1)).bounds{end}{n_ind}(1,1)+1;
                ind2add=param.tracks(inds_alive(n_triples(row_ind,col_ind),1)).bounds{end}{n_ind}(1,1);
                param.tracks(inds_alive(n_triples(row_ind,col_ind),1)).bounds{end}{n_ind}=[param.tracks(inds_alive(n_triples(row_ind,col_ind),1)).bounds{end}{n_ind}(1:ind2add,:) ; uint16(n_triples(row_ind,4:5)) ; param.tracks(inds_alive(n_triples(row_ind,col_ind),1)).bounds{end}{n_ind}((ind2add+1):end,:)];
            end
        end
    end
end

for row_ind=1:size(n_triples,1)
    % first half
    if n_triples(row_ind,2)~=edge_ind
        n_ind=find(param.tracks(inds_alive(n_triples(row_ind,2),1)).neighs{end}==inds_alive(n_triples(row_ind,3),1));
        if isempty(n_ind)
            param.tracks(inds_alive(n_triples(row_ind,2),1)).bounds{end}{end+1}=uint16([1 0;n_triples(row_ind,4:5)]);
            param.tracks(inds_alive(n_triples(row_ind,2),1)).neighs{end}(end+1)=uint32(inds_alive(n_triples(row_ind,3),1));
        else
            param.tracks(inds_alive(n_triples(row_ind,2),1)).bounds{end}{n_ind}(1,1)=param.tracks(inds_alive(n_triples(row_ind,2),1)).bounds{end}{n_ind}(1,1)+1;
            ind2add=param.tracks(inds_alive(n_triples(row_ind,2),1)).bounds{end}{n_ind}(1,1);
            param.tracks(inds_alive(n_triples(row_ind,2),1)).bounds{end}{n_ind}=[param.tracks(inds_alive(n_triples(row_ind,2),1)).bounds{end}{n_ind}(1:ind2add,:) ; uint16(n_triples(row_ind,4:5)) ; param.tracks(inds_alive(n_triples(row_ind,2),1)).bounds{end}{n_ind}((ind2add+1):end,:)];
        end
    end
    % second half
    if n_triples(row_ind,3)~=edge_ind
        n_ind=find(param.tracks(inds_alive(n_triples(row_ind,3),1)).neighs{end}==inds_alive(n_triples(row_ind,2),1));
        if isempty(n_ind)
            param.tracks(inds_alive(n_triples(row_ind,3),1)).bounds{end}{end+1}=uint16([1 0;n_triples(row_ind,4:5)]);
            param.tracks(inds_alive(n_triples(row_ind,3),1)).neighs{end}(end+1)=uint32(inds_alive(n_triples(row_ind,2),1));
        else
            param.tracks(inds_alive(n_triples(row_ind,3),1)).bounds{end}{n_ind}(1,1)=param.tracks(inds_alive(n_triples(row_ind,3),1)).bounds{end}{n_ind}(1,1)+1;
            ind2add=param.tracks(inds_alive(n_triples(row_ind,3),1)).bounds{end}{n_ind}(1,1);
            param.tracks(inds_alive(n_triples(row_ind,3),1)).bounds{end}{n_ind}=[param.tracks(inds_alive(n_triples(row_ind,3),1)).bounds{end}{n_ind}(1:ind2add,:) ; uint16(n_triples(row_ind,4:5)) ; param.tracks(inds_alive(n_triples(row_ind,3),1)).bounds{end}{n_ind}((ind2add+1):end,:)];
        end
    end
end

%     keyboard

% 2x2 pixel squares are problem to the endpoint detection, thus these
% points must be marked as separate end points
square_inds=find(conv2(single(~logical(segment_temp1)),[1 1;1 1],'same')==4);
square_inds=[square_inds;square_inds+1;square_inds+param.img_s(2)+4;square_inds+param.img_s(2)+4+1];

m_points=[segment_temp1(square_inds-param.img_s(2)-4-1) segment_temp1(square_inds-param.img_s(2)-4) segment_temp1(square_inds-param.img_s(2)-4+1)...
    segment_temp1(square_inds-1)  segment_temp1(square_inds) segment_temp1(square_inds+1)...
    segment_temp1(square_inds+param.img_s(2)+4-1) segment_temp1(square_inds+param.img_s(2)+4) segment_temp1(square_inds+param.img_s(2)+4+1)];
m_sort=sort(m_points,2);
m_diff=diff(sort(m_points,2),1,2)>0;

[x_ind y_ind]=find(m_diff');
n_pairs=[reshape(m_sort(sub2ind([size(m_diff,1) 9],y_ind,x_ind+1)),[2 size(m_sort,1)])' [floor(square_inds/(param.img_s(2)+4))+1-2 mod(square_inds-1,(param.img_s(2)+4))+1-2]];

% save edge pixels and neighbours to the tracklet struct
for e_ind=1:size(n_pairs,1)
    % first half
    n_ind=find(param.tracks(inds_alive(n_pairs(e_ind,1),1)).neighs{end}==inds_alive(n_pairs(e_ind,2),1));
    param.tracks(inds_alive(n_pairs(e_ind,1),1)).bounds{end}{n_ind}(1,1)=param.tracks(inds_alive(n_pairs(e_ind,1),1)).bounds{end}{n_ind}(1,1)+1;
    ind2add=param.tracks(inds_alive(n_pairs(e_ind,1),1)).bounds{end}{n_ind}(1,1);
    param.tracks(inds_alive(n_pairs(e_ind,1),1)).bounds{end}{n_ind}=[param.tracks(inds_alive(n_pairs(e_ind,1),1)).bounds{end}{n_ind}(1:ind2add,:) ; uint16(n_pairs(e_ind,3:4)) ; param.tracks(inds_alive(n_pairs(e_ind,1),1)).bounds{end}{n_ind}((ind2add+1):end,:)];
    
    % second half
    n_ind=find(param.tracks(inds_alive(n_pairs(e_ind,2),1)).neighs{end}==inds_alive(n_pairs(e_ind,1),1));
    param.tracks(inds_alive(n_pairs(e_ind,2),1)).bounds{end}{n_ind}(1,1)=param.tracks(inds_alive(n_pairs(e_ind,2),1)).bounds{end}{n_ind}(1,1)+1;
    ind2add=param.tracks(inds_alive(n_pairs(e_ind,2),1)).bounds{end}{n_ind}(1,1);
    param.tracks(inds_alive(n_pairs(e_ind,2),1)).bounds{end}{n_ind}=[param.tracks(inds_alive(n_pairs(e_ind,2),1)).bounds{end}{n_ind}(1:ind2add,:) ; uint16(n_pairs(e_ind,3:4)) ; param.tracks(inds_alive(n_pairs(e_ind,2),1)).bounds{end}{n_ind}((ind2add+1):end,:)];
end

inds_alive(end,:)=[];

end


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global param;

%% get cell locations from current selection
[x,y,c]=ginput(2);

[~,m_ind]=min((handles.selection(:,1)-x(1)).^2+(handles.selection(:,2)-y(1)).^2);
cell_i1=handles.selection(m_ind,3);
cell_i1_temp=cell_i1;

[~,m_ind]=min((handles.selection(:,1)-x(2)).^2+(handles.selection(:,2)-y(2)).^2);
cell_i2=handles.selection(m_ind,3);
cell_i2_temp=cell_i2;

% [cell_i1 cell_i2]

length_of_edge=zeros(handles.time,1,'uint8'); % [cell_id_selection cell_id time]

%% find length of junction
for time=1:length(length_of_edge)
    t_ind1=param.tracks(cell_i1).t==time;
    t_ind2=param.tracks(cell_i2).t==time;
    if sum(t_ind1)>0 & sum(t_ind2)>0        
        n_ind=param.tracks(cell_i1).neighs{t_ind1}==cell_i2;
        if sum(n_ind)>0
            length_of_edge(time)=size(param.tracks(cell_i1).bounds{t_ind1}{n_ind},1)-1+1;
        else
            length_of_edge(time)=1;
        end
    end
end

last_common_t=find(length_of_edge(1:(end-1))>0 & length_of_edge(2:(end))==0);
for time=(last_common_t+1):min(handles.time,last_common_t+15)
    
%     if time==109
%         keyboard
%     end
    
    if length(cell_i1)==1 & param.tracks(cell_i1).t(end)==(time-1) & param.tracks(cell_i1).daughters(1)>0
        cell_i1=param.tracks(cell_i1).daughters;
    end
    if length(cell_i2)==1 & param.tracks(cell_i2).t(end)==(time-1) & param.tracks(cell_i2).daughters(1)>0
        cell_i2=param.tracks(cell_i2).daughters;
    end    

    t_ind1=[];
    t_ind2=[];
    
    l=ones(2)*nan;
    
    for i_ind=1:length(cell_i1)
        ind=find(param.tracks(cell_i1(i_ind)).t==time);        
        if isempty(ind)
            t_ind1(i_ind)=0;
        else
            t_ind1(i_ind)=ind;
        end
    end
    
    for i_ind=1:length(cell_i2)
        ind=find(param.tracks(cell_i2(i_ind)).t==time);
        if isempty(ind)
            t_ind2(i_ind)=0;
        else
            t_ind2(i_ind)=ind;
        end
    end    
    
    for i_ind=find(t_ind1)
        for j_ind=find(t_ind2)
            n_ind=param.tracks(cell_i1(i_ind)).neighs{t_ind1(i_ind)}==cell_i2(j_ind);
            if sum(n_ind)>0
                l(i_ind,j_ind)=size(param.tracks(cell_i1(i_ind)).bounds{t_ind1(i_ind)}{n_ind},1)-1+1;
            else
                l(i_ind,j_ind)=1;
            end            
        end        
    end
    
%     time+14
%     l
    
%     length_of_edge(time)=sum(l(l>1));
%     if max(l(:))==1
%         length_of_edge(time)=1;
%     end    

    length_of_edge(time)=max(l(:));
    if isnan(length_of_edge(time))
        length_of_edge(time)=0;
    end    

end

% keyboard

cell_i1=cell_i1_temp;
cell_i2=cell_i2_temp;
first_common_t=find(length_of_edge(1:(end-1))==0 & length_of_edge(2:(end))>0);
for time=(first_common_t):-1:max(1,first_common_t-14)
    
%     if time==109
%         keyboard
%     end
    
    if param.tracks(cell_i1).t(1)==(time+1) & param.tracks(cell_i1).birth>0
        cell_i1=param.tracks(cell_i1).birth;
    end
    if param.tracks(cell_i2).t(1)==(time+1) & param.tracks(cell_i2).birth>0
        cell_i2=param.tracks(cell_i2).birth;
    end

    t_ind1=[];
    t_ind2=[];
    
    l=ones(2)*nan;
    
    for i_ind=1:length(cell_i1)
        ind=find(param.tracks(cell_i1(i_ind)).t==time);        
        if isempty(ind)
            t_ind1(i_ind)=0;
        else
            t_ind1(i_ind)=ind;
        end
    end
    
    for i_ind=1:length(cell_i2)
        ind=find(param.tracks(cell_i2(i_ind)).t==time);
        if isempty(ind)
            t_ind2(i_ind)=0;
        else
            t_ind2(i_ind)=ind;
        end
    end    
    
    for i_ind=find(t_ind1)
        for j_ind=find(t_ind2)
            n_ind=param.tracks(cell_i1(i_ind)).neighs{t_ind1(i_ind)}==cell_i2(j_ind);
            if sum(n_ind)>0
                l(i_ind,j_ind)=size(param.tracks(cell_i1(i_ind)).bounds{t_ind1(i_ind)}{n_ind},1)-1+1;
            else
                l(i_ind,j_ind)=1;
            end            
        end        
    end
    
%     time+14
%     l
    
%     length_of_edge(time)=sum(l(l>1));
%     if max(l(:))==1
%         length_of_edge(time)=1;
%     end    

    length_of_edge(time)=max(l(:));
    if isnan(length_of_edge(time))
        length_of_edge(time)=0;
    end    

end


%% plot junction length and...

non_smooth_curve=double(length_of_edge);
non_smooth_curve_times=(1:186)+14;

smooth_curve=smooth(non_smooth_curve,30,'lowess');

t_delta=2;

curve_smooth_deriv=[(smooth_curve(t_delta:end)-smooth_curve(1:(end-t_delta+1))).*double(antti_minfilt2(non_smooth_curve>0,[1 30])); 0];

curve_smooth_deriv_thres=(t_delta-1)/15;

shortening=curve_smooth_deriv<=-curve_smooth_deriv_thres;
lengthening=curve_smooth_deriv>=curve_smooth_deriv_thres;

figure(2)
clf
hold off

plot(non_smooth_curve_times,non_smooth_curve,non_smooth_curve_times,smooth_curve,'.r');

hold on;

plot(non_smooth_curve_times(shortening),ones(sum(shortening))*22,'.c');
plot(non_smooth_curve_times(lengthening),ones(sum(lengthening))*22,'.r');

ylim([0 44]);

% keyboard

% figure(3),hold off,plot(non_smooth_curve_times, antti_bfilter2(double(uint8(non_smooth_curve-1))/43,10,[3 1])*43+1,'.g',non_smooth_curve_times,non_smooth_curve);ylim([0 44]);
n=18;

non_smooth_curve_mod=antti_medfilt3(non_smooth_curve,[n 1]);

% non_smooth_curve_mod(non_smooth_curve_mod==1)=-22;
% non_smooth_curve_mod(non_smooth_curve_mod==0)=nan;

non_smooth_curve_mod(non_smooth_curve_mod==1)=-1;
non_smooth_curve_mod(non_smooth_curve==0)=nan;
non_smooth_curve_mod(non_smooth_curve_mod>1)=1;


disappearing_T1=conv(non_smooth_curve_mod,[-ones(n/2,1); ones(n/2,1)]/n*22,'same');

disappearing_T1(1:(n/2))=nan;
disappearing_T1((end-n/2+1):end)=nan;

% keyboard
[vals, X_dis]=findpeaks(disappearing_T1);
X_dis(vals<10)=[];
vals(vals<10)=[];

[vals, X_app]=findpeaks(-disappearing_T1);
X_app(vals<10)=[];
vals(vals<10)=[];

figure(3)
clf
hold on
plot(non_smooth_curve_times,disappearing_T1,'r',non_smooth_curve_times,non_smooth_curve,non_smooth_curve_times,non_smooth_curve_mod*22,'.b');
ylim([-22 44]);

plot(non_smooth_curve_times(X_dis), disappearing_T1(X_dis),'oc',non_smooth_curve_times(X_app), disappearing_T1(X_app),'om');

% plot((1:length(curve_smooth_deriv))+floor(t_delta/2),curve_smooth_deriv)

% keyboard

end

% new data set
% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% open grid file, keep asking, if cancelled
% this can be replaced with better question

while 1
    [filename, pathname, ~]=uigetfile('*.mat','Open grid file');
    if filename~=0
        break
    end
end
load([pathname filename],'timeStructureX','timeStructureY');

% find first image in image sequence    
while 1
    [filename, pathname, ~]=uigetfile('*.tif','Point first image from grid sequence');
    if filename~=0
        break
    end
end        
[~,name,~] = fileparts(filename);
underscore_inds=strfind(name,'_');

first_ind_in_img_seq=str2double(name(underscore_inds(end)+1:end));
image_path=[pathname name(1:underscore_inds(end))];

AllNames=dir([image_path '*']);
last_ind_in_img_seq=length(AllNames);

AnswerOut = inputdlg('Number of time points','Number of time points',1,{num2str(last_ind_in_img_seq)});
% user cancelled the dialog
if isempty(AnswerOut)
    return;
end
% check the index exceed the maximum
if str2num(AnswerOut{:})<=last_ind_in_img_seq
    last_ind_in_img_seq=str2num(AnswerOut{:});
end

handles.time_interval=[first_ind_in_img_seq last_ind_in_img_seq];

% show first image and the grid

figure(2)
clf
img_in=imread([image_path num2str(1+first_ind_in_img_seq-1,'%.04d') '.tif']);
imagesc(img_in);
colormap('gray')
hold on
line(timeStructureX{handles.time_interval(1)}',timeStructureY{handles.time_interval(1)}')
[x, y]=ginput(1);

[~,time_struct_id]=min((mean(timeStructureX{handles.time_interval(1)},2)-x).^2+(mean(timeStructureY{handles.time_interval(1)},2)-y).^2);

close(2)

boundary=50;


%%

%     keyboard

for time=handles.time_interval(1):handles.time_interval(2)
    handles.min_max_mean(time-handles.time_interval(1)+1,1:6)=[min(timeStructureX{time}(time_struct_id,:)) max(timeStructureX{time}(time_struct_id,:)) min(timeStructureY{time}(time_struct_id,:)) max(timeStructureY{time}(time_struct_id,:)) mean(timeStructureX{time}(time_struct_id,:)) mean(timeStructureY{time}(time_struct_id,:))];
end

handles.max_x_span=ceil((ceil(max(handles.min_max_mean(:,2)-handles.min_max_mean(:,1))+boundary*2)/2))*2;
handles.max_y_span=ceil((ceil(max(handles.min_max_mean(:,4)-handles.min_max_mean(:,3))+boundary*2)/2))*2;

handles.o_imgs=zeros([handles.max_y_span+1 handles.max_x_span+1 handles.time_interval(2)-handles.time_interval(1)+1],'uint8');
handles.s_imgs=zeros(size(handles.o_imgs),'uint8');
handles.s_imgs_independent=zeros(size(handles.o_imgs),'uint8');   

for time=handles.time_interval(1):handles.time_interval(2)
%         handles.box(time-handles.time_interval(1)+1,1:4)=[round(handles.min_max_mean(time-handles.time_interval(1)+1,5))-handles.max_x_span/2 round(handles.min_max_mean(time-handles.time_interval(1)+1,5))+handles.max_x_span/2 round(handles.min_max_mean(time-handles.time_interval(1)+1,6))-handles.max_y_span/2 round(handles.min_max_mean(time-handles.time_interval(1)+1,6))+handles.max_y_span/2];
    handles.box(time-handles.time_interval(1)+1,1:4)=[round(mean(handles.min_max_mean(time-handles.time_interval(1)+1,5)))-handles.max_x_span/2 round(mean(handles.min_max_mean(time-handles.time_interval(1)+1,5)))+handles.max_x_span/2 round(mean(handles.min_max_mean(time-handles.time_interval(1)+1,6)))-handles.max_y_span/2 round(mean(handles.min_max_mean(time-handles.time_interval(1)+1,6)))+handles.max_y_span/2];
    if sum(isnan(handles.box(time-handles.time_interval(1)+1,1:4)))>0
        break;
    end
    pixels=cell(1,2);
%         pixels{1}=handles.box(time-handles.time_interval(1)+1,3:4);
%         pixels{2}=handles.box(time-handles.time_interval(1)+1,1:2);   

    pixels{1}=min(max(handles.box(time-handles.time_interval(1)+1,3:4),1),size(img_in,1));
    pixels{2}=min(max(handles.box(time-handles.time_interval(1)+1,1:2),1),size(img_in,2));        

    if handles.box(time-handles.time_interval(1)+1,3)<1
        fy=-handles.box(time-handles.time_interval(1)+1,3)+2;
    else
        fy=1;
    end

    if handles.box(time-handles.time_interval(1)+1,4)>size(img_in,1)
        ly=diff(handles.box(time-handles.time_interval(1)+1,3:4))+1-(handles.box(time-handles.time_interval(1)+1,4)-size(img_in,1));
    else
        ly=diff(handles.box(time-handles.time_interval(1)+1,3:4))+1;
    end        

%         handles.o_imgs(:,:,time-handles.time_interval(1)+1)=imread([image_path num2str(time+first_ind_in_img_seq-1,'%.04d') '.tif'],'pixelregion',pixels);    
    handles.o_imgs(fy:ly,:,time-handles.time_interval(1)+1)=imread([image_path num2str(time,'%.04d') '.tif'],'pixelregion',pixels);
end

handles.f_imgs=zeros(size(handles.o_imgs),'uint8');
imwrite(handles.o_imgs(:,:,1),'temp.tif');
for i_ind=2:size(handles.o_imgs,3)
    imwrite(handles.o_imgs(:,:,i_ind),'temp.tif','writemode','append');
end
%     antti_evaluate_imagej_script_windows('Y:\Antti\programs\semi_automatic_tracker\bandpass_filter_z-stack.ijm', 'Z:\mDrives\raid0\DSLM_WorkInProgress\Exp186\antti\0087_edit_time_independent_and_track\temp.tif','Z:\mDrives\raid0\DSLM_WorkInProgress\Exp186\antti\0087_edit_time_independent_and_track\filt_temp.tif');
antti_evaluate_imagej_script_windows('Y:\Antti\programs\semi_automatic_tracker\bandpass_filter_z-stack.ijm',[pwd filesep 'temp.tif'],[pwd filesep 'filt_temp.tif']);
for i_ind=1:size(handles.o_imgs,3)
    handles.f_imgs(:,:,i_ind)=imread('filt_temp.tif','index',i_ind);
end
delete temp.tif filt_temp.tif

handles.time=[];    
handles.cur_img=imagesc(handles.o_imgs(:,:,1),'parent',handles.axis1);
set(handles.slider2,'Min',handles.time_interval(1),'Max',handles.time_interval(2),'Value',handles.time_interval(1),'SliderStep',[1/(handles.time_interval(2)-handles.time_interval(1)) 10/(handles.time_interval(2)-handles.time_interval(1))]);    
set(handles.text1,'string',num2str(handles.time_interval(1)));

handles.selection=[];
handles.centroid_scatter=[];
handles.centroid_tracks_lines=[];
handles.centroid_neighbour_lines=[];
handles.centroid_polygons=[];

handles.c_ind_max=16;
handles.c_map=jet(handles.c_ind_max);

handles.correspondence=[];

handles.inds_alive=[];
handles.labeled_independent_image=[];
handles.cents=[];

handles.projected_points=[];
handles.projection_scatter=[];

colormap(handles.axis1,'gray');
hold(handles.axis1,'on');
axis(handles.axis1,'off','equal');    
set(handles.cur_img,'ButtonDownFcn',@position_and_button);
set(handles.pushbutton12,'enable','off');
guidata(hObject, handles);

end

% load data set
% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname, index]=uigetfile('*.mat');

% return if user cancelled
if index==0
    return
end

load([pathname filename]);
global param;

handles.time_interval=handles1.time_interval;

handles.min_max_mean=handles1.min_max_mean;
handles.max_x_span=handles1.max_x_span;
handles.max_y_span=handles1.max_y_span;
handles.o_imgs=handles1.o_imgs;
handles.s_imgs=handles1.s_imgs;
handles.s_imgs_independent=handles1.s_imgs_independent;
handles.box=handles1.box;
handles.f_imgs=handles1.f_imgs;
handles.time=handles1.time;
handles.selection=handles1.selection;
handles.centroid_scatter=[];
handles.centroid_tracks_lines=[];
handles.centroid_neighbour_lines=[];
handles.centroid_polygons=[];
handles.c_ind_max=handles1.c_ind_max;
handles.c_map=handles1.c_map;
handles.correspondence=handles1.correspondence;
handles.slider2_value=round(handles1.slider2_value);
handles.projected_points=handles1.projected_points;
handles.projection_scatter=[];

handles.inds_alive=handles1.inds_alive;    
handles.labeled_independent_image=handles1.labeled_independent_image;    
handles.cents=handles1.cents;    

handles.cur_img=imagesc(handles.o_imgs(:,:,handles.slider2_value-handles.time_interval(1)+1),'parent',handles.axis1);    

colormap(handles.axis1,'gray');
hold(handles.axis1,'on');
axis(handles.axis1,'off','equal');    

handles.centroid_scatter=draw_selection(handles);
handles.projection_scatter=draw_projection(handles);
handles.centroid_tracks_lines=draw_tracks(handles);
handles.centroid_neighbour_lines=draw_neighbours(handles);
handles.centroid_polygons=draw_polygons(handles);

guidata(hObject, handles);
set(handles.slider2,'Min',handles.time_interval(1),'Max',handles.time_interval(2),'Value',handles.slider2_value,'SliderStep',[1/(handles.time_interval(2)-handles.time_interval(1)) 10/(handles.time_interval(2)-handles.time_interval(1))]);
set(handles.text1,'string',num2str(handles.slider2_value));

draw_image(handles);
set(handles.cur_img,'ButtonDownFcn',@position_and_button);


end


% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox10
handles.centroid_tracks_lines=draw_tracks(handles);
guidata(handles.figure1, handles);
handles.centroid_neighbour_lines=draw_neighbours(handles);
handles.centroid_polygons=draw_polygons(handles);
guidata(handles.figure1, handles);
end


% --- Executes on button press in checkbox11.
function checkbox11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox11
handles.centroid_tracks_lines=draw_tracks(handles);
guidata(handles.figure1, handles);
handles.centroid_neighbour_lines=draw_neighbours(handles);
handles.centroid_polygons=draw_polygons(handles);
guidata(handles.figure1, handles);
end