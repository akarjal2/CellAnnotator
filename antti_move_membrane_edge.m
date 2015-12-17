function [did_move_happen,add_segment_line,rm_segment_line]=antti_move_membrane_edge(x,y,candidate,rm_ind_cell_i,f_img_m)

% keyboard

global param

cell_i=candidate(1);
cell_i_t=candidate(2);

cell_i2=param.tracks(cell_i).neighs{cell_i_t}(rm_ind_cell_i);
cell_i2_t=find(param.tracks(cell_i2).t==param.tracks(cell_i).t(cell_i_t));
rm_ind_cell_i2=find(param.tracks(cell_i2).neighs{cell_i2_t}==cell_i);

did_move_happen=1;

add_segment_line=[];
rm_segment_line=[];

% find points on outer boundary
e_points=zeros(500,2);
e_points_ind=1;
% first cell
for b_ind=1:length(param.tracks(cell_i).bounds{cell_i_t})
    if b_ind~=rm_ind_cell_i                
        e_points(e_points_ind:(e_points_ind+size(param.tracks(cell_i).bounds{cell_i_t}{b_ind},1)-2),:)=param.tracks(cell_i).bounds{cell_i_t}{b_ind}(2:end,:);
        e_points_ind=e_points_ind+size(param.tracks(cell_i).bounds{cell_i_t}{b_ind},1)-1;
    else
        points_to_add=param.tracks(cell_i).bounds{cell_i_t}{b_ind};
        rm_segment_line=points_to_add((points_to_add(1,1)+2):end,:);
        points_to_add([1 (points_to_add(1,1)+2):end],:)=[];
        e_points(e_points_ind:(e_points_ind+size(points_to_add,1)-1),:)=points_to_add;
        e_points_ind=e_points_ind+size(points_to_add,1);        
    end        
end
% second cell
for b_ind=1:length(param.tracks(cell_i2).bounds{cell_i2_t})
    if b_ind~=rm_ind_cell_i2                
        e_points(e_points_ind:(e_points_ind+size(param.tracks(cell_i2).bounds{cell_i2_t}{b_ind},1)-2),:)=param.tracks(cell_i2).bounds{cell_i2_t}{b_ind}(2:end,:);
        e_points_ind=e_points_ind+size(param.tracks(cell_i2).bounds{cell_i2_t}{b_ind},1)-1;
    else
        points_to_add=param.tracks(cell_i2).bounds{cell_i2_t}{b_ind};
        points_to_add([1 (points_to_add(1,1)+2):end],:)=[];
        e_points(e_points_ind:(e_points_ind+size(points_to_add,1)-1),:)=points_to_add;
        e_points_ind=e_points_ind+size(points_to_add,1);        
    end        
end
e_points(e_points_ind:end,:)=[];

min_e_points=min(e_points);
max_e_points=max(e_points);
div_w_s=fliplr(max_e_points-min_e_points+1);

mask=logical(zeros(div_w_s+[2 2]));
mask(sub2ind(div_w_s+[2 2],e_points(:,2)-min_e_points(2)+2,e_points(:,1)-min_e_points(1)+2))=1;
mask1=imfill(mask,'holes');
mask2=mask1-mask;
[y_c, x_c]=find(mask2);

div_img=zeros(div_w_s+[2 2]);
div_img(sub2ind(div_w_s+[2 2],y_c,x_c))=f_img_m(sub2ind(fliplr(param.img_s),y_c-2+min_e_points(2),x_c-2+min_e_points(1)));
div_img(div_img==(2^8-1))=2^8-2;

div_img(find(mask))=2^8-1;

d_seed1=uint16([round(x(1)) round(y(1))]);
d_seed2=uint16([round(x(2)) round(y(2))]);

% keyboard

daughter_mask=zeros(size(mask2));

radius=sqrt(sum((double(d_seed1)-double(d_seed2)).^2))/3;

daughter_mask=antti_ellipseMatrix(double(d_seed1(2)-min_e_points(2)+2),double(d_seed1(1)-min_e_points(1)+2),radius,radius,0,daughter_mask,1);
daughter_mask=antti_ellipseMatrix(double(d_seed2(2)-min_e_points(2)+2),double(d_seed2(1)-min_e_points(1)+2),radius,radius,0,daughter_mask,1);


daughter_mask=daughter_mask&mask2;

div_img(find(daughter_mask))=0;
% segment
div_img_mod=imhmin(div_img,2^7);
div_segment=watershed(div_img_mod);

% collect results
div_info=regionprops(div_segment,'centroid','area','perimeter','orientation','majoraxislength','minoraxislength');

% keyboard

if (d_seed1(2)-min_e_points(2)+2)<1 | (d_seed1(1)-min_e_points(1)+2)<1 | (d_seed1(2)-min_e_points(2)+2)>size(div_segment,1) | (d_seed1(1)-min_e_points(1)+2)>size(div_segment,2) | (d_seed2(2)-min_e_points(2)+2)<1 | (d_seed2(1)-min_e_points(1)+2)<1 | (d_seed2(2)-min_e_points(2)+2)>size(div_segment,1) | (d_seed2(1)-min_e_points(1)+2)>size(div_segment,2)
    did_move_happen=0;
else
    d_ind1=div_segment(d_seed1(2)-min_e_points(2)+2,d_seed1(1)-min_e_points(1)+2);
    d_ind2=div_segment(d_seed2(2)-min_e_points(2)+2,d_seed2(1)-min_e_points(1)+2);
    
    if d_ind1==1 | d_ind1==0 | d_ind2==1 | d_ind2==0 | d_ind1==d_ind2 | length(div_info)>3
        did_move_happen=0;
    end
end

% 190     collect segment information from segmented image
% replace mother segment with daughter segments

if did_move_happen==0
    return;
end

% keyboard

% create daughters

param.tracks(end+1)=struct('t',uint16(1),'A',uint32(div_info(d_ind1).Area),'cent',uint16(div_info(d_ind1).Centroid+min_e_points-2),'ellipse',int16([div_info(d_ind1).MajorAxisLength div_info(d_ind1).MinorAxisLength div_info(d_ind1).Orientation]),'perim',uint32(div_info(d_ind1).Perimeter),'birth',int32(cell_i),'death',int16(0),'daughters',uint32([0 0]),'neighs',cell(1),'bounds',cell(1));

param.tracks(end+1)=struct('t',uint16(1),'A',uint32(div_info(d_ind2).Area),'cent',uint16(div_info(d_ind2).Centroid+min_e_points-2),'ellipse',int16([div_info(d_ind2).MajorAxisLength div_info(d_ind2).MinorAxisLength div_info(d_ind2).Orientation]),'perim',uint32(div_info(d_ind2).Perimeter),'birth',int32(cell_i),'death',int16(0),'daughters',uint32([0 0]),'neighs',cell(1),'bounds',cell(1));

param.tracks(end-1).bounds=cell(1);
param.tracks(end-1).bounds{end}=cell(0);
param.tracks(end-1).neighs=cell(1);

param.tracks(end).bounds=cell(1);
param.tracks(end).bounds{end}=cell(0);
param.tracks(end).neighs=cell(1);

first_new_cell_ind=cell_i;
second_new_cell_ind=cell_i2;

% update neighbours outer membrane of "old" mother segment
for mom_neigh=1:length(param.tracks(cell_i).bounds{cell_i_t})
    num_boundary=param.tracks(cell_i).bounds{cell_i_t}{mom_neigh}(1,1);
    
    ind1_found=zeros(num_boundary,1);
    ind2_found=zeros(num_boundary,1);
    for end_ind=1:num_boundary
        b_point=double(param.tracks(cell_i).bounds{cell_i_t}{mom_neigh}(1+end_ind,:))-min_e_points+2;
        b_point_neigh=div_segment(b_point(2)-1:b_point(2)+1,b_point(1)-1:b_point(1)+1);
        ind1_found(end_ind)=~isempty(find(b_point_neigh==d_ind1));
        ind2_found(end_ind)=~isempty(find(b_point_neigh==d_ind2));
    end
    
    % if second daughter cell does not have any common boundary endpoints with
    % current segment, must whole corresponding boundary belong to the first
    % daughter cell
    if isempty(find(ind2_found))
        neighbour_track_i=param.tracks(cell_i).neighs{cell_i_t}(mom_neigh);
        if neighbour_track_i>0
            neighbour_bound_i=find(param.tracks(neighbour_track_i).neighs{end}==cell_i);
            param.tracks(neighbour_track_i).neighs{end}(neighbour_bound_i)=first_new_cell_ind;
        end
        param.tracks(length(param.tracks)-1).neighs{end}(end+1)=neighbour_track_i;
        param.tracks(length(param.tracks)-1).bounds{end}{end+1}=param.tracks(cell_i).bounds{cell_i_t}{mom_neigh};
        
        % if first daughter cell does not have any common boundary endpoints with
        % current segment, must whole corresponding boundary belong to the second
        % daughter cell
    elseif isempty(find(ind1_found))
        neighbour_track_i=param.tracks(cell_i).neighs{cell_i_t}(mom_neigh);
        if neighbour_track_i>0
            neighbour_bound_i=find(param.tracks(neighbour_track_i).neighs{end}==cell_i);
            param.tracks(neighbour_track_i).neighs{end}(neighbour_bound_i)=second_new_cell_ind;
        end
        param.tracks(length(param.tracks)).neighs{end}(end+1)=neighbour_track_i;
        param.tracks(length(param.tracks)).bounds{end}{end+1}=param.tracks(cell_i).bounds{cell_i_t}{mom_neigh};
        
        % boundary need to be split into two
    else
        bound_points=zeros(size(param.tracks(cell_i).bounds{cell_i_t}{mom_neigh},1)-1,1);
        for b_p_ind=2:size(param.tracks(cell_i).bounds{cell_i_t}{mom_neigh},1)
            b_point=double(param.tracks(cell_i).bounds{cell_i_t}{mom_neigh}(b_p_ind,:))-min_e_points+2;
            b_point_neigh=div_segment(b_point(2)-1:b_point(2)+1,b_point(1)-1:b_point(1)+1);
            ind1_ok=~isempty(find(b_point_neigh==d_ind1));
            ind2_ok=~isempty(find(b_point_neigh==d_ind2));
            bound_points(b_p_ind-1)=ind1_ok+ind2_ok*2;
        end
        
        % form first half of membrane segment
        b_segment=param.tracks(cell_i).bounds{cell_i_t}{mom_neigh}([find(ind1_found)+1 ; find(bound_points==3)+1],:);
        b_segment=[size(b_segment,1) 0 ; b_segment ; param.tracks(cell_i).bounds{cell_i_t}{mom_neigh}(find((bound_points==1) | (bound_points==0))+1,:)];
        
        neighbour_track_i=param.tracks(cell_i).neighs{cell_i_t}(mom_neigh);
        if neighbour_track_i>0
            neighbour_bound_i=find(param.tracks(neighbour_track_i).neighs{end}==cell_i);
            param.tracks(neighbour_track_i).neighs{end}(neighbour_bound_i)=first_new_cell_ind;
            param.tracks(neighbour_track_i).bounds{end}{neighbour_bound_i}=uint16(b_segment);
        end
        
        param.tracks(length(param.tracks)-1).bounds{end}{end+1}=uint16(b_segment);
        param.tracks(length(param.tracks)-1).neighs{end}(end+1)=neighbour_track_i;
        
        % form second half of membrane segment
        b_segment=param.tracks(cell_i).bounds{cell_i_t}{mom_neigh}([find(ind2_found)+1 ; find(bound_points==3)+1],:);
        b_segment=[size(b_segment,1) 0 ; b_segment ; param.tracks(cell_i).bounds{cell_i_t}{mom_neigh}(find((bound_points==2)  | (bound_points==0))+1,:)];
        
        neighbour_track_i=param.tracks(cell_i).neighs{cell_i_t}(mom_neigh);
        if neighbour_track_i>0
            param.tracks(neighbour_track_i).neighs{end}(end+1)=second_new_cell_ind;
            param.tracks(neighbour_track_i).bounds{end}{end+1}=uint16(b_segment);
        end
        
        param.tracks(length(param.tracks)).bounds{end}{end+1}=uint16(b_segment);
        param.tracks(length(param.tracks)).neighs{end}(end+1)=neighbour_track_i;
    end
end
% find and create membrane between two new daughter segments

b_segment=[0 0];
[ind_y, ind_x]=find(div_segment==0);
for b_ind=1:length(ind_y)
    b_neigh=div_segment(ind_y(b_ind)-1:ind_y(b_ind)+1,ind_x(b_ind)-1:ind_x(b_ind)+1);
    ind1_ok=~isempty(find(b_neigh==d_ind1));
    ind2_ok=~isempty(find(b_neigh==d_ind2));
    indb_ok=~isempty(find(b_neigh==1));
    if ~indb_ok
        b_segment(end+1,:)=[ind_x(b_ind) ind_y(b_ind)];
    elseif ind2_ok & ind1_ok
        b_segment(1,1)=b_segment(1,1)+1;
        b_segment=[b_segment(1:b_segment(1,1),:) ; ind_x(b_ind) ind_y(b_ind) ; b_segment((b_segment(1,1)+1):end,:)];
    end
end
b_segment(2:end,:)=b_segment(2:end,:)+repmat(min_e_points-2,[size(b_segment,1)-1 1]);

param.tracks(length(param.tracks)-1).bounds{end}{end+1}=uint16(b_segment);
param.tracks(length(param.tracks)-1).neighs{end}(end+1)=second_new_cell_ind;

param.tracks(length(param.tracks)).bounds{end}{end+1}=uint16(b_segment);
param.tracks(length(param.tracks)).neighs{end}(end+1)=first_new_cell_ind;

b_segment(find(b_segment(:,1)<1 | b_segment(:,1)>param.img_s(1) | b_segment(:,2)<1 | b_segment(:,2)>param.img_s(2)),:)=[];

add_segment_line=b_segment(:,2)+(b_segment(:,1)-1)*param.img_s(2);

% do not kill mother"s"    
param.tracks(cell_i).A(end)=param.tracks(end-1).A(end);
param.tracks(cell_i).cent(end,:)=param.tracks(end-1).cent(end,:);
param.tracks(cell_i).ellipse(end,:)=param.tracks(end-1).ellipse(end,:);
param.tracks(cell_i).perim(end)=param.tracks(end-1).perim(end);        
param.tracks(cell_i).neighs(end)=param.tracks(end-1).neighs(end);
param.tracks(cell_i).bounds(end)=param.tracks(end-1).bounds(end);

param.tracks(cell_i2).A(end)=param.tracks(end).A(end);
param.tracks(cell_i2).cent(end,:)=param.tracks(end).cent(end,:);
param.tracks(cell_i2).ellipse(end,:)=param.tracks(end).ellipse(end,:);
param.tracks(cell_i2).perim(end)=param.tracks(end).perim(end);        
param.tracks(cell_i2).neighs(end)=param.tracks(end).neighs(end);
param.tracks(cell_i2).bounds(end)=param.tracks(end).bounds(end);

param.tracks((end-1):end)=[];    


