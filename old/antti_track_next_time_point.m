function [s_img]=antti_track_next_time_point(f_img,time,projected_points)

global param

bw_mask=zeros(fliplr(param.img_s));
bw_mask(sub2ind(fliplr(param.img_s),double(projected_points(:,3)),double(projected_points(:,2))))=1;
segment_temp=(watershed(imimposemin(f_img,bw_mask)));

inds_alive=projected_points;

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

