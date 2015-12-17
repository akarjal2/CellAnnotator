function s_points_to_rm=merge_first_time_point_cells(m_cell_ind,dd_cell_ind)

% keyboard

global param;

mother=param.tracks(m_cell_ind);
dead_daughter=param.tracks(dd_cell_ind);

m_cell=mother;

% %                t: 1
% %                A: 377
% %             cent: [12 140]
% %          ellipse: [24 21 -41]
% %            perim: 73
% %            birth: -1
% %            death: -11
% %        daughters: [1490 1491]
% %           neighs: {[3 4 34 41 48 0]}
% %           bounds: {{1x6 cell}}
% %      num_overlap: 0
% %     num_underlap: 0

t_inds_to_add=dead_daughter.t;

m_cell.A(t_inds_to_add)=m_cell.A(t_inds_to_add)+dead_daughter.A;

m_cell.cent(t_inds_to_add,:)=(repmat(double(m_cell.A(t_inds_to_add)'),[1,2]).*double(m_cell.cent(t_inds_to_add,:))+repmat(double(dead_daughter.A'),[1 2]).*double(dead_daughter.cent))./(repmat(double(dead_daughter.A')+double(m_cell.A(t_inds_to_add)'),[1 2]));

% this need to be done properly, now ellipses are really not combined
% ellipses...

m_cell.perim(t_inds_to_add)=m_cell.perim(t_inds_to_add)+dead_daughter.perim;

for i_ind=1:length(t_inds_to_add)
% remove mutual edge    
    
    temp_s_points=m_cell.bounds{t_inds_to_add(i_ind)}{m_cell.neighs{t_inds_to_add(i_ind)}==dd_cell_ind};
    
    s_points_to_rm{t_inds_to_add(i_ind)}=temp_s_points((temp_s_points(1,1)+1+1):end,:);

    m_cell.bounds{t_inds_to_add(i_ind)}(m_cell.neighs{t_inds_to_add(i_ind)}==dd_cell_ind)=[];    
    m_cell.neighs{t_inds_to_add(i_ind)}(m_cell.neighs{t_inds_to_add(i_ind)}==dd_cell_ind)=[];    
    
    dead_daughter.bounds{i_ind}(dead_daughter.neighs{i_ind}==m_cell_ind)=[];
    dead_daughter.neighs{i_ind}(dead_daughter.neighs{i_ind}==m_cell_ind)=[];    
    
    for j_ind=1:length(m_cell.neighs{t_inds_to_add(i_ind)})
% find if daughters had common neighbours and modify accordingly
        n_ind=find(dead_daughter.neighs{i_ind}==m_cell.neighs{t_inds_to_add(i_ind)}(j_ind));
        
%         if para==878
%             keyboard
%         end
        
        if ~isempty(n_ind)
            dead_points=dead_daughter.bounds{i_ind}{n_ind};
            m_cell.bounds{t_inds_to_add(i_ind)}{j_ind}=[m_cell.bounds{t_inds_to_add(i_ind)}{j_ind};dead_points(2:end,:)];            
                                    
            if m_cell.neighs{t_inds_to_add(i_ind)}(j_ind)~=0
                t_ind=find(param.tracks(m_cell.neighs{t_inds_to_add(i_ind)}(j_ind)).t==m_cell.t(t_inds_to_add(i_ind)));
                ind_to_mod=find(param.tracks(m_cell.neighs{t_inds_to_add(i_ind)}(j_ind)).neighs{t_ind}==m_cell_ind);

                param.tracks(m_cell.neighs{t_inds_to_add(i_ind)}(j_ind)).neighs{t_ind}(ind_to_mod)=m_cell_ind;
                param.tracks(m_cell.neighs{t_inds_to_add(i_ind)}(j_ind)).bounds{t_ind}{ind_to_mod}=[m_cell.bounds{t_inds_to_add(i_ind)}{j_ind};dead_points(2:end,:)];

                ind_to_mod=find(param.tracks(m_cell.neighs{t_inds_to_add(i_ind)}(j_ind)).neighs{t_ind}==dd_cell_ind);
                param.tracks(m_cell.neighs{t_inds_to_add(i_ind)}(j_ind)).neighs{t_ind}(ind_to_mod)=[];
                param.tracks(m_cell.neighs{t_inds_to_add(i_ind)}(j_ind)).bounds{t_ind}(ind_to_mod)=[];
            end
            
            dead_daughter.neighs{i_ind}(n_ind)=[];
            dead_daughter.bounds{i_ind}(n_ind)=[];
                        
        elseif m_cell.neighs{t_inds_to_add(i_ind)}(j_ind)~=0
            t_ind=find(param.tracks(m_cell.neighs{t_inds_to_add(i_ind)}(j_ind)).t==m_cell.t(t_inds_to_add(i_ind)));
            ind_to_mod=find(param.tracks(m_cell.neighs{t_inds_to_add(i_ind)}(j_ind)).neighs{t_ind}==m_cell_ind);
            param.tracks(m_cell.neighs{t_inds_to_add(i_ind)}(j_ind)).neighs{t_ind}(ind_to_mod)=m_cell_ind;                
        end
    end

    for j_ind=1:length(dead_daughter.neighs{i_ind})   
        
%         if dead_daughter.neighs{i_ind}(j_ind)==878
%             keyboard
%         end
        
        if dead_daughter.neighs{i_ind}(j_ind)~=0
            t_ind=find(param.tracks(dead_daughter.neighs{(i_ind)}(j_ind)).t==dead_daughter.t((i_ind)));
            ind_to_mod=find(param.tracks(dead_daughter.neighs{(i_ind)}(j_ind)).neighs{t_ind}==dd_cell_ind);
            param.tracks(dead_daughter.neighs{(i_ind)}(j_ind)).neighs{t_ind}(ind_to_mod)=m_cell_ind;                
        end
    end    
    
    m_cell.neighs{t_inds_to_add(i_ind)}=[m_cell.neighs{t_inds_to_add(i_ind)} dead_daughter.neighs{i_ind}];
    m_cell.bounds{t_inds_to_add(i_ind)}=[m_cell.bounds{t_inds_to_add(i_ind)} dead_daughter.bounds{i_ind}];
    
end

% m_cell.num_underlap(t_inds_to_add)=1;

param.tracks(m_cell_ind)=m_cell;

param.tracks(dd_cell_ind).t=[];
param.tracks(dd_cell_ind).A=[];
param.tracks(dd_cell_ind).cent=[];
param.tracks(dd_cell_ind).ellipse=[];
param.tracks(dd_cell_ind).perim=[];
param.tracks(dd_cell_ind).birth=[];
param.tracks(dd_cell_ind).death=[];
param.tracks(dd_cell_ind).daughters=[];
param.tracks(dd_cell_ind).neighs=[];
param.tracks(dd_cell_ind).bounds=[];
% param.tracks(dd_cell_ind).num_overlap=[];
% param.tracks(dd_cell_ind).num_underlap=[];

