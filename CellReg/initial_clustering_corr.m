function [cell_to_index_map, correlation_map,all_neighbor_correlations,all_neighbor_distances,all_assigned_correlations,non_assigned_correlations,num_candidates,is_in_ROI]=initial_clustering_corr(distance_thresh,correlation_thresh,all_filters_corrected,all_centroids_corrected,num_sessions,varargin)
% Initial clustering by correlation threshold:

if ~isempty(varargin)
    cells_in_ROI=varargin{1};
end
initial_cell_num=size(all_filters_corrected{1},1);
cell_to_index_map=zeros(initial_cell_num,num_sessions);
cell_to_index_map(:,1)=1:initial_cell_num;
correlation_map=zeros(initial_cell_num,num_sessions);    
all_sessions_filters=all_filters_corrected{1};
all_sessions_centroids=all_centroids_corrected{1};
count=0;


duplicate_match_count=0;
all_neighbor_correlations=zeros(1,num_sessions^2*initial_cell_num);
all_assigned_correlations=zeros(1,num_sessions^2*initial_cell_num);
non_assigned_correlations=zeros(1,num_sessions^2*initial_cell_num);
all_neighbor_distances=zeros(1,num_sessions^2*initial_cell_num);
is_in_ROI=zeros(1,num_sessions^2*initial_cell_num);

neighbor_count=0;
assigned_count=0;
non_assigned_count=0;
num_candidates=0;
h=waitbar(0,'Registering cells','Units', 'normalized', 'Position',[0.4 0.5 0.2 0.07]);
for n=2:num_sessions;
    waitbar((n-1)/(num_sessions-1),h,['Registering cells in session number ' num2str(n) '/' num2str(num_sessions)])
    new_filters=all_filters_corrected{n};
    new_centroids=all_centroids_corrected{n};
    h2=waitbar(0,'Registering cells','Units', 'normalized', 'Position',[0.4 0.4 0.2 0.07]);
    for k=1:size(new_filters,1)
        waitbar((k-1)/size(new_filters,1),h2,['Registering cell number ' num2str(k) '/' num2str(size(new_filters,1))])
        is_assigned=0;
        new_filter=squeeze(new_filters(k,:,:));
        centroid=repmat(new_centroids(k,:),size(all_sessions_centroids,1),1);
        if ~isempty(varargin)
            is_in_ROI_temp=cells_in_ROI{n}(k);
        end
        distance_vec=sqrt(sum((centroid-all_sessions_centroids).^2,2));
        filters_to_check=find(distance_vec<distance_thresh);
        if ~isempty(filters_to_check)
            corr_vec=zeros(1,length(filters_to_check));
            num_candidates=num_candidates+length(filters_to_check);
            for m=1:length(filters_to_check)
                suspected_filter=squeeze(all_sessions_filters(filters_to_check(m),:,:));
                if sum(sum(suspected_filter))==0 || sum(sum(new_filter))==0
                    corr_vec(m)=0;
                else
                    neighbor_count=neighbor_count+1;                  
                    corr_vec(m)=corr2(suspected_filter,new_filter);
                    all_neighbor_correlations(neighbor_count)=corr_vec(m);
                    all_neighbor_distances(neighbor_count)=distance_vec(filters_to_check(m));
                    if ~isempty(varargin)                        
                        is_in_ROI(neighbor_count)=is_in_ROI_temp;
                    end
                end
            end
            [highest_corr,highest_corr_ind]=max(corr_vec);
            if highest_corr<correlation_thresh
                count=count+1;
                all_sessions_filters(initial_cell_num+count,:,:)=new_filter;
                all_sessions_centroids(initial_cell_num+count,:)=new_centroids(k,:);
                cell_to_index_map(initial_cell_num+count,:)=zeros(1,num_sessions);
                cell_to_index_map(initial_cell_num+count,n)=k;
                correlation_map(initial_cell_num+count,:)=zeros(1,num_sessions);
            else
                index=filters_to_check(highest_corr_ind);
                if cell_to_index_map(index,n)==0
                    cell_to_index_map(index,n)=k;               
                    correlation_map(index,n)=highest_corr;
                    assigned_count=assigned_count+1;
                    is_assigned=1;
                    all_assigned_correlations(1,assigned_count)=highest_corr;
                else
                    duplicate_match_count=duplicate_match_count+1;
                    if highest_corr>correlation_map(index,n) % switch between cells
                        count=count+1;
                        switch_cell=cell_to_index_map(index,n);
                        switch_filter=squeeze(all_filters_corrected{n}(switch_cell,:,:));
                        switch_centroid=(all_centroids_corrected{n}(switch_cell,:));
                        all_sessions_filters(initial_cell_num+count,:,:)=switch_filter;
                        all_sessions_centroids(initial_cell_num+count,:)=switch_centroid;
                        cell_to_index_map(initial_cell_num+count,n)=switch_cell;
                        correlation_map(initial_cell_num+count,:)=zeros(1,num_sessions);
                        cell_to_index_map(index,n)=k;
                        correlation_map(index,n)=highest_corr;
                        assigned_count=assigned_count+1;
                        is_assigned=1;
                        all_assigned_correlations(1,assigned_count)=highest_corr;               
                    else
                        count=count+1;
                        all_sessions_filters(initial_cell_num+count,:,:)=new_filter;
                        all_sessions_centroids(initial_cell_num+count,:)=new_centroids(k,:);
                        cell_to_index_map(initial_cell_num+count,:)=zeros(1,num_sessions);
                        cell_to_index_map(initial_cell_num+count,n)=k;
                        correlation_map(initial_cell_num+count,:)=zeros(1,num_sessions);
                    end
                end
            end
            if is_assigned==0
                non_assigned_count=non_assigned_count+length(filters_to_check);
                non_assigned_correlations(non_assigned_count-length(filters_to_check)+1:non_assigned_count)=corr_vec;
            else
                temp_corr_vec=corr_vec;
                temp_corr_vec(highest_corr_ind)=[];
                non_assigned_count=non_assigned_count+length(filters_to_check)-1;
                non_assigned_correlations(non_assigned_count-length(temp_corr_vec)+1:non_assigned_count)=temp_corr_vec;
            end
        else
            count=count+1;
            all_sessions_filters(initial_cell_num+count,:,:)=new_filter;
            all_sessions_centroids(initial_cell_num+count,:)=new_centroids(k,:);
            cell_to_index_map(initial_cell_num+count,:)=zeros(1,num_sessions);
            cell_to_index_map(initial_cell_num+count,n)=k;
            correlation_map(initial_cell_num+count,:)=zeros(1,num_sessions);
        end
    end
    close(h2);
end
close(h)

all_neighbor_correlations(neighbor_count+1:end)=[];
all_neighbor_distances(neighbor_count+1:end)=[];
all_assigned_correlations(assigned_count+1:end)=[];
non_assigned_correlations(non_assigned_count+1:end)=[];
non_assigned_correlations(non_assigned_correlations<0.01)=[];
if ~isempty(varargin)    
    is_in_ROI(neighbor_count+1:end)=[];
end

end

