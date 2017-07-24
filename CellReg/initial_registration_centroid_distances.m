function [cell_to_index_map,registered_cells_centroid_distances,non_registered_cells_centroid_distances]=initial_registration_centroid_distances(maximal_distance,centroid_distance_threshold,centroid_locations)
% This function performs an initial cell registration across sessios
% based on a chosen centroid distance thresohld. 

% Inputs:
% 1. maximal_distance
% 2. centroid_distance_threshold
% 3. centroid_locations

% Outputs:
% 1. cell_to_index_map
% 2. registered_cells_centroid_distances
% 3. non_registered_cells_centroid_distances

number_of_sessions=size(centroid_locations,2);

% initializing the registration with the cells from session #1:
initial_number_of_cells=size(centroid_locations{1},1);
cell_to_index_map=zeros(initial_number_of_cells,number_of_sessions);
cell_to_index_map(:,1)=1:initial_number_of_cells;
centroid_distance_map=zeros(initial_number_of_cells,number_of_sessions);
registered_centroid_locations=centroid_locations{1};

% allocating space:
count=0;
duplicate_match_count=0;
neighbor_centroid_distances=zeros(1,number_of_sessions^2*5000);
registered_cells_centroid_distances=zeros(1,number_of_sessions^2*5000);
non_registered_cells_centroid_distances=zeros(1,number_of_sessions^2*5000);

neighbor_count=0;
assigned_count=0;
non_assigned_count=0;
num_candidates=0;
% Registering cells according  in the order of the sessions:
for n=2:number_of_sessions;
    new_centroids=centroid_locations{n};
    for k=1:size(new_centroids,1) % for each cell
        is_assigned=0;
        centroid=repmat(new_centroids(k,:),size(registered_centroid_locations,1),1);
        all_distances=sqrt(sum((centroid-registered_centroid_locations).^2,2));
        centroid_locations_to_check=find(all_distances<maximal_distance);
        if ~isempty(centroid_locations_to_check) % finding the best candidate for each cell
            num_candidates=num_candidates+length(centroid_locations_to_check);
            distance_vec=zeros(1,length(centroid_locations_to_check));
            for m=1:length(centroid_locations_to_check)
                neighbor_count=neighbor_count+1;
                distance_vec(m)=all_distances(centroid_locations_to_check(m));
                neighbor_centroid_distances(neighbor_count)=all_distances(centroid_locations_to_check(m));
            end
            [lowest_dist,lowest_dist_ind]=min(distance_vec);
            if lowest_dist>centroid_distance_threshold % no registration - new cell to list
                count=count+1;
                registered_centroid_locations(initial_number_of_cells+count,:)=new_centroids(k,:);
                cell_to_index_map(initial_number_of_cells+count,:)=zeros(1,number_of_sessions);
                cell_to_index_map(initial_number_of_cells+count,n)=k;
                centroid_distance_map(initial_number_of_cells+count,:)=zeros(1,number_of_sessions);
                centroid_distance_map(initial_number_of_cells+count,n)=lowest_dist;
            else % check if there is already a registered cell
                index=centroid_locations_to_check(lowest_dist_ind);
                if cell_to_index_map(index,n)==0 % register cells together
                    cell_to_index_map(index,n)=k;
                    centroid_distance_map(index,n)=lowest_dist;
                    assigned_count=assigned_count+1;
                    is_assigned=1;
                    registered_cells_centroid_distances(1,assigned_count)=lowest_dist;
                else % there is already a registered cell
                    duplicate_match_count=duplicate_match_count+1;
                    if lowest_dist<centroid_distance_map(index,n) % switch between cells
                        count=count+1;
                        switch_cell=cell_to_index_map(index,n);
                        switch_centroid=(centroid_locations{n}(switch_cell,:));
                        registered_centroid_locations(initial_number_of_cells+count,:)=switch_centroid;
                        cell_to_index_map(initial_number_of_cells+count,n)=switch_cell;
                        centroid_distance_map(initial_number_of_cells+count,:)=zeros(1,number_of_sessions);
                        cell_to_index_map(index,n)=k;
                        centroid_distance_map(index,n)=lowest_dist;
                        assigned_count=assigned_count+1;
                        is_assigned=1;
                        registered_cells_centroid_distances(1,assigned_count)=lowest_dist;
                    else % no registration - new cell to list
                        count=count+1;
                        registered_centroid_locations(initial_number_of_cells+count,:)=new_centroids(k,:);
                        cell_to_index_map(initial_number_of_cells+count,:)=zeros(1,number_of_sessions);
                        centroid_distance_map(initial_number_of_cells+count,:)=zeros(1,number_of_sessions);
                        cell_to_index_map(initial_number_of_cells+count,n)=k;
                        centroid_distance_map(initial_number_of_cells+count,:)=zeros(1,number_of_sessions);
                        centroid_distance_map(initial_number_of_cells+count,n)=lowest_dist;
                    end
                end
            end
            if is_assigned==0 % add this value to the non registered cells
                non_assigned_count=non_assigned_count+length(centroid_locations_to_check);
                non_registered_cells_centroid_distances(non_assigned_count-length(centroid_locations_to_check)+1:non_assigned_count)=distance_vec;
            else % add this value to the registered cells
                temp_distance_vec=distance_vec;
                temp_distance_vec(lowest_dist_ind)=[];
                non_assigned_count=non_assigned_count+length(centroid_locations_to_check)-1;
                non_registered_cells_centroid_distances(non_assigned_count-length(temp_distance_vec)+1:non_assigned_count)=temp_distance_vec;
            end
        else % no candidates - new cell to list
            count=count+1;
            registered_centroid_locations(initial_number_of_cells+count,:)=new_centroids(k,:);
            cell_to_index_map(initial_number_of_cells+count,:)=zeros(1,number_of_sessions);
            cell_to_index_map(initial_number_of_cells+count,n)=k;
            centroid_distance_map(initial_number_of_cells+count,:)=zeros(1,number_of_sessions);
        end
    end
end

registered_cells_centroid_distances(assigned_count+1:end)=[];
non_registered_cells_centroid_distances(non_assigned_count+1:end)=[];

end

