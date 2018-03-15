function [all_to_all_indexes,all_to_all_spatial_correlations,all_to_all_centroid_distances,neighbors_spatial_correlations,neighbors_centroid_distances,neighbors_x_displacements,neighbors_y_displacements,NN_spatial_correlations,NNN_spatial_correlations,NN_centroid_distances,NNN_centroid_distances]=compute_data_distribution(spatial_footprints,centroid_locations,maximal_distance,imaging_technique)
% This function computes the distributions of distances and correlations
% for all the neighboring cells pairs acorss sessions with a distance <12 microns

% Inputs:
% 1. spatial_footprints
% 2. centroid_locations
% 3. maximal distance - in pixels 
% 4. imaging_technique


% Outputs:
% 1. all_to_all_indexes - cell-pairs with centroid_distance<maximal distance
% 2. all_to_all_spatial_correlationss - computed for cell-pairs in all_to_all_index
% 3. all_to_all_centroid_distances - computed for cell-pairs in all_to_all_index
% 4. neighbors_spatial_correlations - a concatanted vector with all the values
% 5. neighbors_centroid_distances
% 6. neighbors_x_displacements - decomposition of distance to (x,y)
% 7. neighbors_y_displacements
% 8. NN_spatial_correlations - nearest neighbors only
% 9. NNN_spatial_correlations - non-nearest neighbors only
% 10. NN_centroid_distances
% 11. NNN_centroid_distances

number_of_sessions=size(spatial_footprints,2);
cell_memory_allocation=10000;

% defining the outputs:
all_to_all_indexes=cell(1,number_of_sessions);
all_to_all_spatial_correlations=cell(1,number_of_sessions);
all_to_all_centroid_distances=cell(1,number_of_sessions);
neighbors_spatial_correlations=zeros(1,number_of_sessions^2*cell_memory_allocation);
neighbors_centroid_distances=zeros(1,number_of_sessions^2*cell_memory_allocation);
neighbors_x_displacements=zeros(1,number_of_sessions^2*cell_memory_allocation);
neighbors_y_displacements=zeros(1,number_of_sessions^2*cell_memory_allocation);
NN_spatial_correlations=zeros(1,number_of_sessions^2*cell_memory_allocation);
NNN_spatial_correlations=zeros(1,number_of_sessions^2*cell_memory_allocation);
NN_centroid_distances=zeros(1,number_of_sessions^2*cell_memory_allocation);
NNN_centroid_distances=zeros(1,number_of_sessions^2*cell_memory_allocation);

% calculating spatial correlations and cetnroid ditances for neigboring
% cells-pairs:
neighbor_count=0;
NN_count=0;
NNN_count=0;
disp('Calculating the distributions of cell-pair similarities:')
display_progress_bar('Terminating previous progress bars',true)    
for n=1:number_of_sessions 
    display_progress_bar(['Calculating spatial correlations and centroid distances for session #' num2str(n) ' - '],false)
    new_spatial_footprints=spatial_footprints{n};
    new_centroids=centroid_locations{n};
    number_of_cells=size(new_spatial_footprints,1);
    all_to_all_spatial_correlations{n}=cell(number_of_cells,number_of_sessions);
    all_to_all_centroid_distances{n}=cell(number_of_cells,number_of_sessions);
    all_to_all_indexes{n}=cell(number_of_cells,number_of_sessions);
    sessions_to_compare=1:number_of_sessions;
    sessions_to_compare(n)=[];
    for k=1:number_of_cells % for each cell
        display_progress_bar(100*(k)/(number_of_cells),false)
        new_spatial_footprint=squeeze(new_spatial_footprints(k,:,:));
        for m=1:length(sessions_to_compare)
            this_session=sessions_to_compare(m);
            this_session_centroids=centroid_locations{this_session};
            centroid=repmat(new_centroids(k,:),size(this_session_centroids,1),1);
            this_session_spatial_footprints=spatial_footprints{this_session};
            distance_vec=sqrt(sum((centroid-this_session_centroids).^2,2));
            diff_temp=centroid-this_session_centroids;
            spatial_footprints_to_check=find(distance_vec<maximal_distance);
            distance_vec_x=diff_temp(spatial_footprints_to_check,1);
            distance_vec_y=diff_temp(spatial_footprints_to_check,2);
            this_distance_vec=distance_vec(spatial_footprints_to_check);
            if ~isempty(spatial_footprints_to_check) % all neighboring cells
                corr_vec=zeros(1,length(spatial_footprints_to_check));
                num_empty_spatial_footprints=0;
                for l=1:length(spatial_footprints_to_check)
                    suspected_spatial_footprint=squeeze(this_session_spatial_footprints(spatial_footprints_to_check(l),:,:));
                    if sum(sum(suspected_spatial_footprint))==0 || sum(sum(new_spatial_footprint))==0
                        num_empty_spatial_footprints=num_empty_spatial_footprints+1;
                    else % compute spatial correlation
                        neighbor_count=neighbor_count+1;
                        corr_vec(l)=corr2(suspected_spatial_footprint,new_spatial_footprint);
                        neighbors_spatial_correlations(neighbor_count)=corr_vec(l);
                        neighbors_centroid_distances(neighbor_count)=distance_vec(spatial_footprints_to_check(l));
                        neighbors_x_displacements(neighbor_count)=distance_vec_x(l);
                        neighbors_y_displacements(neighbor_count)=distance_vec_y(l);
                    end
                end
                if num_empty_spatial_footprints<length(spatial_footprints_to_check);
                    NN_count=NN_count+1;
                    NN_spatial_correlations(NN_count)=max(corr_vec);
                    NN_centroid_distances(NN_count)=min(distance_vec(spatial_footprints_to_check));
                    if length(spatial_footprints_to_check)>1 && num_empty_spatial_footprints<length(spatial_footprints_to_check)-1
                        NNN_count=NNN_count+length(spatial_footprints_to_check)-1;
                        [dist_vec_sorted,ind]=sort(distance_vec(spatial_footprints_to_check));
                        dist_vec_sorted(1)=[];
                        NNN_centroid_distances(2+NNN_count-length(spatial_footprints_to_check):NNN_count)=dist_vec_sorted(1:end);
                        [corr_vec_sorted]=corr_vec(ind);
                        corr_vec_sorted(1)=[];
                        NNN_spatial_correlations(2+NNN_count-length(spatial_footprints_to_check):NNN_count)=corr_vec_sorted(1:end);
                    end
                end
                all_to_all_spatial_correlations{n}{k,this_session}=corr_vec;                
                all_to_all_centroid_distances{n}{k,this_session}=this_distance_vec;
                all_to_all_indexes{n}{k,this_session}=spatial_footprints_to_check;
            end
        end
    end
    display_progress_bar(' done',false);
end

% deleting empty elements in the vectors:
neighbors_spatial_correlations(neighbor_count+1:end)=[];
neighbors_centroid_distances(neighbor_count+1:end)=[];
neighbors_x_displacements(neighbor_count+1:end)=[];
neighbors_y_displacements(neighbor_count+1:end)=[];
NN_spatial_correlations(NN_count+1:end)=[];
NNN_spatial_correlations(NNN_count+1:end)=[];
NN_centroid_distances(NN_count+1:end)=[];
NNN_centroid_distances(NNN_count+1:end)=[];

% changing the values slightly if the spatial correlation is exactly 0 or 1
% or if the centroid distance is exactly 0. Otherwise the probabilistic
% modeling will fail:
neighbors_spatial_correlations_temp=neighbors_spatial_correlations;
neighbors_spatial_correlations_temp(neighbors_spatial_correlations_temp==0)=10^-10;
neighbors_spatial_correlations_temp(neighbors_spatial_correlations_temp==1)=1-10^-10;
neighbors_centroid_distances_temp=neighbors_centroid_distances;
neighbors_centroid_distances_temp(neighbors_centroid_distances_temp==0)=10^-10;
neighbors_spatial_correlations_temp(neighbors_centroid_distances>maximal_distance)=[];
neighbors_centroid_distances_temp(neighbors_centroid_distances>maximal_distance)=[];

% the maximal distance should be chosen in a way that neighboring cells
% pairs have non-zero overlap
if sum(neighbors_spatial_correlations_temp<0)/length(neighbors_spatial_correlations_temp)>0.05
    warning('A good fit might not be attainable because some of the cells seem to be smaller than expected. This could also occur if the provided microns per pixel ratio is incorrect or the maximal distance is too large')
end
if strcmp(imaging_technique,'one_photon');
    neighbors_centroid_distances_temp(neighbors_spatial_correlations_temp<0)=[];
    neighbors_spatial_correlations_temp(neighbors_spatial_correlations_temp<0)=[];
end

neighbors_spatial_correlations=neighbors_spatial_correlations_temp;
neighbors_centroid_distances=neighbors_centroid_distances_temp;

NN_spatial_correlations_temp=NN_spatial_correlations;
NNN_spatial_correlations_temp=NNN_spatial_correlations;

NN_centroid_distances_temp=NN_centroid_distances;
NNN_centroid_distances_temp=NNN_centroid_distances;
NN_spatial_correlations_temp(NN_centroid_distances>maximal_distance)=[];
NN_centroid_distances_temp(NN_centroid_distances>maximal_distance)=[];
NNN_spatial_correlations_temp(NNN_centroid_distances>maximal_distance)=[];
NNN_centroid_distances_temp(NNN_centroid_distances>maximal_distance)=[];
if strcmp(imaging_technique,'one_photon');
    NN_centroid_distances_temp(NN_spatial_correlations_temp<0)=[];
    NN_spatial_correlations_temp(NN_spatial_correlations_temp<0)=[];
    NNN_centroid_distances_temp(NNN_spatial_correlations_temp<0)=[];
    NNN_spatial_correlations_temp(NNN_spatial_correlations_temp<0)=[];
end
NN_centroid_distances=NN_centroid_distances_temp;
NNN_centroid_distances=NNN_centroid_distances_temp;
NN_spatial_correlations=NN_spatial_correlations_temp;
NNN_spatial_correlations=NNN_spatial_correlations_temp;

if length(NNN_centroid_distances_temp)<0.1*length(NN_centroid_distances_temp);
    warning('There is insufficient number of non-nearest neighboring cells to estimate the different cells distribution. This could occur if the provided microns per pixel ratio is incorrect or if the data is sparse. You can try increasing the maximal distance')
end

%(x,y) displacements distributions:
neighbors_x_displacements_temp=neighbors_x_displacements;
neighbors_y_displacements_temp=neighbors_y_displacements;
neighbors_x_displacements_temp(abs(neighbors_y_displacements_temp)>maximal_distance)=[];
neighbors_y_displacements_temp(abs(neighbors_y_displacements_temp)>maximal_distance)=[];
neighbors_y_displacements_temp(abs(neighbors_x_displacements_temp)>maximal_distance)=[];
neighbors_x_displacements_temp(abs(neighbors_x_displacements_temp)>maximal_distance)=[];
neighbors_x_displacements=neighbors_x_displacements_temp;
neighbors_y_displacements=neighbors_y_displacements_temp;

end

