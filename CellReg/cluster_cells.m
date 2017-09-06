function [cell_to_index_map,clusters_centroid_locations,varargout]=cluster_cells(cell_to_index_map,all_to_all_p_same,all_to_all_indexes,maximal_distance,registration_threshold,centroid_locations,registration_approach,transform_data)
% This function registers cells with a clustering procedure.
% it uses either the P_same that were calculated for each cell-pair, or
% alternatively their centroid distance or spatial similarity.
% Each cell from a given session is registered to the cluster for which it
% is most similar using one of the clustering criterions mentioned below.

% Inputs:
% 1. cell_to_index_map - list of initial registered cells and their index in each session
% 2. all_to_all_p_same
% 3. all_to_all_indexes
% 4. maximal_distance
% 5. registration_threshold
% 6. centroid_locations
% 7. registration_approach
% 8. transform_data - 'true' to transformed distances to similarities

% Outputs:
% 1. cell_to_index_map - list of registered cells and indexes after clustering
% 2. registered_cells_centroids
% 3. varargout
%   3{1}. cell_scores
%   3{2}. cell_scores_positive
%   3{3}. cell_scores_negative
%   3{4}. cell_scores_exclusive

maximal_number_of_iterations=10; % clustering usually convereges after 1 iteration
cluster_distance_threshold=1.7*maximal_distance; % distance between clusters may be slightly larger than distance between cell-pairs
number_of_sessions=size(centroid_locations,2);

% Different options for clustering criterion:
decision_type='Maximal similarity';
% decision_type='Minimal dissimilarity';
% decision_type='Average similarity';
num_changes_thresh=10;
iteration=1;

% the next variables indicate if clustering converges:
changes_count=zeros(1,maximal_number_of_iterations);
switch_count=zeros(1,maximal_number_of_iterations);
separation_count=zeros(1,maximal_number_of_iterations);
move_count=zeros(1,maximal_number_of_iterations);
delete_count=zeros(1,maximal_number_of_iterations);
changes_count(1)=-1;

% Finding the optimal clustering with an iterative process:
disp('Clustering cells:');
display_progress_bar('Terminating previous progress bars',true)    
while (changes_count(iteration)>num_changes_thresh || changes_count(iteration)==-1) && (iteration)<maximal_number_of_iterations    
    display_progress_bar(['Performing clustering (iteration #' num2str(iteration) ') - '],false)
    
    % Finding the center of mass of all clusters:
    iteration=iteration+1;
    num_clusters=size(cell_to_index_map,1);
    clusters_centroid_locations=zeros(2,num_clusters);
    for n=1:num_clusters
        cells_in_cluster=cell_to_index_map(n,:);
        sessions_ind=find(cells_in_cluster>0);
        N=length(sessions_ind);
        centroids_found=zeros(2,N);
        for k=1:N
            centroids_found(:,k)=centroid_locations{sessions_ind(k)}(cells_in_cluster(sessions_ind(k)),:);
        end
        clusters_centroid_locations(:,n)=mean(centroids_found,2);
    end
    
    % Checking for each cell which cluster is the most correlated:
    for n=1:number_of_sessions % for each session
        num_cells=sum(cell_to_index_map(:,n)>0);
        cluster_ind=find(cell_to_index_map(:,n)>0);
        for k=1:num_cells % for each cell
            this_cell=cell_to_index_map(cluster_ind(k),n);
            this_cell_centroid=centroid_locations{n}(this_cell,:);
            centroid=repmat(this_cell_centroid,num_clusters,1)';
            distance_vec=sqrt(sum((centroid-clusters_centroid_locations).^2));
            clusters_to_check=find(distance_vec<cluster_distance_threshold);
            num_candidates=length(clusters_to_check);
            total_similarity=zeros(1,num_candidates);
            normalization_factor=zeros(1,num_candidates);
            max_in_cluster_similarity=zeros(1,num_candidates);
            min_in_cluster_similarity=ones(1,num_candidates);
            for m=1:num_candidates % going over all candidates
                this_cluster=clusters_to_check(m);
                cells_in_cluster=cell_to_index_map(this_cluster,:);
                % checking if cell from this session is already in the cluster
                sessions_in_cluster=find(cells_in_cluster>0);
                if sum(sessions_in_cluster-n==0)>0 % cell from same session already in cluster
                    sessions_in_cluster(sessions_in_cluster-n==0)=[];
                end
                if ~isempty(sessions_in_cluster)
                    for l=1:length(sessions_in_cluster)
                        in_cluster=cell_to_index_map(this_cluster,sessions_in_cluster(l));
                        in_session=all_to_all_indexes{n}{this_cell,sessions_in_cluster(l)};
                        in_cluster_ind=find(in_session==in_cluster);
                        if ~isempty(all_to_all_p_same{n}{this_cell,sessions_in_cluster(l)}(in_cluster_ind))
                            if transform_data
                                temp_value=all_to_all_p_same{n}{this_cell,sessions_in_cluster(l)}(in_cluster_ind);
                                this_similarity=transform_distance_to_similarity(temp_value,maximal_distance);
                            else
                                this_similarity=all_to_all_p_same{n}{this_cell,sessions_in_cluster(l)}(in_cluster_ind);
                            end
                            total_similarity(m)=total_similarity(m)+this_similarity;
                            normalization_factor(m)=normalization_factor(m)+1;
                            max_in_cluster_similarity(m)=max(max_in_cluster_similarity(m),this_similarity);
                            min_in_cluster_similarity(m)=min(min_in_cluster_similarity(m),this_similarity);
                        else
                            min_in_cluster_similarity(m)=0;
                        end
                    end
                else
                    min_in_cluster_similarity(m)=0;
                end
            end
            
            % different criterions that can be used for clustering
            if strcmp(decision_type,'Minimal dissimilarity')
                [max_similarity,max_similarity_ind]=max(min_in_cluster_similarity);
            elseif strcmp(decision_type,'Average similarity')
                normalized_similarity=total_similarity./(normalization_factor);
                [max_similarity,max_similarity_ind]=max(normalized_similarity);
                average_similarity=max_similarity;
            elseif strcmp(decision_type,'Maximal similarity')
                [max_similarity,max_similarity_ind]=max(max_in_cluster_similarity);
            end
            
            max_cluster=clusters_to_check(max_similarity_ind);
            cells_in_max_cluster=cell_to_index_map(max_cluster,:);
            sessions_in_max_cluster=find(cells_in_max_cluster>0);
            if sum(sessions_in_max_cluster-n==0)>0 % cell from same session already in cluster
                sessions_in_max_cluster(sessions_in_max_cluster-n==0)=[];
            end
            
            if ~isempty(sessions_in_max_cluster)
                cells_in_original_cluster=cell_to_index_map(cluster_ind(k),:);
                num_cells_in_original_cluster=sum(cells_in_original_cluster>0);
                
                if strcmp(decision_type,'Minimal dissimilarity')
                    average_similarity=max_similarity;
                elseif strcmp(decision_type,'Maximal similarity')
                    average_similarity=max_similarity;
                end
                
                if average_similarity<registration_threshold && num_cells_in_original_cluster>1 % split the cell to a new cluster
                    num_clusters=num_clusters+1;
                    cell_to_index_map(num_clusters,:)=zeros(1,number_of_sessions);
                    cell_to_index_map(num_clusters,n)=this_cell;
                    cell_to_index_map(cluster_ind(k),n)=0;
                    clusters_centroid_locations(:,num_clusters)=this_cell_centroid;
                    changes_count(iteration)=changes_count(iteration)+1;
                    separation_count(iteration)=separation_count(iteration)+1;
                elseif average_similarity>=registration_threshold
                    if (clusters_to_check(max_similarity_ind)~=cluster_ind(k)) % need to change clusters
                        this_cluster=clusters_to_check(max_similarity_ind);
                        cells_in_cluster=cell_to_index_map(this_cluster,:);
                        % checking if cell from this session is already in the cluster
                        sessions_in_cluster=find(cells_in_cluster>0);
                        if sum(sessions_in_cluster-n==0)>0 % cell from same session already in cluster - switch
                            temp_cell=cell_to_index_map(this_cluster,n);
                            temp_cell_centroid=centroid_locations{n}(temp_cell,:);
                            temp_similarity=0;
                            % switch only if this cell is more correlated to the cluster
                            cells_in_cluster=cell_to_index_map(this_cluster,:);
                            % checking if cell from this session is already in the cluster
                            sessions_in_cluster=find(cells_in_cluster>0);
                            if sum(sessions_in_cluster-n==0)>0 % cell from same session already in cluster
                                sessions_in_cluster(sessions_in_cluster-n==0)=[];
                            end
                            for l=1:length(sessions_in_cluster)
                                in_cluster=cell_to_index_map(this_cluster,sessions_in_cluster(l));
                                in_session=all_to_all_indexes{n}{temp_cell,sessions_in_cluster(l)};
                                in_cluster_ind=find(in_session==in_cluster);
                                if ~isempty(all_to_all_p_same{n}{temp_cell,sessions_in_cluster(l)}(in_cluster_ind))
                                    if transform_data
                                        temp_value=all_to_all_p_same{n}{temp_cell,sessions_in_cluster(l)}(in_cluster_ind);
                                        this_similarity=transform_distance_to_similarity(temp_value,maximal_distance);
                                        temp_similarity=temp_similarity+this_similarity;
                                    else
                                        temp_similarity=temp_similarity+all_to_all_p_same{n}{temp_cell,sessions_in_cluster(l)}(in_cluster_ind);
                                    end
                                end
                            end
                            if max_similarity>temp_similarity
                                num_clusters=num_clusters+1;
                                cell_to_index_map(this_cluster,n)=this_cell;
                                cell_to_index_map(num_clusters,:)=zeros(1,number_of_sessions);
                                cell_to_index_map(num_clusters,n)=temp_cell;
                                changes_count(iteration)=changes_count(iteration)+1;
                                switch_count(iteration)=switch_count(iteration)+1;
                                cell_to_index_map(cluster_ind(k),n)=0;
                                clusters_centroid_locations(:,num_clusters)=temp_cell_centroid;
                            end
                        else % just need to add the cell to the cluster
                            cell_to_index_map(cluster_ind(k),n)=0;
                            cell_to_index_map(this_cluster,n)=this_cell;
                            changes_count(iteration)=changes_count(iteration)+1;
                            move_count(iteration)=move_count(iteration)+1;
                        end
                    end
                end
            end
        end
    end
    
    % deleting clusters with no cells left
    cell_to_index_map_temp=cell_to_index_map;
    num_clusters=size(cell_to_index_map_temp,1);
    for n=1:num_clusters
        if sum(cell_to_index_map_temp(n,:))==0
            cell_to_index_map(n-delete_count(iteration),:)=[];
            delete_count(iteration)=delete_count(iteration)+1;
        end
    end
    cell_to_index_map_temp=cell_to_index_map;
    
    % Merging clusters if they cross the threshold:
    num_clusters=size(cell_to_index_map_temp,1);
    clusters_centroid_locations=zeros(2,num_clusters);
    for n=1:num_clusters % finding the centroid locations of the clusters
        cells_in_cluster=cell_to_index_map_temp(n,:);
        sessions_ind=find(cells_in_cluster>0);
        N=length(sessions_ind);
        centroids_found=zeros(2,N);
        for k=1:N
            centroids_found(:,k)=centroid_locations{sessions_ind(k)}(cells_in_cluster(sessions_ind(k)),:);
        end
        clusters_centroid_locations(:,n)=mean(centroids_found,2);
    end    
    for n=1:num_clusters % going over each cluster to see if it should merge
        this_cluster_cells=cell_to_index_map_temp(n,:);
        this_cluster_sessions=find(this_cluster_cells>0);
        centroid=repmat(clusters_centroid_locations(:,n),1,num_clusters);
        distance_vec=sqrt(sum((centroid-clusters_centroid_locations).^2));
        clusters_to_check=find(distance_vec<cluster_distance_threshold);
        clusters_to_check=setdiff(clusters_to_check,n);
        num_candidates=length(clusters_to_check);
        for k=1:num_candidates % candidate clusters to merge with
            candidate_cells=cell_to_index_map_temp(clusters_to_check(k),:);
            candidate_sessions=find(candidate_cells>0);
            if isempty(intersect(find(this_cluster_cells>0),find(candidate_cells>0)))
                all_to_all_temp=zeros(sum(this_cluster_cells>0),sum(candidate_cells>0));
                for m=1:sum(this_cluster_cells>0);
                    for l=1:sum(candidate_cells>0)
                        neigbor_cell=all_to_all_indexes{this_cluster_sessions(m)}{this_cluster_cells(this_cluster_sessions(m)),candidate_sessions(l)};
                        this_cell_ind=find(neigbor_cell==candidate_cells(candidate_sessions(l)));
                        if ~isempty(all_to_all_p_same{this_cluster_sessions(m)}{this_cluster_cells(this_cluster_sessions(m)),candidate_sessions(l)}(this_cell_ind))
                            if transform_data
                                temp_value=all_to_all_p_same{this_cluster_sessions(m)}{this_cluster_cells(this_cluster_sessions(m)),candidate_sessions(l)}(this_cell_ind);
                                this_similarity=transform_distance_to_similarity(temp_value,maximal_distance);
                                all_to_all_temp(m,l)=this_similarity;
                            else
                                all_to_all_temp(m,l)=all_to_all_p_same{this_cluster_sessions(m)}{this_cluster_cells(this_cluster_sessions(m)),candidate_sessions(l)}(this_cell_ind);
                            end
                        end
                    end
                end
                % the different criterions to merge clusters
                if strcmp(decision_type,'Average similarity')
                    average_all_to_all_temp_1=mean(all_to_all_temp);
                    average_all_to_all_temp_2=mean(all_to_all_temp,2);
                    over_thresh_1=sum(average_all_to_all_temp_1>registration_threshold);
                    over_thresh_2=sum(average_all_to_all_temp_2>registration_threshold);
                    if over_thresh_2==sum(this_cluster_cells>0) && over_thresh_1==sum(candidate_cells>0)
                        cell_to_index_map_temp(n,candidate_sessions)=cell_to_index_map_temp(clusters_to_check(k),candidate_sessions);
                        cell_to_index_map_temp(clusters_to_check(k),:)=zeros(1,number_of_sessions);
                        changes_count(iteration)=changes_count(iteration)+1;
                        clusters_centroid_locations(:,clusters_to_check(k))=[1000 1000];
                    end
                elseif strcmp(decision_type,'Maximal similarity')
                    max_all_to_all_temp=max(all_to_all_temp(:));
                    if max_all_to_all_temp>registration_threshold;
                        cell_to_index_map_temp(n,candidate_sessions)=cell_to_index_map_temp(clusters_to_check(k),candidate_sessions);
                        cell_to_index_map_temp(clusters_to_check(k),:)=zeros(1,number_of_sessions);
                        changes_count(iteration)=changes_count(iteration)+1;
                        clusters_centroid_locations(:,clusters_to_check(k))=[1000 1000];
                    end
                elseif strcmp(decision_type,'Minimal dissimilarity')
                    average_all_to_all_temp_1=min(all_to_all_temp);
                    average_all_to_all_temp_2=min(all_to_all_temp,2);
                    over_thresh_1=sum(average_all_to_all_temp_1>registration_threshold);
                    over_thresh_2=sum(average_all_to_all_temp_2>registration_threshold);
                    if over_thresh_2==sum(this_cluster_cells>0) && over_thresh_1==sum(candidate_cells>0)
                        cell_to_index_map_temp(n,candidate_sessions)=cell_to_index_map_temp(clusters_to_check(k),candidate_sessions);
                        cell_to_index_map_temp(clusters_to_check(k),:)=zeros(1,number_of_sessions);
                        changes_count(iteration)=changes_count(iteration)+1;
                        clusters_centroid_locations(:,clusters_to_check(k))=[1000 1000];
                    end
                end
            end
        end
    end
    
    cell_to_index_map=cell_to_index_map_temp;
    num_clusters=size(cell_to_index_map_temp,1);
    temp_delete_count=0;
    for n=1:num_clusters
        if sum(cell_to_index_map_temp(n,:))==0
            cell_to_index_map(n-temp_delete_count,:)=[];
            temp_delete_count=temp_delete_count+1;
        end
    end
    delete_count(iteration)=delete_count(iteration)+temp_delete_count;
    display_progress_bar(' done',false)
end

% Calculating the final cell scores:
if strcmp(registration_approach,'Probabilistic')
    [cell_scores,cell_scores_positive,cell_scores_negative,cell_scores_exclusive,p_same_registered_pairs]=compute_scores(cell_to_index_map,all_to_all_indexes,all_to_all_p_same,number_of_sessions);
    varargout{1}=cell_scores;
    varargout{2}=cell_scores_positive;
    varargout{3}=cell_scores_negative;
    varargout{4}=cell_scores_exclusive; 
    varargout{5}=p_same_registered_pairs; 
end

% Calculating neighbor vs corr/distance curves
num_clusters=size(cell_to_index_map,1);
clusters_centroid_locations=zeros(2,num_clusters);
for n=1:num_clusters
    cells_in_cluster=cell_to_index_map(n,:);
    sessions_ind=find(cells_in_cluster>0);
    N=length(sessions_ind);
    centroids_found=zeros(2,N);
    for k=1:N
        centroids_found(:,k)=centroid_locations{sessions_ind(k)}(cells_in_cluster(sessions_ind(k)),:);
    end
    clusters_centroid_locations(:,n)=mean(centroids_found,2);
end

end
