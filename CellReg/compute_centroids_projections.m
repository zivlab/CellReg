function [centroid_projections]=compute_centroids_projections(centroid_locations,spatial_footprints)
% This function projects the centroid locations of all the cells onto the FOV

% Inputs:
% 1. centroid_locations
% 2. spatial_footprints

% Outputs:
% 1. centroid_projections

number_of_sessions=size(centroid_locations,2);

centroid_projections=cell(1,number_of_sessions);
for n=1:number_of_sessions
    this_session_centroids=centroid_locations{n};
    number_of_cells=size(this_session_centroids,1);
    this_session_spatial_footprints=spatial_footprints{n};
    normalized_centroids=zeros(size(this_session_spatial_footprints));
    for k=1:number_of_cells
        if round(this_session_centroids(k,2))>1.5 && round(this_session_centroids(k,1))>1.5 && round(this_session_centroids(k,2))<size(normalized_centroids,2)-1 && round(this_session_centroids(k,1))<size(normalized_centroids,3)-1
            normalized_centroids(k,round(this_session_centroids(k,2))-1:round(this_session_centroids(k,2))+1,round(this_session_centroids(k,1))-1:round(this_session_centroids(k,1))+1)=1/4;
            normalized_centroids(k,round(this_session_centroids(k,2))-1:round(this_session_centroids(k,2))+1,round(this_session_centroids(k,1)))=1/2;
            normalized_centroids(k,round(this_session_centroids(k,2)),round(this_session_centroids(k,1))-1:round(this_session_centroids(k,1))+1)=1/2;
            normalized_centroids(k,round(this_session_centroids(k,2)),round(this_session_centroids(k,1)))=1;
        elseif round(this_session_centroids(k,2))>0 && round(this_session_centroids(k,1))>0 && round(this_session_centroids(k,2))<size(normalized_centroids,2) && round(this_session_centroids(k,1))<size(normalized_centroids,3)
            normalized_centroids(k,round(this_session_centroids(k,2)),round(this_session_centroids(k,1)))=1;
        end
    end
    centroid_projections{n}=squeeze(sum(normalized_centroids,1));
end

end

