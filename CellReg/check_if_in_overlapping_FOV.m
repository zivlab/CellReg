function [is_in_overlapping_FOV] =check_if_in_overlapping_FOV(registered_cells_centroids,overlapping_FOV)
% This function checks which of the registered cells are in the overlapping
% FOV

% Inputs:
% 1. registered_cells_centroids
% 2. overlapping_FOV

% Outputs:
% 1. is_in_overlapping_FOV

number_of_cells=size(registered_cells_centroids,2);
is_in_overlapping_FOV=false(1,number_of_cells);
for n=1:number_of_cells
    if round(registered_cells_centroids(2,n))>0 && round(registered_cells_centroids(2,n))<=size(overlapping_FOV,1) && round(registered_cells_centroids(1,n))>0 && round(registered_cells_centroids(1,n))<=size(overlapping_FOV,2);
        if overlapping_FOV(round(registered_cells_centroids(2,n)),round(registered_cells_centroids(1,n)))>0
            is_in_overlapping_FOV(n)=true;
        end
    end
end

end

