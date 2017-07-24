function [similarity]=transform_distance_to_similarity(measured_distance,maximal_distance)
% This function recieves a centroid distance and transforms it into a 0-1
% similarity measure:

% Inputs:
% 1. measured_distance
% 2. maximal_distance
% 3. microns_per_pixel

% Outputs:
% 1. Similarity

similarity=(maximal_distance-measured_distance)/maximal_distance;


end

