function [spatial_footprints]=load_single_session(file_name)
% This function loads the files with the spatial footprints and centroid
% locations for all the sessions.

% Inputs:
% 1. file_name - cell

% Outputs:
% 1. spatial_footprints - a matrix with the added spatial footprints

temp_data=load(file_name);
if isstruct(temp_data)
    field_name=fieldnames(temp_data);
    spatial_footprints=getfield(temp_data,field_name{1});
else
    spatial_footprints=temp_data;
end


end

