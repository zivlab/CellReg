function [spatial_footprints,number_of_sessions]=load_multiple_sessions(file_names)
% This function loads the files with the spatial footprints and centroid
% locations for all the sessions.

% Inputs: 
% 1. file_names - a cell array

% Outputs: 
% 1. spatial_footprints - cell array with spatial footprints from all sessions
% 2. number_of_sessions

number_of_sessions=size(file_names,2);
if number_of_sessions<2
    error('To add a single session please use the add session button  ')
else        
    spatial_footprints=cell(1,number_of_sessions);
    for n=1:number_of_sessions
        this_file_name=file_names{1,n};
        temp_data=load(this_file_name);
        if isstruct(temp_data)
            field_name=fieldnames(temp_data);
            spatial_footprints{n}=getfield(temp_data,field_name{1});
        else
            spatial_footprints{n}=temp_data;
        end
    end
end

end

