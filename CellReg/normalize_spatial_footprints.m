function [normalized_spatial_footprints]=normalize_spatial_footprints(spatial_footprints)
% This function normalizes the spatial footprints to sum up to 1.
% A threshold for pixel weight can also be used by setting lower pixels to
% zero.

% Inputs:
% 1. spatial_footprints

% Outputs:
% 1. normalized_spatial_footprints

pixel_threshold=0; % to change pixels with low signal to zero (between 0-1)

number_of_sessions=size(spatial_footprints,2);
normalized_spatial_footprints=cell(1,number_of_sessions);
disp('Normalizing spatial footprints:')
display_progress_bar('Terminating previous progress bars',true)    
for n=1:number_of_sessions
    display_progress_bar(['Normalizing spatial footprints for session #' num2str(n) ' - '],false)    
    this_session_spatial_footprints=spatial_footprints{n};
    normalized_spatial_footprints{n}=zeros(size(this_session_spatial_footprints));
    number_of_cells=size(this_session_spatial_footprints,1);
    for k=1:number_of_cells
        display_progress_bar(100*(k)/(number_of_cells),false)
        temp_spatial_footprint=double(squeeze(this_session_spatial_footprints(k,:,:)));
        max_spatial_footprint=max(temp_spatial_footprint(:));
        temp_spatial_footprint(temp_spatial_footprint<pixel_threshold*max_spatial_footprint)=0;        
        temp_spatial_footprint=temp_spatial_footprint./sum(sum(temp_spatial_footprint));
        normalized_spatial_footprints{n}(k,:,:)=temp_spatial_footprint;
    end
    display_progress_bar(' done',false);
end

end

