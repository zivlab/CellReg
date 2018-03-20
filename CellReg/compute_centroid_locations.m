function [centroid_locations]=compute_centroid_locations(spatial_footprints,microns_per_pixel)
% This function computes the centroid locations for all the cells according
% to their spatial footprints (spatial_footprints).

% Inputs:
% 1. spatial_footprints
% 2. microns_per_pixel

% Outputs:
% 1. centroid_locations

typical_cell_size=10; % in micrometers - determines the radius that is used for gaussfit
normalized_typical_cell_size=typical_cell_size/microns_per_pixel;
gaussian_radius=round(2*normalized_typical_cell_size); % for estimation of centroid

number_of_sessions=size(spatial_footprints,2);
centroid_locations=cell(1,number_of_sessions);
disp('Calculating centroid locations:');
display_progress_bar('Terminating previous progress bars',true)    
for n=1:number_of_sessions
    if number_of_sessions>1
        display_progress_bar(['Calculating centroid locations for session #' num2str(n) ' - '],false)
    else
        display_progress_bar('Calculating centroid locations for this session - ',false)
    end
    this_session_spatial_footprints=spatial_footprints{n};
    num_spatial_footprints=size(this_session_spatial_footprints,1);
    centroid_locations{n}=zeros(num_spatial_footprints,2);
    for k=1:num_spatial_footprints
        display_progress_bar(100*(k)/(num_spatial_footprints),false)
        temp_spatial_footprint=squeeze(this_session_spatial_footprints(k,:,:));

        % calculating x and y projections of the spatial footprints for the
        % calculation of the center of mass:
        x_projection=sum(temp_spatial_footprint);
        y_projection=sum(temp_spatial_footprint,2)';
        [~,max_x_ind]=max(x_projection);
        [~,max_y_ind]=max(y_projection);        
        
        % zero padding of x and y projections:
        if max_x_ind>gaussian_radius && max_x_ind<=length(x_projection)-gaussian_radius
        localized_x_projection=x_projection(max_x_ind-gaussian_radius:max_x_ind+gaussian_radius);
        elseif max_x_ind<=gaussian_radius
            zero_padding_size=gaussian_radius-max_x_ind+1;
            localized_x_projection=[zeros(1,zero_padding_size) , x_projection(1:max_x_ind+gaussian_radius)];
        elseif max_x_ind>length(x_projection)-gaussian_radius
            zero_padding_size=gaussian_radius-(length(x_projection)-max_x_ind);
            localized_x_projection=[x_projection(max_x_ind-gaussian_radius:end), zeros(1,zero_padding_size)];
        end
        
        if max_y_ind>gaussian_radius && max_y_ind<=length(y_projection)-gaussian_radius
        localized_y_projection=y_projection(max_y_ind-gaussian_radius:max_y_ind+gaussian_radius);
        elseif max_y_ind<=gaussian_radius
            zero_padding_size=gaussian_radius-max_y_ind+1;
            localized_y_projection=[zeros(1,zero_padding_size) , y_projection(1:max_y_ind+gaussian_radius)];
        elseif max_y_ind>length(y_projection)-gaussian_radius
            zero_padding_size=gaussian_radius-(length(y_projection)-max_y_ind);
            localized_y_projection=[y_projection(max_y_ind-gaussian_radius:end) , zeros(1,zero_padding_size)];
        end
        
        % Calculating the center of mass with a gaussian fit (only for weighted ROIs):
        [~,centroid_x_temp]=gaussfit(-gaussian_radius:gaussian_radius,localized_x_projection./sum(localized_x_projection),0.5*normalized_typical_cell_size,0);
        [~,centroid_y_temp]=gaussfit(-gaussian_radius:gaussian_radius,localized_y_projection./sum(localized_y_projection),0.5*normalized_typical_cell_size,0);    
        centroid_x=max_x_ind+centroid_x_temp;
        centroid_y=max_y_ind+centroid_y_temp;
        centroid_locations{n}(k,:)=[centroid_x,centroid_y];
    end
    display_progress_bar(' done',false);
end

end


