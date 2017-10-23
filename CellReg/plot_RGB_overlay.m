function plot_RGB_overlay(spatial_footprints_projections,RGB_indexes,overlapping_FOV)
% This function plot an RGB overaly of the spatial spatial_footprints projections
% of three specified sessions.

% Inputs:
% 1. spatial_footprints_projections
% 2. RGB_indexes - the three specified sessions
% 3. overlapping_FOV

number_of_sessions=size(spatial_footprints_projections,2);

if max(RGB_indexes)<=number_of_sessions
    y_size=size(spatial_footprints_projections{1},1);
    x_size=size(spatial_footprints_projections{1},2);
    spatial_footprints_projections_rgb=zeros(y_size,x_size,3);
    spatial_footprints_projections_rgb(:,:,1)=spatial_footprints_projections{RGB_indexes(1)};
    spatial_footprints_projections_rgb(:,:,2)=spatial_footprints_projections{RGB_indexes(2)};
    if number_of_sessions>2
        spatial_footprints_projections_rgb(:,:,3)=spatial_footprints_projections{RGB_indexes(3)};
    end
    spatial_footprints_projections_rgb(spatial_footprints_projections_rgb>1)=1;
    spatial_footprints_projections_rgb=spatial_footprints_projections_rgb+0.25*max(max(max(spatial_footprints_projections_rgb)))*repmat(overlapping_FOV,1,1,3);
    spatial_footprints_projections_rgb(spatial_footprints_projections_rgb>1)=1;
    
    imagesc(spatial_footprints_projections_rgb)
    if number_of_sessions>2
        title_string=['RGB overlay - sessions ' num2str(RGB_indexes(1)) ', ' num2str(RGB_indexes(2)) ', and ' num2str(RGB_indexes(3))];
    else
        title_string=['RGB overlay - sessions ' num2str(RGB_indexes(1)) ' and ' num2str(RGB_indexes(2))];
    end
    title(title_string,'FontWeight','Bold','fontsize',10)
    set(gca,'xtick',[])
    set(gca,'ytick',[])
else
    error('The sepcified sessions are not valid')
end

end

