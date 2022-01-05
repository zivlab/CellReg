function [adjusted_spatial_footprints,adjusted_FOV_all_sessions,...
    adjusted_x_size,adjusted_y_size,adjustment_zero_padding]=adjust_FOV_size(spatial_footprints)
% This function adjusts the sizes of the FOV of each session to be the same

% Inputs:
% 1. spatial_footprints

% Outputs:
% 1. adjusted_spatial_footprints
% 2. adjusted_FOV_all_sessions - pixels that are within the FOV are equal to 1
% 3. adjusted_x_size - in pixels (the same for all sessions)
% 4. adjusted_y_size
% 5. adjustment_zero_padding - (x,y) adjustment_zero_padding for each session that were used for adjustment

number_of_sessions=size(spatial_footprints,2);

% Finding a FOV size the fits all sessions:
adjusted_x_size=0;
adjusted_y_size=0;
FOV_all_sessions=cell(1,number_of_sessions);

for n=1:number_of_sessions
    footprint_dat = get_spatial_footprints(spatial_footprints{n});
    FOV_all_sessions{n}=ones(footprint_dat.size(2), footprint_dat.size(3));
    if n==1
        adjusted_x_size = footprint_dat.size(3);
        adjusted_y_size = footprint_dat.size(2);
    else
        adjusted_x_size=max(adjusted_x_size,footprint_dat.size(3));
        adjusted_y_size=max(adjusted_y_size,footprint_dat.size(2));
    end
end

% Correcting the footprints and centroids according to the adjustments:
adjusted_spatial_footprints=cell(1,number_of_sessions);
adjusted_FOV_all_sessions=zeros(adjusted_y_size,adjusted_x_size,number_of_sessions);
adjustment_zero_padding=zeros(2,number_of_sessions);
for n=1:number_of_sessions   
    adjusted_FOV_temp=FOV_all_sessions{n};    
    new_adjusted_FOV=zeros(adjusted_y_size,adjusted_x_size);
    new_adjusted_FOV(1:size(adjusted_FOV_temp,1),1:size(adjusted_FOV_temp,2))=adjusted_FOV_temp;
    adjusted_FOV_all_sessions(:,:,n)=new_adjusted_FOV;
    
    footprints_info = get_spatial_footprints(spatial_footprints{n});
    spatial_footprints_temp = footprints_info.load_footprints;
    spatial_footprints_temp = spatial_footprints_temp.footprints;
    
    [number_of_cells, sz_y, sz_x] = size(spatial_footprints_temp);
    adjusted_spatial_footprints_temp = zeros(number_of_cells,adjusted_y_size,adjusted_x_size);
    adjusted_spatial_footprints_temp(:,1:sz_y,1:sz_x) = spatial_footprints_temp;
    
    if ~isempty(footprints_info.write2path)
        footprint = mat_to_sparse_cell(adjusted_spatial_footprints_temp);
        adjusted_spatial_footprints{n} = [footprints_info.write2path, filesep,...
            'adjusted_spatial_footprints_',num2str(n),'.mat'];
        save([footprints_info.write2path, filesep, 'adjusted_spatial_footprints_',...
            num2str(n),'.mat'],'footprint','-v7.3')
    else
        adjusted_spatial_footprints{n}=adjusted_spatial_footprints_temp;
    end
    adjustment_zero_padding(1,n)=adjusted_x_size-sz_x;
    adjustment_zero_padding(2,n)=adjusted_y_size-sz_y;
    
    % all(abs(adjusted_spatial_footprints{1}(:)-adjusted_spatial_footprints2{1}(:)) < eps)
end

end