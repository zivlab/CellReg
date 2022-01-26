function [spatial_footprints_corrected,centroid_locations_corrected,footprints_projections_corrected,centroid_projections_corrected,maximal_cross_correlation,best_translations,overlapping_area,varargout]=align_images(spatial_footprints,centroid_locations,footprints_projections,centroid_projections,overlapping_area_all_sessions,microns_per_pixel,reference_session_index,alignment_type,sufficient_correlation_centroids,sufficient_correlation_footprints,use_parallel_processing,varargin)

% This function recieves the spatial footprints from different sessions and
% finds the optimal alignment between their FOV's. This is the first step
% after loading the data, and it includes finding the optimal overall
% translations and rotations across sessions.

% Inputs:
% 1. spatial_footprints
% 2. centroid_locations
% 3. footprints_projections
% 4. centroid_projections
% 5. overlapping_area_all_sessions - the overlapping FOV
% 6. microns_per_pixel
% 7. reference_session_index
% 8. alignment_type - 'Translations' or 'Translations and Rotations'
% 9. sufficient_correlation_centroids % correlation between the centroid locations
% 10. sufficient_correlation_footprints % correlation between the spatial footprints
% 11. use_parallel_processing -  'true' for parallel processing
% 12. varargin
%   12{1}. maximal rotation/transformation_smoothness -  if 'Translations and
%   Rotations'/'Non-rigid' is used

% Outputs:
% 1. spatial_footprints_corrected
% 2. centroid_locations_corrected
% 3. footprints_projections_corrected
% 4. centroid_projections_corrected
% 5. maximal_cross_correlation
% 6. best_translations
% 7. overlapping_area
% 8. varargout
% 8{1}. displacement_fields - for non-rigid transormation

rotation_step=0.5; % check rotations every xx degrees
minimal_rotation=0.3; % less than this rotation in degrees does not justify rotating the cells
typical_cell_size=10; % in micrometers - determines the radius that is used for gaussfit
normalized_typical_cell_size=typical_cell_size/microns_per_pixel;

number_of_sessions=size(spatial_footprints,2);
adjusted_x_size=size(spatial_footprints{1},3);
adjusted_y_size=size(spatial_footprints{1},2);

center_of_FOV(1)=size(spatial_footprints{reference_session_index},2)/2;
center_of_FOV(2)=size(spatial_footprints{reference_session_index},3)/2;

% defining the outputs:
spatial_footprints_corrected=spatial_footprints;
centroid_locations_corrected=centroid_locations;
footprints_projections_corrected=footprints_projections;
centroid_projections_corrected=centroid_projections;

if strcmp(alignment_type,'Translations and Rotations') % if correcting for rotations as well
    maximal_rotation=varargin{1};
    rotation_vector=zeros(1,number_of_sessions-1);
    possible_rotations=-maximal_rotation:rotation_step:maximal_rotation;
    all_rotated_projections=cell(1,number_of_sessions);
    centroid_projections_rotated=cell(1,number_of_sessions);
    centroid_locations_rotated=centroid_locations;
    all_rotated_projections{reference_session_index}=footprints_projections_corrected{reference_session_index};
    centroid_projections_rotated{reference_session_index}=centroid_projections_corrected{reference_session_index};
end

maximal_cross_correlation=zeros(1,number_of_sessions-1);
if strcmp(alignment_type,'Translations and Rotations')
    best_rotations=zeros(1,number_of_sessions);
end
best_x_translations=zeros(1,number_of_sessions);
best_y_translations=zeros(1,number_of_sessions);
registration_order=setdiff(1:number_of_sessions,reference_session_index);
overlapping_area=ones(adjusted_y_size,adjusted_x_size);
overlapping_area=overlapping_area.*overlapping_area_all_sessions(:,:,reference_session_index);

% Aligning the images and cells:
display_progress_bar('Terminating previous progress bars',true)
if strcmp(alignment_type,'Non-rigid') % Non-rigid alignment:
    transformation_smoothness=varargin{1};
    best_translations=zeros(2,number_of_sessions);
    displacement_fields=zeros(number_of_sessions,adjusted_y_size,adjusted_x_size,2);
    for n=1:number_of_sessions-1
        disp(['Performing non-rigid transformation for session #' num2str(registration_order(n)) ':'])
        reference_footprints_projections_corrected=footprints_projections_corrected{reference_session_index};
        temp_footprints_projections_corrected=footprints_projections_corrected{registration_order(n)};
        [displacement_field,temp_footprints_projections_non_rigid_corrected]=imregdemons(temp_footprints_projections_corrected,reference_footprints_projections_corrected,'AccumulatedFieldSmoothing',transformation_smoothness);       
        footprints_projections_corrected{registration_order(n)}=temp_footprints_projections_non_rigid_corrected;
        this_session_footprints_unaligned=spatial_footprints_corrected{registration_order(n)};
        this_session_footprints_aligned=zeros(size(this_session_footprints_unaligned));
        number_of_cells=size(this_session_footprints_aligned,1);
        if use_parallel_processing
            disp('Performing non-rigid transformation')
            parfor m=1:number_of_cells
                unaligned_footprint=squeeze(this_session_footprints_unaligned(m,:,:));
                aligned_footprint=imwarp(unaligned_footprint,displacement_field);
                this_session_footprints_aligned(m,:,:)=aligned_footprint;
            end
        else
            display_progress_bar('Aligning spatial footprints: ',false)
            for m=1:number_of_cells
                display_progress_bar(100*(m/number_of_cells),false)
                unaligned_footprint=squeeze(this_session_footprints_unaligned(m,:,:));
                aligned_footprint=imwarp(unaligned_footprint,displacement_field);
                this_session_footprints_aligned(m,:,:)=aligned_footprint;
            end
            display_progress_bar(' done',false)
        end
        spatial_footprints_corrected{registration_order(n)}=this_session_footprints_aligned;
        [centroid_locations_corrected(registration_order(n))]=compute_centroid_locations(spatial_footprints_corrected(registration_order(n)),microns_per_pixel);
        [centroid_projections_corrected(registration_order(n))]=compute_centroids_projections(centroid_locations_corrected(registration_order(n)),spatial_footprints_corrected(registration_order(n)));
        full_FOV_correlation=normxcorr2(footprints_projections_corrected{reference_session_index},footprints_projections_corrected{registration_order(n)});
        maximal_cross_correlation(n)=max(max(full_FOV_correlation));
        displacement_fields(registration_order(n),:,:,:)=displacement_field;
    end
    varargout=cell(1,1);
    varargout{1}=displacement_fields;
else
    display_progress_bar('Terminating previous progress bars',true)
    for n=1:number_of_sessions-1
        overlapping_area_temp=overlapping_area_all_sessions(:,:,n);
        disp(['Aligning session #' num2str(registration_order(n)) ':'])
        if strcmp(alignment_type,'Translations and Rotations')
            if use_parallel_processing
                disp('Checking for rotations')
                temp_correlations_vector=zeros(1,length(possible_rotations));
                reference_centroid_projections_corrected=centroid_projections_corrected{reference_session_index};
                temp_centroid_projections_corrected=centroid_projections_corrected{registration_order(n)};
                parfor r=1:length(possible_rotations)
                    rotated_image=rotate_image_interp(temp_centroid_projections_corrected,possible_rotations(r),[0 0],center_of_FOV);
                    cross_corr=normxcorr2(reference_centroid_projections_corrected,rotated_image);
                    temp_correlations_vector(r)=max(max(cross_corr));
                end
                if max(temp_correlations_vector)<sufficient_correlation_centroids
                    reference_footprints_projections_corrected=footprints_projections_corrected{reference_session_index};
                    temp_footprints_projections_corrected=footprints_projections_corrected{registration_order(n)};
                    parfor r=1:length(possible_rotations)
                        rotated_image=rotate_image_interp(temp_footprints_projections_corrected,possible_rotations(r),[0 0],center_of_FOV);
                        cross_corr=normxcorr2(reference_footprints_projections_corrected,rotated_image);
                        temp_correlations_vector(r)=max(max(cross_corr));
                    end
                end
                [~,ind_best_rotation]=max(temp_correlations_vector);
                
                % finding the best rotation with a gaussian fit:
                rotation_range_to_check=5; % range in degrees to check for the gaussian fit
                normalized_rotation_range_to_check=rotation_range_to_check/rotation_step;
                rotation_range=round(normalized_rotation_range_to_check);
                
                % zero padding:
                if ind_best_rotation>rotation_range && ind_best_rotation<=length(possible_rotations)-rotation_range
                    localized_max_correlation=temp_correlations_vector(ind_best_rotation-rotation_range:ind_best_rotation+rotation_range);
                elseif ind_best_rotation<=rotation_range
                    zero_padding_size=rotation_range-ind_best_rotation+1;
                    localized_max_correlation=[zeros(1,zero_padding_size) , temp_correlations_vector(1:ind_best_rotation+rotation_range)];
                elseif ind_best_rotation>length(possible_rotations)-rotation_range
                    zero_padding_size=rotation_range-(length(possible_rotations)-ind_best_rotation);
                    localized_max_correlation=[temp_correlations_vector(ind_best_rotation-rotation_range:end), zeros(1,zero_padding_size)];
                end
                normalized_localized_max_correlation=localized_max_correlation-min(localized_max_correlation); % transform to zero basline
                sigma_0=0.1*rotation_range;
                [~,best_rotation_temp]=gaussfit(-rotation_range:rotation_range,normalized_localized_max_correlation./sum(normalized_localized_max_correlation),sigma_0,0);
                best_rotation=possible_rotations(ind_best_rotation)+best_rotation_temp;
            else
                display_progress_bar('Checking for rotations: ',false)
                temp_correlations_vector=zeros(1,length(possible_rotations));
                for r=1:length(possible_rotations)
                    display_progress_bar(100*(r)/length(possible_rotations),false)
                    rotated_image=rotate_image_interp(centroid_projections_corrected{registration_order(n)},possible_rotations(r),[0 0],center_of_FOV);
                    cross_corr=normxcorr2(centroid_projections_corrected{reference_session_index},rotated_image);
                    temp_correlations_vector(r)=max(max(cross_corr));
                end
                % finding the best rotation with a gaussian fit:
                [~,ind_best_rotation]=max(temp_correlations_vector);
                rotation_range_to_check=5; % range in degrees to check for the gaussian fit
                normalized_rotation_range_to_check=rotation_range_to_check/rotation_step;
                rotation_range=round(normalized_rotation_range_to_check);
                
                % zero padding:
                if ind_best_rotation>rotation_range && ind_best_rotation<=length(possible_rotations)-rotation_range
                    localized_max_correlation=temp_correlations_vector(ind_best_rotation-rotation_range:ind_best_rotation+rotation_range);
                elseif ind_best_rotation<=rotation_range
                    zero_padding_size=rotation_range-ind_best_rotation+1;
                    localized_max_correlation=[zeros(1,zero_padding_size) , temp_correlations_vector(1:ind_best_rotation+rotation_range)];
                elseif ind_best_rotation>length(possible_rotations)-rotation_range
                    zero_padding_size=rotation_range-(length(possible_rotations)-ind_best_rotation);
                    localized_max_correlation=[temp_correlations_vector(ind_best_rotation-rotation_range:end), zeros(1,zero_padding_size)];
                end
                normalized_localized_max_correlation=localized_max_correlation-min(localized_max_correlation); % transform to zero basline
                sigma_0=0.1*rotation_range;
                [~,best_rotation_temp]=gaussfit(-rotation_range:rotation_range,normalized_localized_max_correlation./sum(normalized_localized_max_correlation),sigma_0,0);
                best_rotation=possible_rotations(ind_best_rotation)+best_rotation_temp;
                display_progress_bar(' done',false);
            end
            rotation_vector(n)=best_rotation;
            if abs(best_rotation)>minimal_rotation
                rotated_projections=rotate_image_interp(footprints_projections_corrected{registration_order(n)}',-best_rotation,[0 0],center_of_FOV);
                overlapping_area_temp=rotate_image_interp(overlapping_area_all_sessions(:,:,n)',-best_rotation,[0 0],center_of_FOV)';
                all_rotated_projections{registration_order(n)}=rotated_projections';
                theta=best_rotation*pi/180;
                transformation=[cos(theta) -sin(theta) ; sin(theta) cos(theta)]';
                trans_inv=transformation^-1;
                centroids_temp=trans_inv*(centroid_locations_corrected{registration_order(n)}-repmat([center_of_FOV(1) ;center_of_FOV(2)],1,size(centroid_locations_corrected{registration_order(n)},1))')'+repmat([center_of_FOV(1) ;center_of_FOV(2)],1,size(centroid_locations_corrected{registration_order(n)},1));
            else
                all_rotated_projections{registration_order(n)}=footprints_projections_corrected{registration_order(n)};
                centroids_temp=centroid_locations_corrected{registration_order(n)}';
            end
            centroid_locations_rotated{registration_order(n)}=centroids_temp';
            this_session_centroids=centroid_locations_rotated{registration_order(n)};
            number_of_cells=size(this_session_centroids,1);
            this_session_spatial_footprints=spatial_footprints_corrected{registration_order(n)};
            normalized_centroids=zeros(size(this_session_spatial_footprints));
            for k=1:number_of_cells
                if round(this_session_centroids(k,2))>1.5 && round(this_session_centroids(k,1))>1.5 && round(this_session_centroids(k,2))<size(normalized_centroids,2)-1 && round(this_session_centroids(k,1))<size(normalized_centroids,3)-1
                    normalized_centroids(k,round(this_session_centroids(k,2))-1:round(this_session_centroids(k,2))+1,round(this_session_centroids(k,1))-1:round(this_session_centroids(k,1))+1)=1/4;
                    normalized_centroids(k,round(this_session_centroids(k,2))-1:round(this_session_centroids(k,2))+1,round(this_session_centroids(k,1)))=1/2;
                    normalized_centroids(k,round(this_session_centroids(k,2)),round(this_session_centroids(k,1))-1:round(this_session_centroids(k,1))+1)=1/2;
                    normalized_centroids(k,round(this_session_centroids(k,2)),round(this_session_centroids(k,1)))=1;
                elseif round(this_session_centroids(k,2))>0.5 && round(this_session_centroids(k,1))>0.5 && round(this_session_centroids(k,2))<size(normalized_centroids,2) && round(this_session_centroids(k,1))<size(normalized_centroids,3)
                    normalized_centroids(k,round(this_session_centroids(k,2)),round(this_session_centroids(k,1)))=1;
                end
            end
            centroid_projections_rotated{registration_order(n)}=squeeze(sum(normalized_centroids,1));
        end
        
        % Finding translations with subpixel resolution:
        if strcmp(alignment_type,'Translations and Rotations')
            full_FOV_correlation=normxcorr2(all_rotated_projections{reference_session_index},all_rotated_projections{registration_order(n)});
            cross_corr_cent=normxcorr2(centroid_projections_rotated{reference_session_index},centroid_projections_rotated{registration_order(n)});
            if max(max(cross_corr_cent))<sufficient_correlation_centroids
                cross_corr_cent=normxcorr2(all_rotated_projections{reference_session_index},all_rotated_projections{registration_order(n)});
            end
        else
            full_FOV_correlation=normxcorr2(footprints_projections_corrected{reference_session_index},footprints_projections_corrected{registration_order(n)});
            cross_corr_cent=normxcorr2(centroid_projections_corrected{reference_session_index},centroid_projections_corrected{registration_order(n)});
            if max(max(cross_corr_cent))<sufficient_correlation_centroids
                cross_corr_cent=normxcorr2(footprints_projections_corrected{reference_session_index},footprints_projections_corrected{registration_order(n)});
            end
        end
        
        cross_corr_size=size(cross_corr_cent);
        partial_FOV_correlation=full_FOV_correlation(round(cross_corr_size(1)/2-cross_corr_size(1)/6):round(cross_corr_size(1)/2+cross_corr_size(1)/6)...
            ,round(cross_corr_size(2)/2-cross_corr_size(2)/6):round(cross_corr_size(2)/2+cross_corr_size(2)/6));
        maximal_cross_correlation(n)=max(max(partial_FOV_correlation));
        
        cross_corr_partial=cross_corr_cent(round(cross_corr_size(1)/2-cross_corr_size(1)/6):round(cross_corr_size(1)/2+cross_corr_size(1)/6)...
            ,round(cross_corr_size(2)/2-cross_corr_size(2)/6):round(cross_corr_size(2)/2+cross_corr_size(2)/6));
        [~,x_ind]=max(max(cross_corr_partial));
        if strcmp(alignment_type,'Translations and Rotations')
            best_rotations(registration_order(n))=best_rotation;
        end
        [~,y_ind]=max(cross_corr_partial(:,x_ind));
        
        % finding the best translation with a gaussian fit:
        gaussian_radius=round(1.5*normalized_typical_cell_size);
        temp_corr_x=cross_corr_partial(y_ind-1:y_ind+1,x_ind-gaussian_radius:x_ind+gaussian_radius);
        sigma_0=normalized_typical_cell_size/5;
        [~,mu_x_1]=gaussfit(-gaussian_radius:gaussian_radius,temp_corr_x(1,:)./sum(temp_corr_x(1,:)),sigma_0,0);
        [~,mu_x_2]=gaussfit(-gaussian_radius:gaussian_radius,temp_corr_x(2,:)./sum(temp_corr_x(2,:)),sigma_0,0);
        [~,mu_x_3]=gaussfit(-gaussian_radius:gaussian_radius,temp_corr_x(3,:)./sum(temp_corr_x(2,:)),sigma_0,0);
        sub_x=mean([mu_x_1 , mu_x_2 , mu_x_2 , mu_x_3]);
        temp_corr_y=cross_corr_partial(y_ind-gaussian_radius:y_ind+gaussian_radius,x_ind-1:x_ind+1);
        [~,mu_y_1]=gaussfit(-gaussian_radius:gaussian_radius,temp_corr_y(:,1)./sum(temp_corr_y(:,1)),sigma_0,0);
        [~,mu_y_2]=gaussfit(-gaussian_radius:gaussian_radius,temp_corr_y(:,2)./sum(temp_corr_y(:,2)),sigma_0,0);
        [~,mu_y_3]=gaussfit(-gaussian_radius:gaussian_radius,temp_corr_y(:,3)./sum(temp_corr_y(:,3)),sigma_0,0);
        sub_y=mean([mu_y_1 , mu_y_2 , mu_y_2 , mu_y_3]);
        x_ind=x_ind+round(cross_corr_size(2)/2-cross_corr_size(2)/6)-1;
        y_ind=y_ind+round(cross_corr_size(1)/2-cross_corr_size(1)/6)-1;
        if abs(sub_x)>1
            warning(['X axis sub-pixel correction was ' num2str(round(100*sub_x)/100) ' for session ' num2str(registration_order(n))])
            warndlg(['X axis sub-pixel correction was ' num2str(round(100*sub_x)/100) ' for session ' num2str(registration_order(n))])
            x_ind_sub=x_ind;
        else
            x_ind_sub=x_ind+sub_x;
        end
        if abs(sub_y)>1
            warning(['Y axis sub-pixel correction was ' num2str(round(100*sub_y)/100) ' for session ' num2str(registration_order(n))])
            warndlg(['Y axis sub-pixel correction was ' num2str(round(100*sub_y)/100) ' for session ' num2str(registration_order(n))])
            y_ind_sub=y_ind;
        else
            y_ind_sub=y_ind+sub_y;
        end
        best_x_translations(registration_order(n))=(x_ind_sub-adjusted_x_size);
        best_y_translations(registration_order(n))=(y_ind_sub-adjusted_y_size);
        
        % aligning projections and centroid locations:
        if maximal_cross_correlation(n)>median(partial_FOV_correlation(:))+sufficient_correlation_footprints
            if strcmp(alignment_type,'Translations and Rotations')
                untranslated_footprints_projections=all_rotated_projections{registration_order(n)};
                untranslated_centroid_projections=centroid_projections_rotated{registration_order(n)};
                centroids_temp=centroid_locations_rotated{registration_order(n)};
            else
                untranslated_footprints_projections=footprints_projections_corrected{registration_order(n)};
                untranslated_centroid_projections=centroid_projections_corrected{registration_order(n)};
                centroids_temp=centroid_locations_corrected{registration_order(n)};
            end
            centroids_temp(:,1)=centroids_temp(:,1)-(x_ind_sub-adjusted_x_size);
            centroids_temp(:,2)=centroids_temp(:,2)-(y_ind_sub-adjusted_y_size);
            centroid_locations_corrected{registration_order(n)}=centroids_temp;
            translated_projections=translate_projections(untranslated_footprints_projections',[y_ind_sub-adjusted_y_size x_ind_sub-adjusted_x_size]);
            translated_centroid_projections=translate_projections(untranslated_centroid_projections',[y_ind_sub-adjusted_y_size x_ind_sub-adjusted_x_size]);
            translated_overlapping_area=translate_projections(overlapping_area_temp',[y_ind_sub-adjusted_y_size x_ind_sub-adjusted_x_size]);
            translated_overlapping_area=translated_overlapping_area';
            footprints_projections_corrected{registration_order(n)}=translated_projections';
            centroid_projections_corrected{registration_order(n)}=translated_centroid_projections';
            
            % Rotating/translating each spatial footprint:
            unaligned_spatial_footprints=spatial_footprints_corrected{registration_order(n)};
            number_of_cells=size(unaligned_spatial_footprints,1);
            if strcmp(alignment_type,'Translations and Rotations') && abs(best_rotation)>minimal_rotation
                % rotating cells
                unrotated_spatial_footprints=unaligned_spatial_footprints;
                rotated_translated_spatial_footprints=zeros(number_of_cells,adjusted_y_size,adjusted_x_size);
                unrotated_centroid_locations=centroid_locations{registration_order(n)};
                if use_parallel_processing
                    disp('Rotating and translating spatial footprints')
                    parfor m=1:number_of_cells
                        unrotated_spatial_footprint=squeeze(unrotated_spatial_footprints(m,:,:));
                        unrotated_centroid=unrotated_centroid_locations(m,:);
                        rotated_translated_spatial_footprint=rotate_spatial_footprint(unrotated_spatial_footprint',-best_rotation,[y_ind_sub-adjusted_y_size x_ind_sub-adjusted_x_size],center_of_FOV,unrotated_centroid,microns_per_pixel);
                        rotated_translated_spatial_footprints(m,:,:)=rotated_translated_spatial_footprint';
                    end
                else
                    display_progress_bar('Rotating spatial footprints: ',false)
                    for m=1:size(centroid_locations_corrected{registration_order(n)},1)
                        display_progress_bar(100*(m/size(unrotated_centroid_locations,1)),false)
                        unrotated_spatial_footprint=squeeze(unrotated_spatial_footprints(m,:,:));
                        unrotated_centroid=unrotated_centroid_locations(m,:);
                        rotated_translated_spatial_footprint=rotate_spatial_footprint(unrotated_spatial_footprint',-best_rotation,[y_ind_sub-adjusted_y_size x_ind_sub-adjusted_x_size],center_of_FOV,unrotated_centroid,microns_per_pixel);
                        rotated_translated_spatial_footprints(m,:,:)=rotated_translated_spatial_footprint';
                    end
                    display_progress_bar(' done',false)
                end
                aligned_spatial_footprints=rotated_translated_spatial_footprints;
            else % no rotations - translating cells
                if strcmp(alignment_type,'Translations and Rotations') && abs(best_rotation)<=minimal_rotation
                    disp('No rotations required') % less than the minimal rotation that justifies rotating each cell - translating cells
                end
                untranslated_spatial_footprints=unaligned_spatial_footprints;
                translated_spatial_footprints=zeros(number_of_cells,adjusted_y_size,adjusted_x_size);
                if strcmp(alignment_type,'Non-rigid') % Non-rigid alignment:
                    untranslated_centroid_locations=centroid_locations_corrected{registration_order(n)};
                else
                    untranslated_centroid_locations=centroid_locations{registration_order(n)};
                end
                if use_parallel_processing
                    disp('Translating spatial footprints')
                    parfor k=1:number_of_cells
                        untranslated_spatial_footprint=squeeze(untranslated_spatial_footprints(k,:,:));
                        untranslated_centroid=untranslated_centroid_locations(k,:);
                        translated_spatial_footprint=translate_spatial_footprint(untranslated_spatial_footprint',[y_ind_sub-adjusted_y_size x_ind_sub-adjusted_x_size],untranslated_centroid,microns_per_pixel);
                        translated_spatial_footprints(k,:,:)=translated_spatial_footprint';
                    end
                else
                    display_progress_bar('Translating spatial footprints: ',false)
                    for k=1:number_of_cells
                        display_progress_bar(100*(k/size(untranslated_centroid_locations,1)),false)
                        untranslated_spatial_footprint=squeeze(untranslated_spatial_footprints(k,:,:));
                        untranslated_centroid=untranslated_centroid_locations(k,:);
                        translated_spatial_footprint=translate_spatial_footprint(untranslated_spatial_footprint',[y_ind_sub-adjusted_y_size x_ind_sub-adjusted_x_size],untranslated_centroid,microns_per_pixel);
                        translated_spatial_footprints(k,:,:)=translated_spatial_footprint';
                    end
                    display_progress_bar(' done',false)
                end
                aligned_spatial_footprints=translated_spatial_footprints;
            end
            spatial_footprints_corrected{registration_order(n)}=aligned_spatial_footprints;
            overlapping_area=overlapping_area.*translated_overlapping_area;
        else % if no appropriate rotations/translations were found
            if strcmp(alignment_type,'Translations and Rotations') % rotating cells
                warning(['Session ' num2str(registration_order(n)) ' does not resemble the reference session - try using non-rigid transformation'])
                warndlg(['Session ' num2str(registration_order(n)) ' does not resemble the reference session - try using non-rigid transformation'])
            elseif strcmp(alignment_type,'Translations') % rotating cells
                warning(['Session ' num2str(registration_order(n)) ' does not resemble the reference session - try using rotations'])
                warndlg(['Session ' num2str(registration_order(n)) ' does not resemble the reference session - try using rotations'])
            else
                warning(['Session ' num2str(registration_order(n)) ' does not resemble the reference session - could not align session'])
                warndlg(['Session ' num2str(registration_order(n)) ' does not resemble the reference session - could not align session'])
            end
        end
    end
    
    best_translations=microns_per_pixel*[best_x_translations ; best_y_translations];
    if strcmp(alignment_type,'Translations and Rotations')
        best_translations=[best_translations ; best_rotations];
    end    
end
