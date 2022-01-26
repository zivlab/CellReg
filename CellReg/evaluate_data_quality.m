function [all_projections_correlations,number_of_cells_per_session]=evaluate_data_quality(spatial_footprints,centroid_projections_corrected,footprints_projections_corrected,maximal_cross_correlation,best_translations,reference_session_index,sufficient_correlation,alignment_type)
% This function assesses the quality of the data and its suitabilty for
% longitudinal analysis.

% Inputs:
% 1. spatial_footprints
% 2. centroid_projections_corrected
% 3. microns_per_pixel
% 4. maximal_cross_correlation - between each session and the reference
% 5. best_translations
% 6. reference_session_index
% 7. sufficient_correlation % smaller correlation imply different optical section or high noise levels
% 8. alignment_type

% Outputs:
% 1. all_projections_correlations - correlations for all pairs of sessions
% 2. number_of_cells_per_session

large_rotation=10; % in degrees
large_translation=50; % in microns
abnormal_number_of_cells_ratio=1.5;
number_of_sessions=size(centroid_projections_corrected,2);
registration_order=setdiff(1:number_of_sessions,reference_session_index);

% calculating correlations between the centroid projections for all pairs of sessions:
if number_of_sessions>2
    all_projections_correlations=zeros(number_of_sessions,number_of_sessions);
    for n=1:number_of_sessions
        all_projections_correlations(n,n)=1;
        for k=n+1:number_of_sessions
            all_projections_correlations(n,k)=corr2(footprints_projections_corrected{n},footprints_projections_corrected{k});            
            all_projections_correlations(k,n)=all_projections_correlations(n,k);
        end
    end
else
    all_projections_correlations=maximal_cross_correlation;
end

number_of_cells_per_session=zeros(1,number_of_sessions);
for n=1:number_of_sessions
    number_of_cells_per_session(n)=size(spatial_footprints{n},1);
end
mean_number_of_cells=mean(number_of_cells_per_session);

% warning for unstable preparation:
for n=1:number_of_sessions-1
    if number_of_cells_per_session(n)>abnormal_number_of_cells_ratio*mean_number_of_cells
        warning(['Extermely high number of cells is observed in session number ' num2str(registration_order(n))])
    end
    if number_of_cells_per_session(n)<mean_number_of_cells/abnormal_number_of_cells_ratio
        warning(['Extermely low number of cells is observed in session number ' num2str(registration_order(n))])
    end    
    if best_translations(1,n)>large_translation
        warning([num2str(best_translations(1,n)) ' micron x translation was found for session number ' num2str(registration_order(n))])
    end
    if best_translations(2,n)>large_translation
        warning([num2str(best_translations(1,n)) ' micron y translation was found for session number ' num2str(registration_order(n))])
    end
    if size(best_translations,1)==3
        if abs(best_translations(3,n))>large_rotation
            warning([num2str(best_translations(2,n)) ' degrees rotation was found for session number ' num2str(registration_order(n))])
        end
    end
end
    if maximal_cross_correlation(n)<sufficient_correlation
        if strcmp(alignment_type,'Translations and Rotations')
            warning(['No appropriate translations/rotations were found for session number ' num2str(registration_order(n)) ' - consider using non-rigid transformation'])
        elseif strcmp(alignment_type,'Translations')
            warning(['No appropriate translations were found for session number ' num2str(registration_order(n)) ' - consider using rotations as well'])
        else
            warning(['No appropriate translations were found for session number ' num2str(registration_order(n))])
        end
    end
end

