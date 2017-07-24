function save_log_file(results_directory,file_names,imaging_technique,microns_per_pixel,adjusted_x_size,adjusted_y_size,alignment_type,reference_session,maximal_distance,number_of_bins,initial_registration_type,initial_threshold,registration_approach,model_type,final_threshold,optimal_cell_to_index_map,cell_registered_struct,comments,varargin)
% This function saves into a log-file all the parameters that were used for
% registration and a summary of the registration results.

% Inputs:
% 1. results_directory
% 2. file_names - a cell array
% 3. imaging_technique
% 4. microns_per_pixel
% 5. adjusted_x_size
% 6. adjusted_y_size
% 7. alignment_type
% 8. reference_session
% 9. maximal_distance
% 10. number_of_bins
% 11. initial_registration_type
% 12. initial_threshold
% 13. registration_approach
% 14. model_type
% 15. final_threshold
% 16. optimal_cell_to_index_map
% 17. cell_registered_struct
% 18. comments
% 19. varargin
% 19{1}. uncertain_fraction
% 19{2}. false_positive
% 19{3}. true_positive
% 19{4}. model_MSE

number_of_cells=size(optimal_cell_to_index_map,1);
if ~isempty(varargin)
    uncertain_fraction=varargin{1};
    false_positive=varargin{2};
    true_positive=varargin{3};
    model_MSE=varargin{4};
end

logFile=fopen(fullfile(results_directory,['logFile_' datestr(clock,'yyyymmdd_HHMMss') '.txt']), 'wt' );

% General data parameters:
number_of_sessions=size(file_names,2);
fprintf(logFile,'Sessions list:\n');
for n=1:number_of_sessions
    if n==1
        fprintf(logFile,'-----------------\n');
    end
    fprintf(logFile,'%s',['Session ' num2str(n) ' - ' file_names{1,n}]);
    fprintf(logFile,'\n');
end
fprintf(logFile,'\n\nGeneral data parameters:\n------------------------\n');
if strcmp(imaging_technique,'one_photon')
    fprintf(logFile,'Imaging technique - 1-photon\n');
else
    fprintf(logFile,'Imaging technique - 2-photon\n');
end
fprintf(logFile,['Number of sessions - ' num2str(number_of_sessions) '\n']);
FOV_size=[adjusted_x_size*microns_per_pixel adjusted_y_size*microns_per_pixel];
fprintf(logFile,['FOV size - ' num2str(FOV_size(1)) 'X' num2str(FOV_size(2)) ' [microns] ; (x,y)\n']);
image_size=[adjusted_x_size adjusted_y_size];
fprintf(logFile,['Image size - ' num2str(image_size(1)) 'X' num2str(image_size(2)) ' [pixels] ; (x,y)\n\n']);

% Image alignment parameters:
fprintf(logFile,'\nSession alignment parameters:\n----------------------------------\n');
fprintf(logFile,['Reference session - ' num2str(reference_session) '\n']);
fprintf(logFile,['Alignment Type - ' alignment_type '\n\n']);

% Probabilistic model parameters:
fprintf(logFile,'\nCell registration parameters:\n-----------------------------\n');
fprintf(logFile,['Model maximal distance - ' num2str(maximal_distance) ' [microns]\n']);
if strcmp(registration_approach,'Probabilistic')
    fprintf(logFile,['Number of bins - ' num2str(number_of_bins) '\n']);
end

% Initial alignment parameters:
if strcmp(initial_registration_type,'Spatial correlation')
    correlation_thresh=round(100*initial_threshold)/100;
    fprintf(logFile,['Initial registration type - ' initial_registration_type '\n']);
    fprintf(logFile,['Initial threshold - ' num2str(correlation_thresh) '\n']);
else
    distance_thresh=round(10*initial_threshold)/10;
    fprintf(logFile,['initial registration type - ' initial_registration_type '\n']);
    fprintf(logFile,['Initial threshold - ' num2str(distance_thresh) ' [microns]\n']);
end

% Final alignment parameters:
if strcmp(registration_approach,'Probabilistic')
    num_bins_p_same=length(true_positive);
    decision_thresh=round(100*final_threshold)/100;
    fprintf(logFile,'Registration approach - Probabilistic modeling\n');
    if strcmp(model_type,'Spatial correlation')
        final_type='Spatial correlations';        
    elseif strcmp(model_type,'Centroid distance')
        final_type='Centroid distances';        
    end
    fprintf(logFile,['Final registration type - ' final_type '\n']);
    fprintf(logFile,['P_same threshold - ' num2str(decision_thresh) ' \n\n']);
    uncertain_pairs_fraction=round(100*uncertain_fraction)/100;
    cumsum_false_positive=cumsum(false_positive);
    cumsum_false_split=1-cumsum(true_positive);
    false_positives=round(100*cumsum_false_positive(round(num_bins_p_same*(1-decision_thresh))))/100;
    false_negatives=round(100*cumsum_false_split(round(num_bins_p_same*(1-decision_thresh))))/100;
else
    fprintf(logFile,'Registration approach - Simple threshold\n');
    if strcmp(model_type,'Spatial correlation')
        fprintf(logFile,['Final threshold - ' num2str(final_threshold) ' \n\n']);
    elseif strcmp(model_type,'Centroid distance')
        fprintf(logFile,['Final threshold - ' num2str(final_threshold) ' [microns] \n\n']);
    end
end

% General results:
fprintf(logFile,'\nGeneral results:\n----------------\n');
fprintf(logFile,['Final number of cells - ' num2str(number_of_cells) '\n']);
average_fraction_active_cells=round(100*sum(sum(optimal_cell_to_index_map>0))/number_of_sessions/number_of_cells)/100;
fprintf(logFile,['Average fraction of active cells - ' num2str(average_fraction_active_cells) '\n\n']);

if strcmp(registration_approach,'Probabilistic')
    fprintf(logFile,['Fraction of false positives - ' num2str(false_positives) '\n']);
    fprintf(logFile,['Fraction of false negatives - ' num2str(false_negatives) '\n']);
    fprintf(logFile,['Fraction of uncertain cell-pairs - ' num2str(uncertain_pairs_fraction) '\n']);
end

if strcmp(registration_approach,'Probabilistic')
    register_scores= cell_registered_struct.cell_scores';
    cell_scores_positive=cell_registered_struct.true_positive_scores';
    cell_scores_negative=cell_registered_struct.true_negative_scores';
    cell_scores_exclusive=cell_registered_struct.exclusivity_scores';    
    average_score=round(100*mean(register_scores))/100;
    fprintf(logFile,['\nAverage register score - ' num2str(average_score) '\n']);
    average_score_positive=round(100*mean(cell_scores_positive(cell_scores_positive>=0)))/100;
    fprintf(logFile,['Average true positive score - ' num2str(average_score_positive) '\n']);
    average_score_negative=round(100*mean(cell_scores_negative(cell_scores_negative>=0)))/100;
    fprintf(logFile,['Average true negative score - ' num2str(average_score_negative) '\n']);
    average_score_exclusive=round(100*mean(cell_scores_exclusive(cell_scores_exclusive>=0)))/100;
    fprintf(logFile,['Average exclusivity score - ' num2str(average_score_exclusive) '\n']);
    fprintf(logFile,['\nDiscrepancy of the model - ' num2str(round(100*model_MSE)/100) '\n']);   
end

fprintf(logFile,'\n\nComments:\n---------\n');
num_rows=size(comments,1);
for n=1:num_rows
    fprintf(logFile,comments(n,:));
    fprintf(logFile,'\n');
end
fclose(logFile);

end

