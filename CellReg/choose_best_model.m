function [best_model_string]=choose_best_model(MSE_centroid_distances_model,centroid_distances_model_same_cells,centroid_distances_model_different_cells,p_same_given_centroid_distance,varargin)
% This function uses the registration uncertainty levels and mean squared
% error for the fit of each model to choose the best model

% Inputs:
% 1. MSE_centroid_distances_model
% 2. centroid_distances_model_same_cells
% 3. centroid_distances_model_different_cells
% 4. p_same_given_centroid_distance
% 5. varargin
%   5{1}. MSE_spatial_correlations_model
%   5{2}. spatial_correlations_model_same_cells
%   5{3}. spatial_correlations_model_different_cells
%   5{4}. p_same_given_spatial_correlation

% Outputs:
% 1. best_model_string

MSE_spatial_correlations_model=varargin{1};
spatial_correlations_model_same_cells=varargin{2};
spatial_correlations_model_different_cells=varargin{3};
p_same_given_spatial_correlation=varargin{4};

[~,ind_05_correlation]=min(abs(0.5-(p_same_given_spatial_correlation)));
false_positive_correlation=sum(spatial_correlations_model_different_cells(ind_05_correlation:end))/sum(spatial_correlations_model_different_cells);
false_negative_correlation=sum(spatial_correlations_model_same_cells(1:ind_05_correlation))/sum(spatial_correlations_model_same_cells);
[~,ind_05_distance]=min(abs(0.5-(p_same_given_centroid_distance)));
false_positive_distance=sum(centroid_distances_model_different_cells(1:ind_05_correlation))/sum(centroid_distances_model_different_cells);
false_negative_distance=sum(centroid_distances_model_same_cells(ind_05_correlation:end))/sum(centroid_distances_model_same_cells);

cost_vector=[false_positive_distance+false_negative_distance+MSE_centroid_distances_model,false_positive_correlation+false_negative_correlation+MSE_spatial_correlations_model];
if MSE_centroid_distances_model>0.1
    warning('There is large discrepancy between the centroid distances model and the data')
end
if MSE_spatial_correlations_model>0.1
    warning('There is large discrepancy between the spatial correlations model and the data')
end
[~,best_model]=min(cost_vector);
if best_model==1
    disp('The centroid distances model is best suited for the data');
    best_model_string='Centroid distance';
else
    disp('The spatial correlations model is best suited for the data');
    best_model_string='Spatial correlation';
end

end

