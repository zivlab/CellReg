function [best_model_string]=choose_best_model(uncertain_fraction_centroid_distances,MSE_centroid_distances_model,imaging_technique,varargin)
% This function uses the registration uncertainty levels and mean squared
% error for the fit of each model to choose the best model

% Inputs:
% 1. uncertain_fraction_centroid_distances
% 2. MSE_centroid_distances_model
% 3. imaging_technique
% 4. varargin
%   4{1}. uncertain_fraction_spatial_correlations
%   4{2}. MSE_spatial_correlations_model
 
% Outputs:
% 1. best_model_string

if strcmp(imaging_technique,'one_photon');
    uncertain_fraction_spatial_correlations=varargin{1};
    MSE_spatial_correlations_model=varargin{2};
    cost_vector=[uncertain_fraction_centroid_distances+MSE_centroid_distances_model,uncertain_fraction_spatial_correlations+MSE_spatial_correlations_model];
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
else
    best_model_string='Centroid distance';
    if MSE_centroid_distances_model>0.1        
        warning('There is large discrepancy between the centroid distances model and the data')
    end
end

end

