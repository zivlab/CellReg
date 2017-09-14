function [p_same_centers_of_bins,uncertain_fraction_centroid_distances,cdf_p_same_centroid_distances,false_positive_per_distance_threshold,true_positive_per_distance_threshold,varargout]=estimate_registration_accuracy(centroid_distances_model_parameters,p_same_certainty_threshold,neighbors_centroid_distances,centroid_distances_model_same_cells,centroid_distances_model_different_cells,p_same_given_centroid_distance,centers_of_bins,varargin)
% This function estiamtes the fraction of cell-pairs with less than 95%
% confidence of registration, and the false positive and false negative
% rates.

% Inputs:
% 1. centroid_distances_model_parameters
% 2. p_same_certainty_threshold
% 3. neighbors_centroid_distances
% 4. centroid_distances_model_same_cells
% 5. centroid_distances_model_different_cells
% 6. p_same_given_centroid_distance
% 7. centers_of_bins
% 8. varargin
%   8{1}. spatial_correlations_model_parameters
%   8{2}. neighbors_spatial_correlations
%   8{3}. spatial_correlations_model_same_cells
%   8{4}. spatial_correlations_model_different_cells
%   8{5}. p_same_given_spatial_correlation

% Outputs:
% 1. p_same_centers_of_bins - evenly spaced between 0-1
% 2. uncertain_fraction_centroid_distances - cell pairs with 0.05<P_Same<0.95
% 3. cdf_p_same_centroid_distances
% 4. false_positive_per_distance_threshold
% 5. true_positive_per_distance_threshold
% 6. varargout
%   6{1}. uncertain_fraction_spatial_correlations
%   6{2}. cdf_p_same_centroid_distances
%   6{3}. false_positive_per_distance_threshold
%   6{4}. true_positive_per_distance_threshold

p_same_certainty_threshold=1-p_same_certainty_threshold; % it should be 1 minus the value

number_of_p_same_bins=1000;

p_same_centers_of_bins_temp=linspace(0,1,2*number_of_p_same_bins+1);
p_same_centers_of_bins=p_same_centers_of_bins_temp(2:2:end);
if ~isempty(varargin)
    cumulative_true_positive_per_correlation_threshold=zeros(1,length(p_same_centers_of_bins));
    cumulative_false_positive_per_correlation_threshold=zeros(1,length(p_same_centers_of_bins));
    cdf_p_same_spatial_correlations=zeros(1,length(p_same_centers_of_bins));
end
cumulative_true_positive_per_distance_threshold=zeros(1,length(p_same_centers_of_bins));
cumulative_false_positive_per_distance_threshold=zeros(1,length(p_same_centers_of_bins));
cdf_p_same_centroid_distances=zeros(1,length(p_same_centers_of_bins));

step=p_same_centers_of_bins(2)-p_same_centers_of_bins(1);
if ~isempty(varargin)
    spatial_correlations_model_parameters=varargin{1};
    neighbors_spatial_correlations=varargin{2};
    spatial_correlations_model_same_cells=varargin{3};
    spatial_correlations_model_different_cells=varargin{4};
    p_same_given_spatial_correlation=varargin{5};
    p_same_weight_spatial_correlations=spatial_correlations_model_parameters(1);
    [spatial_correlations_distribution,~]=hist(neighbors_spatial_correlations,centers_of_bins{2});
end
[centroid_distances_distribution,~]=hist(neighbors_centroid_distances,centers_of_bins{1});
p_same_weight_centroid_distances=centroid_distances_model_parameters(1);
for n=1:length(p_same_centers_of_bins) 
    if ~isempty(varargin) % if the spatial correlations model is also used
        spatial_correlations_model_same_cells_for_ROC=spatial_correlations_model_same_cells./sum(spatial_correlations_model_same_cells);
        spatial_correlations_model_different_cells_for_ROC=spatial_correlations_model_different_cells./sum(spatial_correlations_model_different_cells);
        cumulative_true_positive_per_correlation_threshold(n)=1-p_same_weight_spatial_correlations*sum(spatial_correlations_model_same_cells_for_ROC(1-p_same_given_spatial_correlation>=p_same_centers_of_bins(n)-step/2))...
            /(p_same_weight_spatial_correlations*sum(spatial_correlations_model_same_cells_for_ROC(1-p_same_given_spatial_correlation>=p_same_centers_of_bins(n)-step/2))+(1-p_same_weight_spatial_correlations)*sum(spatial_correlations_model_different_cells_for_ROC(1-p_same_given_spatial_correlation>=p_same_centers_of_bins(n)-step/2)));
        cumulative_false_positive_per_correlation_threshold(n)=(1-p_same_weight_spatial_correlations)*sum(spatial_correlations_model_different_cells_for_ROC(1-p_same_given_spatial_correlation<=p_same_centers_of_bins(n)+step/2))...
            /((1-p_same_weight_spatial_correlations)*sum(spatial_correlations_model_different_cells_for_ROC(1-p_same_given_spatial_correlation<=p_same_centers_of_bins(n)+step/2))+p_same_weight_spatial_correlations*sum(spatial_correlations_model_same_cells_for_ROC(1-p_same_given_spatial_correlation<=p_same_centers_of_bins(n)+step/2)));
        cdf_p_same_spatial_correlations(n)=sum(spatial_correlations_distribution(p_same_given_spatial_correlation<=p_same_centers_of_bins(n)+step/2 & p_same_given_spatial_correlation>=p_same_centers_of_bins(n)-step/2));
    end
    centroid_distances_model_same_cells_for_ROC=centroid_distances_model_same_cells./sum(centroid_distances_model_same_cells);
    centroid_distances_model_different_cells_for_ROC=centroid_distances_model_different_cells./sum(centroid_distances_model_different_cells);
    cumulative_true_positive_per_distance_threshold(n)=1-p_same_weight_centroid_distances*sum(centroid_distances_model_same_cells_for_ROC(1-p_same_given_centroid_distance>=p_same_centers_of_bins(n)-step/2))...
        /(p_same_weight_centroid_distances*sum(centroid_distances_model_same_cells_for_ROC(1-p_same_given_centroid_distance>=p_same_centers_of_bins(n)-step/2))+(1-p_same_weight_centroid_distances)*sum(centroid_distances_model_different_cells_for_ROC(1-p_same_given_centroid_distance>=p_same_centers_of_bins(n)-step/2)));
    cumulative_false_positive_per_distance_threshold(n)=(1-p_same_weight_centroid_distances)*sum(centroid_distances_model_different_cells_for_ROC(1-p_same_given_centroid_distance<=p_same_centers_of_bins(n)+step/2))...
        /((1-p_same_weight_centroid_distances)*sum(centroid_distances_model_different_cells_for_ROC(1-p_same_given_centroid_distance<=p_same_centers_of_bins(n)+step/2))+p_same_weight_centroid_distances*sum(centroid_distances_model_same_cells_for_ROC(1-p_same_given_centroid_distance<=p_same_centers_of_bins(n)+step/2)));
    cdf_p_same_centroid_distances(n)=sum(centroid_distances_distribution(p_same_given_centroid_distance<=p_same_centers_of_bins(n)+step/2 & p_same_given_centroid_distance>=p_same_centers_of_bins(n)-step/2));
end
false_positive_per_distance_threshold=diff([0, cumulative_false_positive_per_distance_threshold]);
true_positive_per_distance_threshold=diff([0, cumulative_true_positive_per_distance_threshold]);

if ~isempty(varargin) % if the spatial correlations model is also used:
    false_positive_per_correlation_threshold=diff([0, cumulative_false_positive_per_correlation_threshold]);
    true_positive_per_correlation_threshold=diff([0, cumulative_true_positive_per_correlation_threshold]);
    uncertain_fraction_spatial_correlations=(sum(cdf_p_same_spatial_correlations(1-p_same_centers_of_bins>p_same_certainty_threshold & 1-p_same_centers_of_bins<1-p_same_certainty_threshold)))/sum(cdf_p_same_spatial_correlations);
    varargout{1}=uncertain_fraction_spatial_correlations;
    varargout{2}=cdf_p_same_spatial_correlations;
    varargout{3}=false_positive_per_correlation_threshold;
    varargout{4}=true_positive_per_correlation_threshold;
end
uncertain_fraction_centroid_distances=(sum(cdf_p_same_centroid_distances(1-p_same_centers_of_bins>p_same_certainty_threshold & 1-p_same_centers_of_bins<1-p_same_certainty_threshold)))/sum(cdf_p_same_centroid_distances);

end

