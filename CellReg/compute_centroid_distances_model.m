function [centroid_distances_model_parameters,p_same_given_centroid_distance,centroid_distances_distribution,centroid_distances_model_same_cells,centroid_distances_model_different_cells,centroid_distances_model_weighted_sum,MSE_centroid_distances_model,centroid_distance_intersection]=compute_centroid_distances_model(neighbors_centroid_distances,microns_per_pixel,centers_of_bins)
% This function recieves the distribution of centroid distances for all
% neighboring cells-pairs across sessions, and computes a probabilistic
% model by finding the weighted sum of same cells and different cells that
% best fit the data. 
% Different cells are modeled by a multipication of a linear function by a
% sigmoid function. Same cells are models as a lognormal distribution.

% Inputs:
% 1. neighbors_centroid_distances - all neighboring cell-pairs
% 2. microns_per_pixel
% 3. maximal_distance
% 4. centers_of_bins

% Outputs:
% 1. centroid_distances_model_parameters
% 2. p_same_given_centroid_distance
% 3. centroid_distances_distribution
% 4. centroid_distances_model_same_cells
% 5. centroid_distances_model_different_cells
% 6. centroid_distances_model_weighted_sum
% 7. MSE_centroid_distances_model
% 8. centroid_distance_intersection - the distance for which same cells are
% equal to different cells (P_same=0.5)

number_of_bins=length(centers_of_bins{1});
centroid_distances_centers=centers_of_bins{1};

[centroid_distances_distribution,~]=hist(neighbors_centroid_distances,centroid_distances_centers);
centroid_distances_distribution=centroid_distances_distribution./sum(centroid_distances_distribution)*(number_of_bins/(microns_per_pixel*(centroid_distances_centers(2)-centroid_distances_centers(1))+microns_per_pixel*(centroid_distances_centers(end)-centroid_distances_centers(1))));

% finding initial parameters for lsqcurvefit function:
maximal_distance_to_fit=9; % this distance is used only to find the initial paramters
data_to_fit=neighbors_centroid_distances(microns_per_pixel*neighbors_centroid_distances<maximal_distance_to_fit);
parmhat=lognfit(data_to_fit);
optimal_delta=length(data_to_fit)/length(neighbors_centroid_distances);
p_0=optimal_delta;
c_0=6;
a_0=1;
b_0=(centroid_distances_distribution(end)-centroid_distances_distribution(round(number_of_bins/2)))/(microns_per_pixel*(centroid_distances_centers(end)-centroid_distances_centers(round(number_of_bins/2))));
b_0=b_0/(1-p_0);
initial_parameters=[p_0 parmhat a_0 c_0 b_0];
F = @(x,xdata)...
    x(1)*(1./(xdata.*x(3).*sqrt(2*pi)).*exp(-(log(xdata)-x(2)).^2./(2*x(3)^2))...
    + (1-x(1)).*x(6).*xdata./(1+exp(-x(4).*(xdata-x(5)))));

lb = [0 -inf 0 0 0 0];
ub = [1 Inf Inf Inf Inf inf];
options = statset('MaxIter',1000, 'MaxFunEvals',2000);

% finding the parameters that best fit the data:
centroid_distances_model_parameters=lsqcurvefit(F,initial_parameters,microns_per_pixel*centroid_distances_centers,centroid_distances_distribution,lb,ub,options);

% calculating the distribution for same cells:
centroid_distances_model_same_cells=lognpdf(microns_per_pixel*centroid_distances_centers,centroid_distances_model_parameters(2),centroid_distances_model_parameters(3));
centroid_distances_model_same_cells=centroid_distances_model_same_cells./sum(centroid_distances_model_same_cells)*(number_of_bins/(microns_per_pixel*(centroid_distances_centers(2)-centroid_distances_centers(1))+microns_per_pixel*(centroid_distances_centers(end)-centroid_distances_centers(1))));
% calculating the distribution for different cells:
centroid_distances_model_different_cells=centroid_distances_model_parameters(6)*microns_per_pixel*centroid_distances_centers./(1+exp(-centroid_distances_model_parameters(4)*(microns_per_pixel*centroid_distances_centers-centroid_distances_model_parameters(5))));
centroid_distances_model_different_cells=centroid_distances_model_different_cells./sum(centroid_distances_model_different_cells)*(number_of_bins/(microns_per_pixel*(centroid_distances_centers(2)-centroid_distances_centers(1))+microns_per_pixel*(centroid_distances_centers(end)-centroid_distances_centers(1))));
% calculating the weighted sum:
centroid_distances_model_weighted_sum=centroid_distances_model_parameters(1)*centroid_distances_model_same_cells+(1-centroid_distances_model_parameters(1))*centroid_distances_model_different_cells;
% findind the intersection between same cells and different cells:
index_range_of_intersection=find(centroid_distances_centers>1/microns_per_pixel & centroid_distances_centers<10/microns_per_pixel);
[~,index_of_intersection]=min(abs(centroid_distances_model_parameters(1)*centroid_distances_model_same_cells(index_range_of_intersection)-(1-centroid_distances_model_parameters(1))*centroid_distances_model_different_cells(index_range_of_intersection)));
centroid_distance_intersection=round(100*microns_per_pixel*centroid_distances_centers(index_of_intersection+index_range_of_intersection(1)-1))/100;

% calculating the discrepancy of the model (normalized mean sqaured error)
MSE_centroid_distances_model=sum(abs(((centroid_distances_distribution-centroid_distances_model_weighted_sum))*(microns_per_pixel*(centroid_distances_centers(2)-centroid_distances_centers(1))+microns_per_pixel*(centroid_distances_centers(end)-centroid_distances_centers(1)))/number_of_bins))/2;

% calculating the P_same of the model:
p_same_given_centroid_distance=centroid_distances_model_parameters(1).*centroid_distances_model_same_cells./(centroid_distances_model_parameters(1).*centroid_distances_model_same_cells+(1-centroid_distances_model_parameters(1)).*centroid_distances_model_different_cells);
p_same_given_centroid_distance(1)=p_same_given_centroid_distance(2); % avoid p_same going to 0 at 0 distance because of the lognormal distribution

end

