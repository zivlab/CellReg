function [spatial_correlations_model_parameters,p_same_given_spatial_correlation,spatial_correlations_distribution,spatial_correlations_model_same_cells,spatial_correlations_model_different_cells,spatial_correlations_model_weighted_sum,MSE_spatial_correlations_model,spatial_correlation_intersection]=compute_spatial_correlations_model(neighbors_spatial_correlations,centers_of_bins)
% This function recieves the distribution of spatial correlations for all
% neighboring cells-pairs across sessions, and computes a probabilistic
% model by finding the weighted sum of same cells and different cells that
% best fit the data. Same cells are modeled as a lognormal distribution.
% Different cells are modeled as a betha distribution.

% Inputs:
% 1. neighbors_spatial_correlations - all neighboring cell-pairs
% 2. centers_of_bins

% Outputs:
% 1. spatial_correlations_model_parameters
% 2. p_same_given_spatial_correlation
% 3. spatial_correlations_distribution
% 4. spatial_correlations_model_same_cells
% 5. spatial_correlations_model_different_cells,
% 6. spatial_correlations_model_weighted_sum
% 7. MSE_spatial_correlations_model
% 8. spatial_correlation_intersection - the correlation for which same cells are
% equal to different cells (P_same=0.5)

number_of_bins=length(centers_of_bins{2});
spatial_correlations_centers=centers_of_bins{2};

% the maximal distance should be chosen in a way that neighboring cells
% pairs have non-zero overlap
[spatial_correlations_distribution,~]=hist(neighbors_spatial_correlations,spatial_correlations_centers);
if sum(spatial_correlations_distribution(spatial_correlations_centers<0.05))>0.05*sum(spatial_correlations_distribution)
    warning('A good fit might not be attainable because some of the cells seem to be smaller than expected. This could also occur if the provided microns per pixel ratio is incorrect or the maximal distance is too large')
end

% finding initial parameters for mle function:
pdf_betamixture = @(x,p,A1,A2,B1,B2) ...
    p*lognpdf(x,A1,B1) + (1-p)*betapdf(x,A2,B2);
phat_same=lognfit(1-neighbors_spatial_correlations(neighbors_spatial_correlations>=0.7));
phat_diff=betafit(neighbors_spatial_correlations(neighbors_spatial_correlations<0.7));
pinitial_parameters=0.5;
Ainitial_parameters=[phat_same(1) phat_diff(1)];
Binitial_parameters=[phat_same(2) phat_diff(2)];
initial_parameters=[pinitial_parameters Ainitial_parameters Binitial_parameters];
lb = [0 -inf 0 0 0];
ub = [1 Inf Inf Inf Inf];
options = statset('MaxIter',1000, 'MaxFunEvals',2000);

% finding the parameters that best fit the data:
spatial_correlations_model_parameters=mle(1-neighbors_spatial_correlations,'pdf',pdf_betamixture, 'start',initial_parameters, 'lower',lb, 'upper',ub,'options',options);

% calculating the distribution for same cells:
spatial_correlations_model_same_cells=lognpdf(1-spatial_correlations_centers,spatial_correlations_model_parameters(2),spatial_correlations_model_parameters(4));
spatial_correlations_model_same_cells=spatial_correlations_model_same_cells./sum(spatial_correlations_model_same_cells)*(number_of_bins/((spatial_correlations_centers(2)-spatial_correlations_centers(1))+(spatial_correlations_centers(end)-spatial_correlations_centers(1))));
% the same cells model is multiplied by a sigmoid function because 
%the lognormal ditribution goes to infinity but the correlation is bounded:
sigmoid_function=@(x,ac)1./(1+exp(-ac(1)*(x-ac(2)))); % defining the sigmoid function - (sigmf requires Fuzzy Logic Toolbox)
smoothing_func=sigmoid_function(spatial_correlations_centers,[20 min(spatial_correlations_centers)+0.5]);
spatial_correlations_model_same_cells=spatial_correlations_model_same_cells.*smoothing_func;
spatial_correlations_model_same_cells(1:round(number_of_bins/10:end))=0;
% calculating the distribution for different cells:
spatial_correlations_model_different_cells=betapdf(1-spatial_correlations_centers,spatial_correlations_model_parameters(3),spatial_correlations_model_parameters(5));
spatial_correlations_model_different_cells=spatial_correlations_model_different_cells./sum(spatial_correlations_model_different_cells)*(number_of_bins/((spatial_correlations_centers(2)-spatial_correlations_centers(1))+(spatial_correlations_centers(end)-spatial_correlations_centers(1))));
% calculating the weighted sum:
spatial_correlations_model_weighted_sum=spatial_correlations_model_parameters(1)*spatial_correlations_model_same_cells+(1-spatial_correlations_model_parameters(1))*spatial_correlations_model_different_cells;
[spatial_correlations_distribution,~]=hist(neighbors_spatial_correlations,spatial_correlations_centers);
spatial_correlations_distribution=spatial_correlations_distribution./sum(spatial_correlations_distribution)*(number_of_bins/((spatial_correlations_centers(2)-spatial_correlations_centers(1))+(spatial_correlations_centers(end)-spatial_correlations_centers(1))));

% calculating the discrepancy of the model (normalized mean sqaured error)
MSE_spatial_correlations_model=sum(abs(((spatial_correlations_distribution-spatial_correlations_model_weighted_sum))*((spatial_correlations_centers(2)-spatial_correlations_centers(1))+(spatial_correlations_centers(end)-spatial_correlations_centers(1)))/number_of_bins))/2;

% calculating the P_same of the model:
minimal_p_same_threshold=0.001;
p_same_given_spatial_correlation=spatial_correlations_model_parameters(1).*spatial_correlations_model_same_cells./(spatial_correlations_model_parameters(1).*spatial_correlations_model_same_cells+(1-spatial_correlations_model_parameters(1)).*spatial_correlations_model_different_cells);
indexes_to_smooth=find(spatial_correlations_model_same_cells<minimal_p_same_threshold*max(spatial_correlations_model_same_cells));
sigmoid_function=@(x,ac)1./(1+exp(-ac(1)*(x-ac(2)))); % defining a sigmoid function
smoothing_func=sigmoid_function(1:length(indexes_to_smooth),[0.05*length(indexes_to_smooth) 0.8*length(indexes_to_smooth)]);
p_same_given_spatial_correlation(indexes_to_smooth)=p_same_given_spatial_correlation(indexes_to_smooth).*smoothing_func;

% finding the intersection between same cells and different cells:
index_range_of_intersection=find(spatial_correlations_model_same_cells>minimal_p_same_threshold*max(spatial_correlations_model_same_cells));
index_range_of_intersection(end)=[];
[~,index_of_intersection]=min(abs(spatial_correlations_model_parameters(1)*spatial_correlations_model_same_cells(index_range_of_intersection)-(1-spatial_correlations_model_parameters(1))*spatial_correlations_model_different_cells(index_range_of_intersection)));
spatial_correlation_intersection=round(100*spatial_correlations_centers(index_of_intersection+index_range_of_intersection(1)-1))/100;

end

