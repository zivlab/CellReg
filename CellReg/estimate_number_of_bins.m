function [number_of_bins,centers_of_bins]=estimate_number_of_bins(spatial_footprints,maximal_distance)
% This function estiamtes the number of bins to use for the probability
% distributions that are modeled. and then creates the centerss with the
% centers of the bins

% Inputs:
% 1. spatial_footprints

% Outputs:
% 1. number_of_bins
% 2. centers_of_bins

minimal_number_of_bins=40;
maximal_number_of_bins=100;
number_of_cells_factor=600;

number_of_sessions=size(spatial_footprints,2);
number_of_bins=10*round((minimal_number_of_bins+length(spatial_footprints{1})/number_of_cells_factor*(number_of_sessions)^2)/10);
if number_of_bins>maximal_number_of_bins
    number_of_bins=maximal_number_of_bins;
end

centers_of_bins=cell(1,2);
distances_centers_temp=linspace(0,maximal_distance,2*number_of_bins+1);
distances_centers=distances_centers_temp(2:2:end);
centers_of_bins{1}=distances_centers;
corrlations_centers_temp=linspace(0,1,2*number_of_bins+1);
corrlations_centers=corrlations_centers_temp(2:2:end);
centers_of_bins{2}=corrlations_centers;

end

