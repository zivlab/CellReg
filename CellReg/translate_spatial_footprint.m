function [translated_spatial_footprint]=translate_spatial_footprint(original_spatial_footprint,translations,centroid_location,microns_per_pixel)
% This function translates a spatial footprint according to the provided
% translations.

% Inputs:
% 1. original_spatial_footprint
% 2. translations
% 3. centroid_location - is used to translate only values in a certain radius
% 4. microns_per_pixel 

% Outputs:
% 1.translated_spatial_footprint

maximal_cell_radius=25; % in microns

N=size(original_spatial_footprint,1);
M=size(original_spatial_footprint,2);
a=translations(2);
b=translations(1);
translated_spatial_footprint=zeros(size(original_spatial_footprint));
for p=1:N
    for q=1:M
        wanted_coords=[p ;q]+[a; b];
        if sqrt(sum(((wanted_coords-centroid_location').^2)))>maximal_cell_radius/microns_per_pixel
            translated_spatial_footprint(p,q)=0;
        else
            wanted_coords(wanted_coords<0)=0;
            if wanted_coords(1)>N+1;
                wanted_coords(1)=N+1;
            end
            if wanted_coords(2)>M+1;
                wanted_coords(1)=M+1;
            end
            translated_spatial_footprint(p,q)=interpolate_pixel_value(original_spatial_footprint,wanted_coords);
        end
    end
end
end

