function [rotated_spatial_footprint]=rotate_cell(original_spatial_footprint,theta,translations,center_of_FOV,centroid_location,microns_per_pixel)
% This function rotates a spatial footprint according to the
% angle theta, and translates the footprint in (x,y) according to translations. 

% Inputs:
% 1. original_spatial_footprint
% 2. theta - angle in degrees
% 3. translations
% 4. center_of_FOV - axis of rotation
% 5. centroid_location - is used to rotate only values in a certain radius
% 6. microns_per_pixel

% Outputs:
% 1. rotated_spatial_footprint

theta=theta*pi/180;
N=size(original_spatial_footprint,1);
M=size(original_spatial_footprint,2);
a=translations(2);
b=translations(1);
transformation=[cos(theta) -sin(theta) ; sin(theta) cos(theta)]';
trans_inv=transformation^-1;
rotated_spatial_footprint=zeros(size(original_spatial_footprint));
for p=1:N
    for q=1:M
        wanted_coords=trans_inv*[p-center_of_FOV(1)-a ;q-center_of_FOV(2)-b]+[center_of_FOV(1) ;center_of_FOV(2)];
        if sqrt(sum(((wanted_coords-centroid_location').^2)))>25/microns_per_pixel
            rotated_spatial_footprint(p,q)=0;
        else
            wanted_coords(wanted_coords<0)=0;
            if wanted_coords(1)>N+1;
                wanted_coords(1)=N+1;
            end
            if wanted_coords(2)>M+1;
                wanted_coords(1)=M+1;
            end
            rotated_spatial_footprint(p,q)=interpolate_pixel_value(original_spatial_footprint,wanted_coords);
        end
    end
end
end

