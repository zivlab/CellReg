function [rotated_image]=rotate_image_interp(original_image,theta,translations,center_of_FOV)
% This function rotates an image according to the
% angle theta, and translates the image in (x,y) according to translations. 

% Inputs:
% 1. original_image
% 2. theta - in degrees
% 3. translations
% 4. center_of_FOV

% Outputs:
% 1. rotated_image

theta=theta*pi/180;
N=size(original_image,1);
M=size(original_image,2);
a=translations(2);
b=translations(1);
transformation=[cos(theta) -sin(theta) ; sin(theta) cos(theta)]';
trans_inv=transformation^-1;
rotated_image=zeros(size(original_image));
for p=1:N
    for q=1:M
        wanted_coords=trans_inv*[p-center_of_FOV(1)-a ;q-center_of_FOV(2)-b]+[center_of_FOV(1) ;center_of_FOV(2)];
        wanted_coords(wanted_coords<0)=0;
        if wanted_coords(1)>N+1;
        wanted_coords(1)=N+1;
        end
        if wanted_coords(2)>M+1;
        wanted_coords(1)=M+1;
        end
        rotated_image(p,q)=interpolate_pixel_value(original_image,wanted_coords);
    end
end
end

