function [interpolated_pixel_value]=interpolate_pixel_value(image,coordinates)
% This function computes the interpolated pixel value of an image in 
% requested coordinates.

% Inputs:
% 1. image
% 2. coordinates

% Outputs:
% 1. interpolated_pixel_value

N=size(image,1);
M=size(image,2);
bilinear_mat=zeros(2,2);
if floor(coordinates(1))>=1 && floor(coordinates(1))<=N && floor(coordinates(2))>=1 && floor(coordinates(2))<=M
    bilinear_mat(1,1)=image(floor(coordinates(1)),floor(coordinates(2)));
else
    bilinear_mat(1,1)=0;
end
if ceil(coordinates(1))>=1 && ceil(coordinates(1))<=N && floor(coordinates(2))>=1 && floor(coordinates(2))<=M
    bilinear_mat(1,2)=image(ceil(coordinates(1)),floor(coordinates(2)));
else
    bilinear_mat(1,2)=0;
end
if floor(coordinates(1))>=1 && floor(coordinates(1))<=N && ceil(coordinates(2))>=1 && ceil(coordinates(2))<=M
    bilinear_mat(2,1)=image(floor(coordinates(1)),ceil(coordinates(2)));
else
    bilinear_mat(2,1)=0;
end
if ceil(coordinates(1))>=1 && ceil(coordinates(1))<=N && ceil(coordinates(2))>=1 && ceil(coordinates(2))<=M
    bilinear_mat(2,2)=image(ceil(coordinates(1)),ceil(coordinates(2)));
else
    bilinear_mat(2,2)=0;
end

alpha=coordinates(2)-floor(coordinates(2));
betha=coordinates(1)-floor(coordinates(1));
interpolated_pixel_value=[1-alpha alpha]*bilinear_mat*[1-betha betha]';
end

