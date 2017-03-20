function [value]=coord_value(image,coords)
% This function receives an image and x',y' coordinates where the image
% should be estimated. The function returns the image value in the
% requested coordinates.
N=size(image,1);
M=size(image,2);
bilinear_mat=zeros(2,2);
if floor(coords(1))>=1 && floor(coords(1))<=N && floor(coords(2))>=1 && floor(coords(2))<=M
    bilinear_mat(1,1)=image(floor(coords(1)),floor(coords(2)));
else
    bilinear_mat(1,1)=0;
end
if ceil(coords(1))>=1 && ceil(coords(1))<=N && floor(coords(2))>=1 && floor(coords(2))<=M
    bilinear_mat(1,2)=image(ceil(coords(1)),floor(coords(2)));
else
    bilinear_mat(1,2)=0;
end
if floor(coords(1))>=1 && floor(coords(1))<=N && ceil(coords(2))>=1 && ceil(coords(2))<=M
    bilinear_mat(2,1)=image(floor(coords(1)),ceil(coords(2)));
else
    bilinear_mat(2,1)=0;
end
if ceil(coords(1))>=1 && ceil(coords(1))<=N && ceil(coords(2))>=1 && ceil(coords(2))<=M
    bilinear_mat(2,2)=image(ceil(coords(1)),ceil(coords(2)));
else
    bilinear_mat(2,2)=0;
end

alpha=coords(2)-floor(coords(2));
betha=coords(1)-floor(coords(1));
value=[1-alpha alpha]*bilinear_mat*[1-betha betha]';
end

