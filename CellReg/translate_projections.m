function [translated_image]=translate_projections(original_image,translation)
% This function receives an image. It rotates the image according to the
% angle theta and translates the image in xy plane according to translated. The
% function returns the rotated and translated image.
N=size(original_image,1);
M=size(original_image,2);
a=translation(2);
b=translation(1);
translated_image=zeros(size(original_image));
for p=1:N
    for q=1:M
        wanted_coords=[p ;q]+[a; b];
        wanted_coords(wanted_coords<0)=0;
        if wanted_coords(1)>N+1;
            wanted_coords(1)=N+1;
        end
        if wanted_coords(2)>M+1;
            wanted_coords(1)=M+1;
        end
        translated_image(p,q)=coord_value(original_image,wanted_coords);
    end
end
end

