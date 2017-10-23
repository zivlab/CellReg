function plot_x_y_displacements_GUI(x_y_displacements,microns_per_pixel,centers_of_bins,maximal_distance,number_of_bins)
% This function plots the (x,y) distribution of cell-pair displacements in
% the GUI

% Inputs:
% 1. x_y_displacements
% 2. microns_per_pixel
% 3. centers_of_bins
% 4. maximal_distance
% 5. number_of_bins

imagesc(log(1+x_y_displacements)./max(max(log(1+x_y_displacements))))
title('Centroid distances','fontsize',14)
colormap('jet');
freezeColors
y=round(linspace(1,number_of_bins,9));
y_label=round(linspace(microns_per_pixel*max(centers_of_bins{1}),-microns_per_pixel*max(centers_of_bins{1}),9));
x=round(linspace(1,number_of_bins,9));
x_label=round(linspace(-microns_per_pixel*max(centers_of_bins{1}),microns_per_pixel*max(centers_of_bins{1}),9));
set(gca,'YTick',y)
set(gca,'YTickLabel',y_label,'fontsize',14)
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',14)
xlabel('x displacement (\mum)','FontWeight','Bold','fontsize',14)
ylabel('y displacement (\mum)','FontWeight','Bold','fontsize',14)
x_for_circle=0:0.01:2*pi;
hold on
plot(number_of_bins/2+number_of_bins/2*4/maximal_distance/microns_per_pixel*sin(x_for_circle),number_of_bins/2+number_of_bins/2*4/maximal_distance/microns_per_pixel*cos(x_for_circle),':','color',[1 1 1],'linewidth',4);
hold on
plot(number_of_bins/2+number_of_bins/2*8/maximal_distance/microns_per_pixel*sin(x_for_circle),number_of_bins/2+number_of_bins/2*8/maximal_distance/microns_per_pixel*cos(x_for_circle),'--','color',[1 1 1],'linewidth',4);

end

