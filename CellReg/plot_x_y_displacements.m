function [x_y_displacements]=plot_x_y_displacements(neighbors_x_displacements,neighbors_y_displacements,microns_per_pixel,maximal_distance,number_of_bins,centers_of_bins,figures_directory,figures_visibility)
% This function plots the (x,y) distribution of cell-pair displacements

% Inputs:
% 1. neighbors_x_displacements
% 2 .neighbors_y_displacements
% 3. microns_per_pixel
% 4. maximal_distance
% 5. number_of_bins
% 6. centers_of_bins
% 7. figures_directory
% 8. figures_visibility

centers_of_bins_xy=cell(1,2);
xout_temp_2=linspace(0,maximal_distance,number_of_bins+1);
xout_2=xout_temp_2(2:2:end);
centers_of_bins_xy{1}=[-flip(xout_2), xout_2];
centers_of_bins_xy{2}=[-flip(xout_2), xout_2];

x_y_displacements=hist3([neighbors_x_displacements;neighbors_y_displacements]',centers_of_bins_xy);
x_y_displacements=flipud(x_y_displacements);
x_y_displacements=fliplr(x_y_displacements);

figure('units','normalized','outerposition',[0.35 0.25 0.3 0.5],'Visible',figures_visibility)
set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')')
set(gcf,'PaperOrientation','portrait');
size_x=0.75;
size_y=0.75;
axes('position',[0.12 0.15 size_x size_y])
imagesc(log(1+x_y_displacements)./max(max(log(1+x_y_displacements))))
axis square
cmap_jet=colormap('jet');
y=round(linspace(1,number_of_bins,9));
y_label=round(linspace(microns_per_pixel*max(centers_of_bins{1}),-microns_per_pixel*max(centers_of_bins{1}),9));
x=round(linspace(1,number_of_bins,9));
x_label=round(linspace(-microns_per_pixel*max(centers_of_bins{1}),microns_per_pixel*max(centers_of_bins{1}),9));set(gca,'YTick',y)
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
axes('position',[0.855 0.15 0.02 size_y])
num_colors=size(cmap_jet,1);
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim([0 1])
ylim([0 1])
for n=1:num_colors
    hold on
    p=patch([0 1 1 0],[n/num_colors n/num_colors (n-1)/num_colors (n-1)/num_colors],cmap_jet(n,:));
    set(p,'FaceAlpha',1,'EdgeColor','none');
end
text(3.5,0.5,'Number of cell-pairs (log)','fontsize',14,'fontweight','bold','rotation',90,'HorizontalAlignment','Center')
text(1.5,0,'0','fontsize',14,'fontweight','bold','HorizontalAlignment','Left')
text(1.5,1,'Max','fontsize',14,'fontweight','bold','HorizontalAlignment','Left')
set(gca,'fontsize',14)
set(gcf,'PaperPositionMode','auto')
savefig(fullfile(figures_directory,'Stage 3 - (x,y) displacements.fig'))
saveas(gcf,fullfile(figures_directory,'Stage 3 - (x,y) displacements'),'png')

end

