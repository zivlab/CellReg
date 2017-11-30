function plot_initial_registration(cell_to_index_map,number_of_bins,spatial_footprints,initial_registration_type,figures_directory,figures_visibility,varargin)
% This function plots the results for the initial registration stage

% Inputs:
% 1. cell_to_index_map
% 2. number_of_bins
% 3. spatial_footprints
% 4. initial_registration_type
% 5. figures_directory
% 6. figures_visibility
% 7. varargin
%   7{1}. registered_cells
%   7{2}. non_registered_cells
%   7{3}. pixel_to_mic
%   7{4}. maximal_distance

registered_cells=varargin{1};
non_registered_cells=varargin{2};

if strcmp(initial_registration_type,'Spatial correlation') % if spatial correlations are used 
    figure('units','normalized','outerposition',[0.3 0.25 0.4 0.5],'Visible',figures_visibility)
    set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')')
    xout=linspace(0,1,number_of_bins);
    [n1,~]=hist(registered_cells,xout);
    [n2,~]=hist(non_registered_cells,xout);
    bar(xout+0.25/number_of_bins,n1,'g','EdgeColor','none','barwidth',0.5);
    hold on
    bar(xout-0.25/number_of_bins,n2,'r','EdgeColor','none','barwidth',0.5);
    xlim([0 1])
    x_label=linspace(0,1,6);
    x=linspace(0,1,6);
    set(gca,'fontsize',14)
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',14)
    xlabel('Spatial correlation','FontWeight','Bold','fontsize',14)
    ylabel('Number of cell-pairs','FontWeight','Bold','fontsize',14)
    set(gca,'fontsize',14)
    legend('Same Cell','Different Cells','location','northwest')
    legend('boxoff')
else
    pixel_to_mic=varargin{3};
    maximal_distance=varargin{4};    
    figure('units','normalized','outerposition',[0.3 0.25 0.4 0.5],'Visible',figures_visibility)
    set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')')
    xout=linspace(0,maximal_distance,number_of_bins);
    [n1,~]=hist(registered_cells,xout);
    [n2,~]=hist(non_registered_cells,xout);
    bar(pixel_to_mic*xout+0.25*pixel_to_mic*maximal_distance/number_of_bins,n1,'g','EdgeColor','none','barwidth',0.5);
    hold on
    bar(pixel_to_mic*xout-0.25*pixel_to_mic*maximal_distance/number_of_bins,n2,'r','EdgeColor','none','barwidth',0.5);
    xlim([0 pixel_to_mic*maximal_distance])
    x_label=0:3:pixel_to_mic*maximal_distance;
    x=0:3:pixel_to_mic*maximal_distance;
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',14)
    xlabel('Centroids distance (\mum)','FontWeight','Bold','fontsize',14)
    ylabel('Number of cell-pairs','FontWeight','Bold','fontsize',14)
    set(gca,'fontsize',14)
    legend('Same Cell','Different Cells','location','northwest')
    legend('boxoff')
end
set(gcf,'PaperPositionMode','auto')
savefig(fullfile(figures_directory,'Stage 4 - same versus different cells.fig'))
saveas(gcf,fullfile(figures_directory,'Stage 4 - same versus different cells'),'png')

% Plotting the registration results with the cell maps from all sessions:
is_initial_stage=true;
plot_all_registered_projections(spatial_footprints,cell_to_index_map,figures_directory,figures_visibility,is_initial_stage)

end

