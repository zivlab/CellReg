function plot_models(centroid_distances_model_parameters,NN_centroid_distances,NNN_centroid_distances,centroid_distances_distribution,centroid_distances_model_same_cells,centroid_distances_model_different_cells,centroid_distances_model_weighted_sum,centroid_distance_intersection,centers_of_bins,microns_per_pixel,maximal_distance,figures_directory,figures_visibility,varargin)
% This function plots the results of the probabilistic models for the
% centroid distances and spatial correlations.

if ~isempty(varargin)
    spatial_correlations_model_parameters=varargin{1};
    NN_spatial_correlations=varargin{2};
    NNN_spatial_correlations=varargin{3};
    spatial_correlations_distribution=varargin{4};    
    spatial_correlations_model_same_cells=varargin{5};
    spatial_correlations_model_different_cells=varargin{6};
    spatial_correlations_model_weighted_sum=varargin{7};
    spatial_correlation_intersection=varargin{8};
end

number_of_bins=length(centers_of_bins{1});

figure('units','normalized','outerposition',[0.2 0.1 0.6 0.8],'Visible',figures_visibility)
set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')')
if ~isempty(varargin)
    subplot(2,2,2)
    [n1,~]=hist(NN_spatial_correlations,centers_of_bins{2});
    [n2,~]=hist(NNN_spatial_correlations,centers_of_bins{2});
    bar(centers_of_bins{2}+0.25/number_of_bins,n1,'g','EdgeColor','none','barwidth',0.5);
    hold on
    bar(centers_of_bins{2}-0.25/number_of_bins,n2,'r','EdgeColor','none','barwidth',0.5);
    xlim([0 1])
    legend('Nearest neighbors','Other neighbors','location','northwest')
    legend('boxoff')
    x_label=linspace(0,1,6);
    x=linspace(0,1,6);
    set(gca,'fontsize',14)
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',14)
    set(gcf,'PaperPositionMode','auto')
    xlabel('Spatial correlation','FontWeight','Bold','fontsize',14)
    ylabel('Number of cell-pairs','FontWeight','Bold','fontsize',14)
    set(gca,'fontsize',14)
end
subplot(2,2,1)
[n1,~]=hist(NN_centroid_distances,centers_of_bins{1});
[n2,~]=hist(NNN_centroid_distances,centers_of_bins{1});
bar(microns_per_pixel*centers_of_bins{1}+0.25*microns_per_pixel*maximal_distance/number_of_bins,n1,'g','EdgeColor','none','barwidth',0.5);
hold on
bar(microns_per_pixel*centers_of_bins{1}-0.25*microns_per_pixel*maximal_distance/number_of_bins,n2,'r','EdgeColor','none','barwidth',0.5);
xlim([0 microns_per_pixel*maximal_distance])
x_label=0:3:microns_per_pixel*maximal_distance;
x=0:3:microns_per_pixel*maximal_distance;
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',14)
xlabel('Centroids distance (\mum)','FontWeight','Bold','fontsize',14)
ylabel('Number of cell-pairs','FontWeight','Bold','fontsize',14)
set(gca,'fontsize',14)
if ~isempty(varargin)
    subplot(2,2,4)
    bar(centers_of_bins{2},spatial_correlations_distribution,'FaceColor','b','EdgeColor','none','barwidth',1);
    hold on
    plot(centers_of_bins{2},spatial_correlations_model_parameters(1)*spatial_correlations_model_same_cells,'--','linewidth',3,'color','g');
    hold on
    plot(centers_of_bins{2},(1-spatial_correlations_model_parameters(1))*spatial_correlations_model_different_cells,'--','linewidth',3,'color','r');
    hold on
    plot(centers_of_bins{2},spatial_correlations_model_weighted_sum,'linewidth',3,'color','k');
    hold on
    plot(centers_of_bins{2},spatial_correlations_model_parameters(1)*spatial_correlations_model_same_cells,'--','linewidth',3,'color','g');
    hold on
    plot(centers_of_bins{2},(1-spatial_correlations_model_parameters(1))*spatial_correlations_model_different_cells,'--','linewidth',3,'color','r');
    hold on
    plot([spatial_correlation_intersection spatial_correlation_intersection],[0 max(spatial_correlations_distribution)],'--','linewidth',2,'color','k')
    xlim([0 1])
    xlabel('Spatial correlation','fontsize',14,'fontweight','bold')
    ylabel('Probability density','fontsize',14,'fontweight','bold')
    hold on
    x_label=linspace(0,1,6);
    x=linspace(0,1,6);
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',14)
    legend('Observed data','Same cell model','Different cells model','Overall model','location','northwest')
    legend('boxoff')
    xlim([0 1])
    set(gca,'fontsize',14)
    normalized_same=spatial_correlations_model_same_cells./sum(spatial_correlations_model_same_cells);
    normalized_diff=spatial_correlations_model_different_cells./sum(spatial_correlations_model_different_cells);
    same_more_than_thresh=sum(normalized_same(centers_of_bins{2}>spatial_correlation_intersection));
    diff_more_than_thresh=sum(normalized_diff(centers_of_bins{2}>spatial_correlation_intersection));
    text(spatial_correlation_intersection+0.1,0.9*max(spatial_correlations_distribution),[num2str(round(100*same_more_than_thresh)) '%'],'fontsize',14,'fontweight','bold','HorizontalAlignment','Center','color','g')
    text(spatial_correlation_intersection-0.1,0.9*max(spatial_correlations_distribution),[num2str(round(100*(1-same_more_than_thresh))) '%'],'fontsize',14,'fontweight','bold','HorizontalAlignment','Center','color','g')
    text(spatial_correlation_intersection+0.1,0.8*max(spatial_correlations_distribution),[num2str(round(100*(diff_more_than_thresh))) '%'],'fontsize',14,'fontweight','bold','HorizontalAlignment','Center','color','r')
    text(spatial_correlation_intersection-0.1,0.8*max(spatial_correlations_distribution),[num2str(round(100*(1-diff_more_than_thresh))) '%'],'fontsize',14,'fontweight','bold','HorizontalAlignment','Center','color','r')
end
subplot(2,2,3)
bar(microns_per_pixel*centers_of_bins{1},centroid_distances_distribution,'FaceColor','b','EdgeColor','none','barwidth',1);
hold on
plot(microns_per_pixel*centers_of_bins{1},centroid_distances_model_parameters(1)*centroid_distances_model_same_cells,'--','color','g','linewidth',3)
hold on
plot(microns_per_pixel*centers_of_bins{1},(1-centroid_distances_model_parameters(1))*centroid_distances_model_different_cells,'--','color','r','linewidth',3)
hold on
plot(microns_per_pixel*centers_of_bins{1},centroid_distances_model_weighted_sum,'color','r','linewidth',3,'color','k')
hold on
plot(microns_per_pixel*centers_of_bins{1},centroid_distances_model_parameters(1)*centroid_distances_model_same_cells,'--','color','g','linewidth',3)
hold on
plot(microns_per_pixel*centers_of_bins{1},(1-centroid_distances_model_parameters(1))*centroid_distances_model_different_cells,'--','color','r','linewidth',3)
hold on
plot([centroid_distance_intersection centroid_distance_intersection],[0 max(centroid_distances_distribution)],'--','linewidth',2,'color','k')
xlim([0 microns_per_pixel*maximal_distance])
xlabel('Centroids distance (\mum)','FontWeight','Bold','fontsize',14)
ylabel('Probability density','FontWeight','Bold','fontsize',14)
x_label=0:3:microns_per_pixel*maximal_distance;
x=0:3:microns_per_pixel*maximal_distance;
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',14)
set(gca,'fontsize',14)
normalized_same=centroid_distances_model_same_cells./sum(centroid_distances_model_same_cells);
normalized_diff=centroid_distances_model_different_cells./sum(centroid_distances_model_different_cells);
same_more_than_thresh=sum(normalized_same(centers_of_bins{1}*microns_per_pixel>centroid_distance_intersection));
diff_more_than_thresh=sum(normalized_diff(centers_of_bins{1}*microns_per_pixel>centroid_distance_intersection));
text(centroid_distance_intersection+1,0.9*max(centroid_distances_distribution),[num2str(round(100*same_more_than_thresh)) '%'],'fontsize',14,'fontweight','bold','HorizontalAlignment','Center','color','g')
text(centroid_distance_intersection-1,0.9*max(centroid_distances_distribution),[num2str(round(100*(1-same_more_than_thresh))) '%'],'fontsize',14,'fontweight','bold','HorizontalAlignment','Center','color','g')
text(centroid_distance_intersection+1,0.8*max(centroid_distances_distribution),[num2str(round(100*(diff_more_than_thresh))) '%'],'fontsize',14,'fontweight','bold','HorizontalAlignment','Center','color','r')
text(centroid_distance_intersection-1,0.8*max(centroid_distances_distribution),[num2str(round(100*(1-diff_more_than_thresh))) '%'],'fontsize',14,'fontweight','bold','HorizontalAlignment','Center','color','r')
set(gcf,'PaperPositionMode','auto')
savefig(fullfile(figures_directory,'Stage 3 - model.fig'))
saveas(gcf,fullfile(figures_directory,'Stage 3 - model'),'png')

end

