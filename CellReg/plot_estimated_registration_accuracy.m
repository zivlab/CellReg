function plot_estimated_registration_accuracy(p_same_centers_of_bins,p_same_certainty_threshold,p_same_given_centroid_distance,centroid_distances_distribution,cdf_p_same_centroid_distances,uncertain_fraction_centroid_distances,true_positive_per_distance_threshold,false_positive_per_distance_threshold,centers_of_bins,maximal_distance,microns_per_pixel,imaging_technique,figures_directory,figures_visibility,varargin)
% This function plots the estiamted fraction of uncertain cell pairs,
% and the false positive and false negative rates for the different models.

number_of_p_same_bins=length(p_same_centers_of_bins);
number_of_bins=length(centers_of_bins{1});
p_same_certainty_threshold=1-p_same_certainty_threshold;

figure('units','normalized','outerposition',[0.2 0.1 0.6 0.8],'Visible',figures_visibility)
set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')')
size_x=0.78;
size_y=0.78;
if strcmp(imaging_technique,'one_photon');
    p_same_given_spatial_correlation=varargin{1};
    spatial_correlations_distribution=varargin{2};
    cdf_p_same_spatial_correlations=varargin{3};
    uncertain_fraction_spatial_correlations=varargin{4};
    true_positive_per_correlation_threshold=varargin{5};
    false_positive_per_correlation_threshold=varargin{6};
        
    % spatial correlation model:
    axes('position',[0.08 0.58 size_x/2.3 size_y/2.85/2.3])
    start_x=0;
    end_x=1;
    y_vec=repmat(spatial_correlations_distribution,[2 1]);
    y_vec=y_vec(:);
    x_vec=(centers_of_bins{2}(2:end)+centers_of_bins{2}(1:end-1))/2;
    x_vec=repmat(x_vec,[2 1]);
    x_vec=[start_x; x_vec(:); end_x];
    for run_bins=1:length(x_vec)/2
        current_color=1-p_same_given_spatial_correlation(run_bins)*[1 1 1];
        patch(x_vec([1 1 2 2]+2*(run_bins-1)),[0 [1 1]*spatial_correlations_distribution(run_bins) 0],current_color,'EdgeColor',current_color)
        hold on
    end
    plot(x_vec,y_vec,'k-','linewidth',2)
    xlabel('Spatial correlation','fontsize',14,'fontweight','bold')
    x_label=linspace(0,1,6);
    x=linspace(0,1,6);
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',14,'fontweight','bold')
    xlim([0 1])
    [~,ind_005]=min(abs(0.05-(p_same_given_spatial_correlation)));
    p_005=centers_of_bins{2}(ind_005);
    [~,ind_05]=min(abs(0.5-(p_same_given_spatial_correlation)));
    p_05=centers_of_bins{2}(ind_05);
    [~,ind_095]=min(abs(0.95-(p_same_given_spatial_correlation)));
    p_095=centers_of_bins{2}(ind_095);
    hold on
    plot([p_005 p_005],[0 max(spatial_correlations_distribution)],'--','linewidth',3,'color','k')
    hold on
    plot([p_05 p_05],[0 max(spatial_correlations_distribution)],'--','linewidth',3,'color','k')
    hold on
    plot([p_095 p_095],[0 max(spatial_correlations_distribution)],'--','linewidth',3,'color','k')
    text(p_005,1.1*max(spatial_correlations_distribution),'0.05','fontsize',14,'fontweight','bold','HorizontalAlignment','Center')
    text(p_095,1.1*max(spatial_correlations_distribution),'0.95','fontsize',14,'fontweight','bold','HorizontalAlignment','Center')
    text(p_05,1.1*max(spatial_correlations_distribution),'0.5','fontsize',14,'fontweight','bold','HorizontalAlignment','Center')
end

% centroid distance model:
axes('position',[0.08 0.8 size_x/2.3 size_y/2.85/2.3])
start_x=0;
end_x=maximal_distance*microns_per_pixel;
y_vec=repmat(centroid_distances_distribution,[2 1]);
y_vec=y_vec(:);
x_vec=(microns_per_pixel*centers_of_bins{1}(2:end)+microns_per_pixel*centers_of_bins{1}(1:end-1))/2;
x_vec=repmat(x_vec,[2 1]);
x_vec=[start_x; x_vec(:); end_x];
for run_bins=1:length(x_vec)/2
    current_color=1-p_same_given_centroid_distance(run_bins)*[1 1 1];
    patch(x_vec([1 1 2 2]+2*(run_bins-1)),[0 [1 1]*centroid_distances_distribution(run_bins) 0],current_color,'EdgeColor',current_color)
    hold on
end
plot(x_vec,y_vec,'k-','linewidth',2)
xlim([0 microns_per_pixel*maximal_distance])
xlabel('Centroids distance (\mum)','FontWeight','Bold','fontsize',14)
x_label=0:3:microns_per_pixel*maximal_distance;
x=0:3:microns_per_pixel*maximal_distance;
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',14,'fontweight','bold')
[~,ind_005]=min(abs(0.05-(p_same_given_centroid_distance)));
p_005=microns_per_pixel*centers_of_bins{1}(ind_005);
[~,ind_05]=min(abs(0.5-(p_same_given_centroid_distance)));
p_05=microns_per_pixel*centers_of_bins{1}(ind_05);
[~,ind_095]=min(abs(0.95-(p_same_given_centroid_distance)));
p_095=microns_per_pixel*centers_of_bins{1}(ind_095);
hold on
plot([p_005 p_005],[0 max(centroid_distances_distribution)],'--','linewidth',3,'color','k')
hold on
plot([p_05 p_05],[0 max(centroid_distances_distribution)],'--','linewidth',3,'color','k')
hold on
plot([p_095 p_095],[0 max(centroid_distances_distribution)],'--','linewidth',3,'color','k')
text(p_005,1.1*max(centroid_distances_distribution),'0.05','fontsize',14,'fontweight','bold','HorizontalAlignment','Center')
text(p_095,1.1*max(centroid_distances_distribution),'0.95','fontsize',14,'fontweight','bold','HorizontalAlignment','Center')
text(p_05,1.1*max(centroid_distances_distribution),'0.5','fontsize',14,'fontweight','bold','HorizontalAlignment','Center')
axes('position',[0.44 0.58 0.03/2.3 size_y/2.3])
color_vec=linspace(1,0,number_of_bins);
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim([0 1])
ylim([0 1])
for n=1:number_of_bins
    hold on
    p=patch([0 1 1 0],[n/number_of_bins n/number_of_bins (n-1)/number_of_bins (n-1)/number_of_bins],[color_vec(n) color_vec(n) color_vec(n)]);
    set(p,'FaceAlpha',1,'EdgeColor','none');
end
text(2.5,0.5,'P_s_a_m_e','fontsize',14,'fontweight','bold','rotation',90,'HorizontalAlignment','Center')
text(1.5,0,'0','fontsize',14,'fontweight','bold','HorizontalAlignment','Left')
text(1.5,1,'1','fontsize',14,'fontweight','bold','HorizontalAlignment','Left')
text(-31.5,0.5,'Probability density','fontsize',14,'fontweight','bold','rotation',90,'HorizontalAlignment','Center')
plot([0 1 1 0 0],[0 0 1 1 0],'color','k')

% Cumulative distribution of P_same and uncertainty fraction:
axes('position',[0.08 0.08 size_x/2.3 size_y/2.3])
plot([0,cumsum(cdf_p_same_centroid_distances)/sum(cdf_p_same_centroid_distances),1],[0,p_same_centers_of_bins,1],'linewidth',2,'color','r')
hold on
if strcmp(imaging_technique,'one_photon');
    plot([0,cumsum(cdf_p_same_spatial_correlations)/sum(cdf_p_same_spatial_correlations),1],[0,p_same_centers_of_bins,1],'linewidth',2,'color','b')
    hold on
end
if strcmp(imaging_technique,'one_photon');
    plot([1-sum(cdf_p_same_spatial_correlations(1-p_same_centers_of_bins<p_same_certainty_threshold))/sum(cdf_p_same_spatial_correlations),1-sum(cdf_p_same_spatial_correlations(1-p_same_centers_of_bins<p_same_certainty_threshold))/sum(cdf_p_same_spatial_correlations)],[0,1],'--','linewidth',2,'color','k')
    hold on
    plot([sum(cdf_p_same_spatial_correlations(1-p_same_centers_of_bins>1-p_same_certainty_threshold))/sum(cdf_p_same_spatial_correlations),sum(cdf_p_same_spatial_correlations(1-p_same_centers_of_bins>1-p_same_certainty_threshold))/sum(cdf_p_same_spatial_correlations)],[0,1],'--','linewidth',2,'color','k')
end
xlabel('Fraction of cell pairs ','fontsize',14,'fontweight','bold')
ylabel('P_s_a_m_e','fontsize',14,'fontweight','bold')
hold on
plot([0 1],[p_same_certainty_threshold p_same_certainty_threshold],'--','linewidth',2,'color','k')
hold on
plot([0 1],[1-p_same_certainty_threshold 1-p_same_certainty_threshold],'--','linewidth',2,'color','k')
if strcmp(imaging_technique,'one_photon');
    legend('Dist.','Corr.','Location','Northeast')
end
legend('boxoff')
y=linspace(0,1,6);
y_label=linspace(0,1,6);
x=linspace(0,1,6);
x_label=linspace(0,1,6);
set(gca,'YTick',y)
set(gca,'YTickLabel',y_label,'fontsize',14,'fontweight','bold')
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',14,'fontweight','bold')
axes('position',[0.3 0.15 size_x/2.3/3 size_y/2.3/3.5])
labels=cell(1,3);
labels{1}='Dist.';
labels{2}='Corr.';
if strcmp(imaging_technique,'one_photon');
    bar([1 2],[uncertain_fraction_centroid_distances; uncertain_fraction_spatial_correlations],0.8,'FaceColor','none','EdgeColor','k')
    xtick_vec=1:2;
else
    bar(1,uncertain_fraction_centroid_distances,0.8,'FaceColor','none','EdgeColor','k')
    xtick_vec=1;
end
box off
ylabel('Uncertain fraction ','FontSize',10,'fontweight','bold')
set(gca,'XTick',xtick_vec)
set(gca,'XTickLabel',labels,'FontSize',10,'fontweight','bold')

% Estimated ROC curves:
axes('position',[0.59 0.08 size_x/2.3 size_y/2.3])
if strcmp(imaging_technique,'one_photon');
    plot(cumsum(false_positive_per_correlation_threshold),cumsum(true_positive_per_correlation_threshold),'linewidth',2,'color','b')
    hold on
end
plot([0,cumsum(false_positive_per_distance_threshold)],[0,cumsum(true_positive_per_distance_threshold)],'linewidth',2,'color','r')
hold on
plot([0 0.265],[0.9 0.14],'--','color',[0.8 0.8 0.8],'linewidth',2)
hold on
plot([0.1 0.77],[1 0.63],'--','color',[0.8 0.8 0.8],'linewidth',2)
ylabel('True positive rate','fontsize',14,'fontweight','bold')
xlabel('False positive rate','fontsize',14,'fontweight','bold')
ylim([0 1])
xlim([0 1])
p=patch([0 0.1 0.1 0],[1 1 0.9 0.9],[0.8 0.8 0.8]);
set(p,'FaceAlpha',0.3,'EdgeColor','none');
y=linspace(0,1,6);
y_label=linspace(0,1,6);
x=linspace(0,1,6);
x_label=linspace(0,1,6);
set(gca,'YTick',y)
set(gca,'YTickLabel',y_label,'fontsize',14,'fontweight','bold')
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',14,'fontweight','bold')
axes('position',[0.682 0.125 size_x/2.3/2 size_y/2.3/2])
plot([0,cumsum(false_positive_per_distance_threshold)],[0,cumsum(true_positive_per_distance_threshold)],'linewidth',2,'color','r')
if strcmp(imaging_technique,'one_photon');
    hold on
    plot(cumsum(false_positive_per_correlation_threshold),cumsum(true_positive_per_correlation_threshold),'linewidth',2,'color','b')
end
hold on
cum_sum_temp=[0,cumsum(false_positive_per_distance_threshold)];
cum_sum_temp_2=[0,cumsum(true_positive_per_distance_threshold)];
plot(cum_sum_temp(round(number_of_p_same_bins/2)),cum_sum_temp_2(round(number_of_p_same_bins/2)),'*','markersize',8,'linewidth',2,'color','k')
if strcmp(imaging_technique,'one_photon');
    hold on
    cum_sum_temp=[0,cumsum(false_positive_per_correlation_threshold)];
    cum_sum_temp_2=[0,cumsum(true_positive_per_correlation_threshold)];
    plot(cum_sum_temp(round(number_of_p_same_bins/2)),cum_sum_temp_2(round(number_of_p_same_bins/2)),'*','markersize',8,'linewidth',2,'color','k')
end

xlim([0 0.1])
ylim([0.9 1])
x_label=0:0.1:0.1;
x=0:0.1:0.2;
y=0.9:0.1:1;
y_label=0.9:0.1:1;
set(gca,'YTick',y)
set(gca,'YTickLabel',y_label,'fontsize',14,'fontweight','bold')
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',14,'fontweight','bold')
p=patch([0 1 1 0],[1 1 0 0],[0.8 0.8 0.8]);
set(p,'FaceAlpha',0.3,'EdgeColor','none');
set(gcf,'PaperPositionMode','auto')
savefig(fullfile(figures_directory,'Stage 3 - registration certainty.fig'))
saveas(gcf,fullfile(figures_directory,'Stage 3 - registration certainty'),'png')

end

