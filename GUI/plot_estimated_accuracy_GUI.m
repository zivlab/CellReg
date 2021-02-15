function plot_estimated_accuracy_GUI(handles,p_val_vec,p_thresh,p_value_dist,normalized_distance,num_pairs_dist,true_merge_dist,false_merge_dist,centers_of_bins,maximal_distance,microns_per_pixel,varargin)
% This function plots in the GUI the estiamted fraction of uncertain cell pairs,
% and the false positive and false negative rates for the different models.

axes(handles.axes3)
start_x=0;
end_x=maximal_distance*microns_per_pixel;
y_vec=repmat(normalized_distance,[2 1]);
y_vec=y_vec(:);
x_vec=(microns_per_pixel*centers_of_bins{1}(2:end)+microns_per_pixel*centers_of_bins{1}(1:end-1))/2;
x_vec=repmat(x_vec,[2 1]);
x_vec=[start_x; x_vec(:); end_x];
for run_bins=1:length(x_vec)/2
    current_color=1-p_value_dist(run_bins)*[1 1 1];
    patch(x_vec([1 1 2 2]+2*(run_bins-1)),[0 [1 1]*normalized_distance(run_bins) 0],current_color,'EdgeColor',current_color)
    hold on
end
plot(x_vec,y_vec,'k-','linewidth',2)
xlim([0 microns_per_pixel*maximal_distance])
xlabel('Centroids distance (\mum)','FontWeight','Bold','fontsize',12)
x_label=0:3:microns_per_pixel*maximal_distance;
x=0:3:microns_per_pixel*maximal_distance;
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
[~,ind_005]=min(abs(0.05-(p_value_dist)));
p_005=microns_per_pixel*centers_of_bins{1}(ind_005);
[~,ind_05]=min(abs(0.5-(p_value_dist)));
p_05=microns_per_pixel*centers_of_bins{1}(ind_05);
[~,ind_095]=min(abs(0.95-(p_value_dist)));
p_095=microns_per_pixel*centers_of_bins{1}(ind_095);
hold on
plot([p_005 p_005],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
hold on
plot([p_05 p_05],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
hold on
plot([p_095 p_095],[0 max(normalized_distance)],'--','linewidth',3,'color','k')
text(p_005,1.1*max(normalized_distance),'0.05','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
text(p_095,1.1*max(normalized_distance),'0.95','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
text(p_05,1.1*max(normalized_distance),'0.5','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')

p_value_corr=varargin{1};
n_corr=varargin{2};
num_pairs_corr=varargin{3};
true_merge_corr=varargin{4};
false_merge_corr=varargin{5};

axes(handles.axes4)
start_x=0;
end_x=1;
y_vec=repmat(n_corr,[2 1]);
y_vec=y_vec(:);
x_vec=(centers_of_bins{2}(2:end)+centers_of_bins{2}(1:end-1))/2;
x_vec=repmat(x_vec,[2 1]);
x_vec=[start_x; x_vec(:); end_x];
for run_bins=1:length(x_vec)/2
    current_color=1-p_value_corr(run_bins)*[1 1 1];
    patch(x_vec([1 1 2 2]+2*(run_bins-1)),[0 [1 1]*n_corr(run_bins) 0],current_color,'EdgeColor',current_color)
    hold on
end
plot(x_vec,y_vec,'k-','linewidth',2)
xlabel('Spatial correlation','fontsize',12,'fontweight','bold')
x_label=linspace(0,1,6);
x=linspace(0,1,6);
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
xlim([0 1])
[~,ind_005]=min(abs(0.05-(p_value_corr)));
p_005=centers_of_bins{2}(ind_005);
[~,ind_05]=min(abs(0.5-(p_value_corr)));
p_05=centers_of_bins{2}(ind_05);
[~,ind_095]=min(abs(0.95-(p_value_corr)));
p_095=centers_of_bins{2}(ind_095);
hold on
plot([p_005 p_005],[0 max(n_corr)],'--','linewidth',3,'color','k')
hold on
plot([p_05 p_05],[0 max(n_corr)],'--','linewidth',3,'color','k')
hold on
plot([p_095 p_095],[0 max(n_corr)],'--','linewidth',3,'color','k')
text(p_005,1.1*max(n_corr),'0.05','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
text(p_095,1.1*max(n_corr),'0.95','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
text(p_05,1.1*max(n_corr),'0.5','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')

axes(handles.axes5)
cla;
plot([0,cumsum(num_pairs_dist)/sum(num_pairs_dist),1],[0,p_val_vec,1],'linewidth',2,'color','r')
hold on
plot([0,cumsum(num_pairs_corr)/sum(num_pairs_corr),1],[0,p_val_vec,1],'linewidth',2,'color','b')
hold on
plot([1-sum(num_pairs_corr(1-p_val_vec<p_thresh))/sum(num_pairs_corr),1-sum(num_pairs_corr(1-p_val_vec<p_thresh))/sum(num_pairs_corr)],[0,1],'--','linewidth',2,'color','k')
hold on
plot([sum(num_pairs_corr(1-p_val_vec>1-p_thresh))/sum(num_pairs_corr),sum(num_pairs_corr(1-p_val_vec>1-p_thresh))/sum(num_pairs_corr)],[0,1],'--','linewidth',2,'color','k')
xlabel('Fraction of cell pairs ','fontsize',14,'fontweight','bold')
ylabel('P_s_a_m_e','fontsize',14,'fontweight','bold')
hold on
plot([0 1],[p_thresh p_thresh],'--','linewidth',2,'color','k')
hold on
plot([0 1],[1-p_thresh 1-p_thresh],'--','linewidth',2,'color','k')
y=linspace(0,1,6);
y_label=linspace(0,1,6);
x=linspace(0,1,6);
x_label=linspace(0,1,6);
set(gca,'YTick',y)
set(gca,'YTickLabel',y_label,'fontsize',14)
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',14)

axes(handles.axes6)
plot([0,cumsum(false_merge_dist)],[0,cumsum(true_merge_dist)],'linewidth',2,'color','r')
hold on
plot(cumsum(false_merge_corr),cumsum(true_merge_corr),'linewidth',2,'color','b')
ylabel('True positive rate','fontsize',14,'fontweight','bold')
xlabel('False positive rate','fontsize',14,'fontweight','bold')
ylim([0 1])
xlim([0 1])
set(gca,'fontsize',14)
set(gcf,'PaperPositionMode','auto')
legend('Centroid distance','Spatial correlation','Location',[0.89 0.3 0.01 0.01])
legend('boxoff')

end

