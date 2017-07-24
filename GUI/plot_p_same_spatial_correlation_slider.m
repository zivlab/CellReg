function plot_p_same_spatial_correlation_slider(spatial_correlations_distribution,p_same_given_spatial_correlation,spatial_correlation_threshold,centers_of_bins)
% This function plots the value of the p_same slider in the GUI
% on top of the used model.

% Inputs:
% 1. spatial_correlations_distribution
% 2. p_same_given_spatial_correlation
% 3. spatial_correlation_threshold
% 4. centers_of_bins


plot([],[])
set(gca,'xtick',[])
set(gca,'ytick',[])
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
xlabel('Spatial correlation','fontsize',12,'fontweight','bold')
x_label=linspace(0,1,6);
x=linspace(0,1,6);
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
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
hold on
plot([spatial_correlation_threshold spatial_correlation_threshold],[0 max(spatial_correlations_distribution)],'linewidth',2,'color','r')
text(p_005,1.1*max(spatial_correlations_distribution),'0.05','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
text(p_095,1.1*max(spatial_correlations_distribution),'0.95','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
text(p_05,1.1*max(spatial_correlations_distribution),'0.5','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')


end

