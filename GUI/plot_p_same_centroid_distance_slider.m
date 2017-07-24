function plot_p_same_centroid_distance_slider(centroid_distances_distribution,p_same_given_centroid_distance,centroid_distance_threshold,centers_of_bins,maximal_distance,microns_per_pixel)
% This function plots the value of the p_same slider in the GUI
% on top of the used model.

% Inputs:
% 1. centroid_distances_distribution
% 2. p_same_given_centroid_distance
% 3. centroid_distance_threshold
% 4. centers_of_bins
% 5. maximal_distance
% 6. microns_per_pixel

plot([],[])
set(gca,'xtick',[])
set(gca,'ytick',[])
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
xlabel('centroid distance (\mum)','FontWeight','Bold','fontsize',12)
x_label=0:3:microns_per_pixel*maximal_distance;
x=0:3:microns_per_pixel*maximal_distance;
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',12,'fontweight','bold')
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
hold on
plot([centroid_distance_threshold centroid_distance_threshold],[0 max(centroid_distances_distribution)],'linewidth',2,'color','r')
text(p_005,1.1*max(centroid_distances_distribution),'0.05','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
text(p_095,1.1*max(centroid_distances_distribution),'0.95','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')
text(p_05,1.1*max(centroid_distances_distribution),'0.5','fontsize',10,'fontweight','bold','HorizontalAlignment','Center')

end

