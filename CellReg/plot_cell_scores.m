function plot_cell_scores(cell_scores_positive,cell_scores_negative,cell_scores_exclusive,cell_scores,p_same_registered_pairs,figures_directory,figures_visibility)
% This function plots the distributions of false postive, false negative,
% exclusivity, and cell scores for all registered cells according to the
% clustering procedure

% Inputs:
% 1. cell_scores_positive
% 2. cell_scores_negative
% 3. cell_scores_exclusive
% 4. cell_scores
% 5. p_same_registered_pairs
% 6. figures_directory
% 7. figures_visibility

number_of_clusters=size(cell_scores,2);

number_of_bins=41;
xout_temp=linspace(0,1,number_of_bins);
xout=xout_temp(2:2:end);
figure('units','normalized','outerposition',[0.25 0.2 0.5 0.6],'Visible',figures_visibility)
set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')')
set(gcf,'PaperOrientation','portrait');
size_x=0.65;
size_y=0.65;

% cell scores:
axes('position',[0.6 0.1 size_x/2 size_y/2])
[n1,~]=hist(cell_scores,xout);
n1=n1./sum(n1);
bar(xout,n1,1)
xlim([0 1])
ylim([0 1])
xlabel('Overall cell scores','fontsize',14,'fontweight','bold')
ylabel('Probability','fontsize',14,'fontweight','bold')
x_label=linspace(0,1,6);
x=linspace(0,1,6);
set(gca,'fontsize',14)
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',14,'fontweight','bold')
h=axes('position',[0.68 0.25 size_x/6 size_y/6]);
plot(flip(xout),cumsum(flip(n1)),'linewidth',2)
ylim([0 1])
x_label=linspace(0,1,3);
x=linspace(0,1,3);
y=linspace(0,1,3);
y_label=linspace(0,1,3);
set(gca,'YTick',y)
set(gca,'YTickLabel',y_label,'fontsize',14,'fontweight','bold')
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',14,'fontweight','bold')
set(h, 'Xdir', 'reverse')
xlabel('Score','fontsize',14,'fontweight','bold')
ylabel('Cum. fraction','fontsize',14,'fontweight','bold')

% exclusivity scores:
axes('position',[0.12 0.1 size_x/2 size_y/2])
[n1,~]=hist(cell_scores_exclusive,xout);
n1=n1./sum(n1);
bar(xout,n1,1)
xlim([0 1])
ylim([0 1])
xlabel('Exclusivity cell scores','fontsize',14,'fontweight','bold')
ylabel('Probability','fontsize',14,'fontweight','bold')
x_label=linspace(0,1,6);
x=linspace(0,1,6);
set(gca,'fontsize',14)
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',14,'fontweight','bold')
h=axes('position',[0.2 0.25 size_x/6 size_y/6]);
plot(flip(xout),cumsum(flip(n1)),'linewidth',2)
ylim([0 1])
x_label=linspace(0,1,3);
x=linspace(0,1,3);
y=linspace(0,1,3);
y_label=linspace(0,1,3);
set(gca,'YTick',y)
set(gca,'YTickLabel',y_label,'fontsize',14,'fontweight','bold')
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',14,'fontweight','bold')
set(h, 'Xdir', 'reverse')
xlabel('Score','fontsize',14,'fontweight','bold')
ylabel('Cum. fraction','fontsize',14,'fontweight','bold')

% true positive scores:
axes('position',[0.6 0.58 size_x/2 size_y/2])
[n1,~]=hist(cell_scores_positive,xout);
n1=n1./sum(n1);
bar(xout,n1,1)
xlim([0 1])
ylim([0 1])
xlabel('True positive scores','fontsize',14,'fontweight','bold')
ylabel('Probability','fontsize',14,'fontweight','bold')
x_label=linspace(0,1,6);
x=linspace(0,1,6);
set(gca,'fontsize',14)
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',14,'fontweight','bold')
text(-0.25,1.2,[num2str(number_of_clusters) ' registered cells'],'fontsize',16,'fontweight','bold','HorizontalAlignment','Center')
h=axes('position',[0.68 0.73 size_x/6 size_y/6]);
plot(flip(xout),cumsum(flip(n1)),'linewidth',2)
ylim([0 1])
x_label=linspace(0,1,3);
x=linspace(0,1,3);
y=linspace(0,1,3);
y_label=linspace(0,1,3);
set(gca,'YTick',y)
set(gca,'YTickLabel',y_label,'fontsize',14,'fontweight','bold')
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',14,'fontweight','bold')
set(h, 'Xdir', 'reverse')
xlabel('Score','fontsize',14,'fontweight','bold')
ylabel('Cum. fraction','fontsize',14,'fontweight','bold')

% true negative scores:
axes('position',[0.12 0.58 size_x/2 size_y/2])
[n1,~]=hist(cell_scores_negative,xout);
n1=n1./sum(n1);
bar(xout,n1,1)
xlim([0 1])
ylim([0 1])
xlabel('True negative scores','fontsize',14,'fontweight','bold')
ylabel('Probability','fontsize',14,'fontweight','bold')
x_label=linspace(0,1,6);
x=linspace(0,1,6);
set(gca,'fontsize',14)
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',14,'fontweight','bold')
h=axes('position',[0.2 0.73 size_x/6 size_y/6]);
plot(flip(xout),cumsum(flip(n1)),'linewidth',2)
ylim([0 1])
x_label=linspace(0,1,3);
x=linspace(0,1,3);
y=linspace(0,1,3);
y_label=linspace(0,1,3);
set(gca,'YTick',y)
set(gca,'YTickLabel',y_label,'fontsize',14,'fontweight','bold')
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',14,'fontweight','bold')
set(h, 'Xdir', 'reverse')
xlabel('Score','fontsize',14,'fontweight','bold')
ylabel('Cum. fraction','fontsize',14,'fontweight','bold')
set(gcf,'PaperPositionMode','auto')
savefig(fullfile(figures_directory,'Stage 5 - cell scores.fig'))
saveas(gcf,fullfile(figures_directory,'Stage 5 - cell scores'),'png')

% plotting the distribution of P_same for all registered cell-pairs
number_of_sessions=size(p_same_registered_pairs{1},1);
all_pairs_p_same=zeros(1,number_of_clusters*number_of_sessions*number_of_sessions);
pair_count=0;
for n=1:number_of_clusters
    for k=1:number_of_sessions
        for m=k+1:number_of_sessions
            if ~isnan(p_same_registered_pairs{n}(k,m))
                pair_count=pair_count+1;
                all_pairs_p_same(pair_count)=p_same_registered_pairs{n}(k,m);
            end
        end
    end
end
all_pairs_p_same(pair_count+1:end)=[];

figure('units','normalized','outerposition',[0.35 0.3 0.3 0.4],'Visible',figures_visibility)
[n1,~]=hist(all_pairs_p_same,xout);
n1=n1./sum(n1);
bar(xout,n1,1)
xlim([0 1])
ylim([0 1])
xlabel('Registered pairs P_s_a_m_e','fontsize',14,'fontweight','bold')
ylabel('Probability','fontsize',14,'fontweight','bold')
x_label=linspace(0,1,6);
x=linspace(0,1,6);
set(gca,'fontsize',14)
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',14,'fontweight','bold')
h=axes('position',[0.3 0.5 0.3 0.3]);
plot(flip(xout),cumsum(flip(n1)),'linewidth',2)
ylim([0 1])
x_label=linspace(0,1,3);
x=linspace(0,1,3);
y=linspace(0,1,3);
y_label=linspace(0,1,3);
set(gca,'YTick',y)
set(gca,'YTickLabel',y_label,'fontsize',14,'fontweight','bold')
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',14,'fontweight','bold')
set(h, 'Xdir', 'reverse')
xlabel('P_s_a_m_e','fontsize',14,'fontweight','bold')
ylabel('Cum. fraction','fontsize',14,'fontweight','bold')
set(gcf,'PaperPositionMode','auto')
savefig(fullfile(figures_directory,'Stage 5 - Registered pairs P_same.fig'))
saveas(gcf,fullfile(figures_directory,'Stage 5 - Registered pairs P_same'),'png')

end

