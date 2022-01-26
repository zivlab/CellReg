function plot_alignment_results(spatial_footprints,centroid_locations,spatial_footprints_corrected,centroid_locations_corrected,footprints_projections,footprints_projections_corrected,reference_session_index,all_centroid_projections_correlations,maximal_cross_correlation,best_translations,overlapping_FOV,alignment_type,number_of_cells_per_session,figures_directory,figures_visibility,varargin)
% This function plots all the results for the image alignment step

number_of_sessions=size(footprints_projections,2);
adjusted_x_size=size(spatial_footprints{1},3);
adjusted_y_size=size(spatial_footprints{1},2);
registration_order=setdiff(1:number_of_sessions,reference_session_index);

best_x_translation=best_translations(1,:);
best_y_translation=best_translations(2,:);
if size(best_translations,1)==3
    best_rotation=best_translations(3,:);
end

% plotting 3 sessions by RGB overlay before and after alignment:
rgb_ind=[1 2 3];
footprints_projections_rgb=zeros(adjusted_y_size,adjusted_x_size,3);
footprints_projections_rgb(:,:,1)=footprints_projections{rgb_ind(1)};
footprints_projections_rgb(:,:,2)=footprints_projections{rgb_ind(2)};
if number_of_sessions>2
    footprints_projections_rgb(:,:,3)=footprints_projections{rgb_ind(3)};
end
footprints_projections_rgb(footprints_projections_rgb>1)=1;

footprints_projections_corrected_rgb=zeros(adjusted_y_size,adjusted_x_size,3);
footprints_projections_corrected_rgb(:,:,1)=footprints_projections_corrected{rgb_ind(1)};
footprints_projections_corrected_rgb(:,:,2)=footprints_projections_corrected{rgb_ind(2)};
if number_of_sessions>2
    footprints_projections_corrected_rgb(:,:,3)=footprints_projections_corrected{rgb_ind(3)};
end
footprints_projections_corrected_rgb(footprints_projections_corrected_rgb>1)=1;
footprints_projections_corrected_rgb=footprints_projections_corrected_rgb+0.25*max(max(max(footprints_projections_corrected_rgb)))*repmat(overlapping_FOV,1,1,3);

figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8],'Visible',figures_visibility)
set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')')
subplot(2,2,1)
imshow(footprints_projections_rgb)
if number_of_sessions>2
    title('RGB overlay (sessions 1-3): Pre-alignment','FontWeight','Bold','fontsize',14)
else
    title('RGB overlay: Pre-alignment','FontWeight','Bold','fontsize',14)
end
subplot(2,2,2)
imshow(footprints_projections_corrected_rgb)
if number_of_sessions>3
    title('RGB overlay (sessions 1-3): Post-alignment','FontWeight','Bold','fontsize',14)
else
    title('RGB overlay: Post-alignment','FontWeight','Bold','fontsize',14)
end
legend_strings=cell(1,number_of_sessions);
if number_of_sessions>16
    color=rand(number_of_sessions,3);
else
    color=[1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 0; 0 0 0; 0.4 0.7 0.4; 1 0.5 0; 0.8 0.8 0.8;...
        0.2 0.1 0.5; 0.4 0.3 0.2; 0.7 0 0.3; 0.6 0.9 0.2 ; 0.1 0.2 0.3; 0.3 0.8 0.3];
end

% plotting the centroid locations of all cells from all sessions on the same FOV:
subplot(2,2,3)
for n=1:number_of_sessions
    centroids=centroid_locations{n};
    h=scatter(centroids(:,1),adjusted_y_size-centroids(:,2),10);
    set(h,'MarkerFaceColor',color(n,:));
    hold on
    legend_strings{n}=['S. ' num2str(n)];
end
title('Centroid locations: Pre-alignment','FontWeight','Bold','FontSize',14)
xlim([0 size(spatial_footprints{1},3)])
ylim([0 size(spatial_footprints{1},2)])
set(gca,'xtick',[])
set(gca,'ytick',[])
legend(legend_strings)
legend('boxoff')

subplot(2,2,4)
for n=1:number_of_sessions
    centroids=centroid_locations_corrected{n};
    h=scatter(centroids(:,1),adjusted_y_size-centroids(:,2),15);
    set(h,'MarkerFaceColor',color(n,:));
    hold on
end
title('Centroid locations: Post-alignment','FontWeight','Bold','FontSize',14)
xlim([0 size(spatial_footprints_corrected{1},3)])
ylim([0 size(spatial_footprints_corrected{1},2)])
set(gca,'xtick',[])
set(gca,'ytick',[])
legend(legend_strings)
legend('boxoff')
set(gcf,'PaperPositionMode','auto')
savefig(fullfile(figures_directory,'Stage 2 - pre vs post alignment.fig'))
saveas(gcf,fullfile(figures_directory,'Stage 2 - pre vs post alignment'),'png')

% Plotting measurments of preparation stability:
if number_of_sessions>2
    % plottig the FOV correlations between all pairs of sessions:
    figure('units','normalized','outerposition',[0.2 0.3 0.6 0.4],'Visible',figures_visibility)
    set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')')
    subplot(1,2,1)
    imagesc(all_centroid_projections_correlations)
    colormap('jet')
    colorbar
    caxis([0 1])
    x_label=cell(1,number_of_sessions);
    if number_of_sessions<8
        for n=1:number_of_sessions
            x_label{n}=num2str(n);
        end
    else
        for n=1:3:number_of_sessions
            x_label{n}=num2str(n);
        end
    end
    x=1:number_of_sessions;
    y_label=cell(1,number_of_sessions);
    if number_of_sessions<8
        for n=1:number_of_sessions
            y_label{n}=num2str(n);
        end
    else
        for n=1:3:number_of_sessions
            y_label{n}=num2str(n);
        end
    end
    y=1:number_of_sessions;
    hold on
    plot([0.5 0.5 number_of_sessions+0.5 number_of_sessions+0.5 0.5],[reference_session_index-0.5 reference_session_index+0.5 reference_session_index+0.5 reference_session_index-0.5 reference_session_index-0.5],'linewidth',3,'color','k')
    set(gca,'fontsize',14)
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',14,'fontweight','bold')
    set(gca,'YTick',y)
    set(gca,'YTickLabel',y_label,'fontsize',14,'fontweight','bold')
    title('Maximal correlation','fontsize',14,'fontweight','bold')
    xlabel('Session number','FontWeight','Bold','fontsize',14)
    ylabel('Session number','FontWeight','Bold','fontsize',14)
    set(gca,'fontsize',14)
    subplot(1,2,2)
    for n=1:number_of_sessions
        ind_to_plot=setdiff(1:number_of_sessions,n);
        plot(1:number_of_sessions-1,all_centroid_projections_correlations(n,ind_to_plot),'*--','linewidth',2,'markersize',8,'color',color(n,:));
        hold on
    end
    ylim([0 1])
    xlim([0 number_of_sessions])
    legend_strings=cell(1,number_of_sessions);
    for n=1:number_of_sessions
        legend_strings{n}=['Ref. S. - ' num2str(n)];
    end
    legend(legend_strings)
    legend('boxoff')
    if number_of_sessions<8
        for n=1:number_of_sessions
            x_label{n}=num2str(n);
        end
    else
        for n=1:3:number_of_sessions
            x_label{n}=num2str(n);
        end
    end
    x=1:number_of_sessions-1;
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',14,'fontweight','bold')
    set(gca,'fontsize',14)
    xlabel('Session number','FontWeight','Bold','fontsize',14)
    ylabel('Maximal correlation','FontWeight','Bold','fontsize',14)
    set(gcf,'PaperPositionMode','auto')
    savefig(fullfile(figures_directory,'Stage 2 - abnormalities test - Correlations.fig'))
    saveas(gcf,fullfile(figures_directory,'Stage 2 - abnormalities test - Correlations'),'png')
end

figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8],'Visible',figures_visibility)
set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')')
% plotting the translations for each session compared to the reference
subplot(2,2,1)
for n=1:number_of_sessions
    plot(best_x_translation(n),best_y_translation(n),'*','markersize',8,'linewidth',2,'color',color(n,:));
    hold on
end
legend(legend_strings,'Location','NorthWest')
legend('boxoff')
x_for_circle=0:0.01:2*pi;
plot(25*sin(x_for_circle),25*cos(x_for_circle),'--','color','k','linewidth',1);
hold on
plot(50*sin(x_for_circle),50*cos(x_for_circle),'--','color','k','linewidth',1);
hold on
plot(75*sin(x_for_circle),75*cos(x_for_circle),'--','color','k','linewidth',1);
hold on
plot(100*sin(x_for_circle),100*cos(x_for_circle),'--','color','k','linewidth',1);
ylim([-100 100])
xlim([-100 100])
axis square
x_label=-100:50:100;
x=-100:50:100;
y=-100:50:100;
y_label=-100:50:100;
set(gca,'YTick',y)
set(gca,'YTickLabel',y_label,'fontsize',14,'fontweight','bold')
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',14,'fontweight','bold')
xlabel('Y translation (\mum)','FontWeight','Bold','fontsize',14)
ylabel('X translation (\mum)','FontWeight','Bold','fontsize',14)

mean_number_of_cells=mean(number_of_cells_per_session);
is_high_number_of_cells=zeros(1,number_of_sessions);
is_low_number_of_cells=zeros(1,number_of_sessions);
for n=1:number_of_sessions
    if number_of_cells_per_session(n)>1.5*mean_number_of_cells
        is_high_number_of_cells(n)=1;
    end
    if number_of_cells_per_session(n)<2/3*mean_number_of_cells
        is_low_number_of_cells(n)=1;
    end
end

if strcmp(alignment_type,'Translations and Rotations')
    % plotting the rotations for each session compared to the reference
    subplot(2,2,2)
    plot(1:number_of_sessions,best_rotation,'*','linewidth',2,'markersize',8,'color','b')
    xlim([0 number_of_sessions+1])
    ylim([-30 30])
    hold on
    plot([0 number_of_sessions+1],[-10 -10],'--','color','k','linewidth',2)
    hold on
    plot([0 number_of_sessions+1],[-20 -20],'--','color','k','linewidth',2)
    hold on
    plot([0 number_of_sessions+1],[10 10],'--','color','k','linewidth',2)
    hold on
    plot([0 number_of_sessions+1],[20 20],'--','color','k','linewidth',2)
    if number_of_sessions<8
        x_label=cell(1,number_of_sessions);
        for n=1:number_of_sessions
            x_label{n}=num2str(n);
        end
    else
        x_label=cell(1,length(1:3:number_of_sessions));
        for n=1:3:number_of_sessions
            x_label{n}=num2str(n);
        end
    end
    x_label{reference_session_index}='Ref.';
    x=1:number_of_sessions;
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',14,'fontweight','bold')
    y_label=-30:10:30;
    y=-30:10:30;
    set(gca,'YTick',y)
    set(gca,'YTickLabel',y_label,'fontsize',14,'fontweight','bold')
    xlabel('Session number','FontWeight','Bold','fontsize',14)
    ylabel('Rotation (degrees)','FontWeight','Bold','fontsize',14)
    set(gca,'fontsize',14)
end

% plotting the FOV correlation between each session and the reference
subplot(2,2,3)
plot(registration_order,maximal_cross_correlation,'*--','markersize',8,'linewidth',2,'color','b')
hold on
plot([1 number_of_sessions],[mean(maximal_cross_correlation) mean(maximal_cross_correlation)],'--','linewidth',2,'color','r')
xlim([0 number_of_sessions+1])
ylim([0 max([2*mean(maximal_cross_correlation),max(maximal_cross_correlation)])])
x_label=cell(1,number_of_sessions);
if number_of_sessions<8
    for n=1:number_of_sessions
        x_label{n}=num2str(n);
    end
else
    for n=1:3:number_of_sessions
        x_label{n}=num2str(n);
    end
end
x_label{reference_session_index}='Ref.';
x=1:number_of_sessions;
set(gca,'fontsize',14)
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',14,'fontweight','bold')
xlabel('Session number','FontWeight','Bold','fontsize',14)
ylabel('Max. projection correlation','FontWeight','Bold','fontsize',14)
set(gca,'fontsize',14)

% plotting the number of identified neurons for each session
subplot(2,2,4)
plot(1:number_of_sessions,number_of_cells_per_session,'*--','markersize',8,'linewidth',2,'color','b')
hold on
plot([1 number_of_sessions],[mean_number_of_cells mean_number_of_cells],'--','linewidth',2,'color','r')
xlim([0 number_of_sessions+1])
ylim([0 max([2*mean_number_of_cells,max(number_of_cells_per_session)])])
xlabel('Session number','fontsize',14,'fontweight','bold')
ylabel('Number of cells','fontsize',14,'fontweight','bold')
if number_of_sessions<8
    x_label=1:number_of_sessions;
else
    x_label=cell(1,number_of_sessions);
    for n=1:3:number_of_sessions
        x_label{n}=num2str(n);
    end
end
x=1:number_of_sessions;
set(gca,'XTick',x)
set(gca,'XTickLabel',x_label,'fontsize',14,'fontweight','bold')
set(gca,'fontsize',14)
set(gcf,'PaperPositionMode','auto')
savefig(fullfile(figures_directory,'Stage 2 - abnormalities test - general.fig'))
saveas(gcf,fullfile(figures_directory,'Stage 2 - abnormalities test - general'),'png')

% if non-rigid transformation was used:
if strcmp(alignment_type,'Non-rigid')
    displacement_fields=varargin{1};
    adjusted_x_size=size(displacement_fields,3);
    adjusted_y_size=size(displacement_fields,2);
    number_of_arrows=25;
    x_vector=round(adjusted_x_size/number_of_arrows/2):round(adjusted_x_size/number_of_arrows):adjusted_x_size;
    y_vector=round(adjusted_y_size/number_of_arrows/2):round(adjusted_y_size/number_of_arrows):adjusted_y_size;
    [x,y]=meshgrid(x_vector,y_vector);
    subx=4;
    suby=ceil(number_of_sessions/subx);
    if number_of_sessions>4
        figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8],'Visible',figures_visibility)
        set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')')
        for n=1:number_of_sessions
            subplot(suby,subx,n)
            quiver(x,y,squeeze(displacement_fields(n,(y_vector),(x_vector),1)),squeeze(flipud(displacement_fields(n,(y_vector),(x_vector),2))))
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            xlim([0 adjusted_x_size])
            ylim([0 adjusted_y_size])
            colormap('gray')
            title(['Session ' num2str(n)],'fontsize',14,'fontweight','bold')
        end
    else
        figure('units','normalized','outerposition',[0.1 0.2 0.8 0.5],'Visible',figures_visibility)
        set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')')
        for n=1:number_of_sessions
            subplot(1,number_of_sessions,n)
            quiver(x,y,squeeze(displacement_fields(n,(y_vector),(x_vector),1)),squeeze(flipud(displacement_fields(n,(y_vector),(x_vector),2))))
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            xlim([0 adjusted_x_size])
            ylim([0 adjusted_y_size])
            colormap('gray')
            title(['Session ' num2str(n)],'fontsize',14,'fontweight','bold')
        end
    end
    set(gcf,'PaperPositionMode','auto')
    savefig(fullfile(figures_directory,'Stage 1 - Non-rigid transformations.fig'))
    saveas(gcf,fullfile(figures_directory,'Stage 1 - Non-rigid transformations'),'png')    
end
    
end

