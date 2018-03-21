function plot_all_sessions_projections(footprints_projections,figures_directory,figures_visibility)
% This function plots the projections of all the spatial footprints
% for all the sessions.

% Inputs:
% 1. footprints_projections - projection of all spatial footprints 
% 2. figures directory - where to save the figure
% 3. figures_visibility

number_of_sessions=size(footprints_projections,2);
subx=4;
suby=ceil(number_of_sessions/subx);
if number_of_sessions>4
    figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8],'Visible',figures_visibility)
    set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')')
    for n=1:number_of_sessions
        subplot(suby,subx,n)
        imagesc(footprints_projections{n},[0 2])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        colormap('gray')
        title(['Session ' num2str(n)],'fontsize',14,'fontweight','bold')
    end
else
    figure('units','normalized','outerposition',[0.1 0.2 0.8 0.5],'Visible',figures_visibility)
    set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')')
    for n=1:number_of_sessions
        subplot(1,number_of_sessions,n)
        imagesc(footprints_projections{n},[0 2])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        colormap('gray')
        title(['Session ' num2str(n)],'fontsize',14,'fontweight','bold')
    end
end
set(gcf,'PaperPositionMode','auto')
savefig(fullfile(figures_directory,'Stage 1 - spatial footprints projections.fig'))
saveas(gcf,fullfile(figures_directory,'Stage 1 - spatial footprints projections'),'png')

end

