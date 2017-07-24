function plot_single_session_projections(footprints_projection,session_number,figures_visibility)
% This function plots the projections of all the spatial footprints
% for all the sessions.

% Inputs:
% 1. footprints_projection - projection of all spatial footprints 
% 2. session_number - index of the session
% 3. figures_visibility - either 'on' or 'off'

figure('Visible',figures_visibility)
set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')')
imagesc(footprints_projection{1})
set(gca,'xtick',[])
set(gca,'ytick',[])
colormap('gray')
title(['Session ' num2str(session_number)],'fontsize',14,'fontweight','bold')

end

