function plot_all_registered_projections(spatial_footprints,cell_to_index_map,figures_directory,figures_visibility,varargin)
% This function plots the projections of all the cells for all the
% sessions. Green cells are those who were active in all sessions.

% Inputs:
% 1. spatial_footprints
% 2. cell_to_index_map
% 3. figures_directory
% 4. figures_visibility
% 5. varargin
%   5{1}. initial_stage

number_of_sessions=size(spatial_footprints,2);

pixel_weight_threshold=0.5; % for better visualization of cells
all_projections_partial=cell(1,number_of_sessions);
mutual_projections_partial=cell(1,number_of_sessions);
cells_in_all_days=find(sum(cell_to_index_map'>0)==number_of_sessions);
other_cells=cell(1,number_of_sessions);
for n=1:number_of_sessions
    logical_1=sum(cell_to_index_map'>0)<number_of_sessions;
    other_cells{n}=find(cell_to_index_map(:,n)'>0 & logical_1);
end
 
disp('Calculating spatial footprints projections:')
for n=1:number_of_sessions
    display_progress_bar('Terminating previous progress bars',true)    
    display_progress_bar(['Calculating projections for session #' num2str(n) ' - '],false)
    this_session_spatial_footprints=spatial_footprints{n};
    num_spatial_footprints=size(this_session_spatial_footprints,1);
    normalized_spatial_footprints=zeros(size(this_session_spatial_footprints));
    for k=1:num_spatial_footprints
        display_progress_bar(100*(k)/(num_spatial_footprints),false)
        this_spatial_footprint=this_session_spatial_footprints(k,:,:);
        this_spatial_footprint(this_spatial_footprint<pixel_weight_threshold*max(max(this_spatial_footprint)))=0;
        if max(max(this_spatial_footprint))>0
            normalized_spatial_footprints(k,:,:)=this_spatial_footprint/max(max(this_spatial_footprint));
        end
    end
    display_progress_bar(' done',false);
    
    all_projections_partial{n}=zeros(size(this_spatial_footprint,2),size(this_spatial_footprint,3),3);
    mutual_projections_partial{n}=zeros(size(this_spatial_footprint,2),size(this_spatial_footprint,3),3);
    all_projections_partial{n}(:,:,2)=squeeze(sum(normalized_spatial_footprints(cell_to_index_map(cells_in_all_days,n),:,:),1));
    all_projections_partial{n}(:,:,1)=squeeze(sum(normalized_spatial_footprints(cell_to_index_map(other_cells{n},n),:,:),1));
    all_projections_partial{n}(:,:,2)=squeeze(sum(normalized_spatial_footprints(cell_to_index_map(other_cells{n},n),:,:),1))+squeeze(sum(normalized_spatial_footprints(cell_to_index_map(cells_in_all_days,n),:,:),1));
    all_projections_partial{n}(:,:,3)=squeeze(sum(normalized_spatial_footprints(cell_to_index_map(other_cells{n},n),:,:),1));
    mutual_projections_partial{n}(:,:,2)=squeeze(sum(normalized_spatial_footprints(cell_to_index_map(cells_in_all_days,n),:,:),1));
    all_projections_partial{n}(all_projections_partial{n}>1)=1;
end

subx=4;
suby=ceil(number_of_sessions/subx);
if number_of_sessions>4
    figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8],'Visible',figures_visibility)
    set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')')
    for n=1:number_of_sessions
        subplot(suby,subx,n)
        imagesc(all_projections_partial{n})
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        colormap('gray')
        title(['Session ' num2str(n)],'fontsize',14,'fontweight','bold')
        if n==1
            text(0.01*size(all_projections_partial{n},1),0.02*size(all_projections_partial{n},2),'Detected in','fontsize',14,'color','g','fontweight','bold')
            text(0.01*size(all_projections_partial{n},1),0.06*size(all_projections_partial{n},2),'all sessions','fontsize',14,'color','g','fontweight','bold')
        end
    end
else
    figure('units','normalized','outerposition',[0.1 0.2 0.8 0.5],'Visible',figures_visibility)
    set(gcf,'CreateFcn','set(gcf,''Visible'',''on'')')
    for n=1:number_of_sessions
        subplot(1,number_of_sessions,n)
        imagesc(all_projections_partial{n})
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        colormap('gray')
        title(['Session ' num2str(n)],'fontsize',14,'fontweight','bold')
        if n==1
            text(0.01*size(all_projections_partial{n},1),0.02*size(all_projections_partial{n},2),'Detected in','fontsize',14,'color','g','fontweight','bold')
            text(0.01*size(all_projections_partial{n},1),0.06*size(all_projections_partial{n},2),'all sessions','fontsize',14,'color','g','fontweight','bold')
        end
    end
end
set(gcf,'PaperPositionMode','auto')

initial_stage=false;
if ~isempty(varargin)
    if varargin{1}
        initial_stage=true;
    end
end

if initial_stage
    savefig(fullfile(figures_directory,'Stage 4 - projcetions - initial registration.fig'))
    saveas(gcf,fullfile(figures_directory,'Stage 4 - projcetions - initial registration'),'png')
else
    savefig(fullfile(figures_directory,'Stage 5 - projcetions - final registration.fig'))
    saveas(gcf,fullfile(figures_directory,'Stage 5 - projcetions - final registration'),'png')    
end

end

