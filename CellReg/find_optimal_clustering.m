function [cell_to_index_map,mean_cell_scores,mean_cell_scores_positive,mean_cell_scores_negative,mean_cell_scores_exclusive,cell_scores,changes_count,switch_count,separation_count,move_count,delete_count,decision_type,cell_scores_positive,cell_scores_negative,cell_scores_exclusive,all_clusters_centroids]=find_optimal_clustering(cell_to_index_map,all_to_all_p_value_multi,all_to_all_matrix_multi,max_iterations,cluster_distance_thresh,correlation_thresh,all_centroids_corrected,num_sessions,decision_type,~,is_figure,results_dir,figures_dir)

% Finding the optimal clustering with an iterative process:
num_changes_thresh=10;
num_clusters=size(cell_to_index_map,1);
cell_scores=zeros(1,num_clusters);
cell_scores_positive=zeros(1,num_clusters);
cell_scores_negative=zeros(1,num_clusters);
cell_scores_exclusive=zeros(1,num_clusters);
active_sessions=zeros(1,num_clusters);
iteration=1;
for n=1:num_clusters
    good_pairs=0;
    good_pairs_positive=0;
    good_pairs_negative=0;
    good_pairs_exclusive=0;
    num_comparisons=0;
    num_comparisons_positive=0;
    num_comparisons_negative=0;
    cells_in_cluster=cell_to_index_map(n,:);
    for m=1:num_sessions
        for k=1:num_sessions
            if k~=m && cells_in_cluster(m)>0
                this_cell=cell_to_index_map(n,m);
                num_comparisons=num_comparisons+1;
                cells_to_check=all_to_all_matrix_multi{m}{this_cell,k};
                if cell_to_index_map(n,k)==0
                    num_comparisons_negative=num_comparisons_negative+1;
                    if isempty(cells_to_check)
                        good_pairs=good_pairs+1;
                        good_pairs_negative=good_pairs_negative+1;
                    else
                        this_p_value=all_to_all_p_value_multi{m}{this_cell,k};
                        if max(this_p_value)<0.05
                            good_pairs=good_pairs+1;
                            good_pairs_negative=good_pairs_negative+1;
                        end
                    end
                else
                    num_comparisons_positive=num_comparisons_positive+1;
                    this_p_value=all_to_all_p_value_multi{m}{this_cell,k};
                    clustered_cell=cell_to_index_map(n,k);
                    clustered_ind=find(cells_to_check==clustered_cell);
                    if this_p_value(clustered_ind)>0.95
                        this_p_value(clustered_ind)=[];
                        good_pairs_positive=good_pairs_positive+1;
                        if isempty(this_p_value)
                            good_pairs=good_pairs+1;
                            good_pairs_exclusive=good_pairs_exclusive+1;
                        else
                            if max(this_p_value)<0.05
                                good_pairs=good_pairs+1;
                                good_pairs_exclusive=good_pairs_exclusive+1;
                            end
                        end
                    else
                        this_p_value(clustered_ind)=[];
                        if isempty(this_p_value)
                            good_pairs_exclusive=good_pairs_exclusive+1;
                        else
                            if max(this_p_value)<0.05
                                good_pairs_exclusive=good_pairs_exclusive+1;
                            end
                        end
                    end
                end
            end
        end
    end
    cell_scores_positive(n)=good_pairs_positive/num_comparisons_positive;
    cell_scores_negative(n)=good_pairs_negative/num_comparisons_negative;
    cell_scores_exclusive(n)=good_pairs_exclusive/num_comparisons_positive;
    cell_scores(n)=good_pairs/num_comparisons;
    active_sessions(n)=sum(cells_in_cluster>0);
end
if is_figure
    xout_temp=linspace(0,1,num_sessions*2+1);
    xout=xout_temp(2:2:end);
    figure
    fig_size_y=20;
    fig_size_x=25;
    set(gcf,'PaperUnits','centimeters','PaperPosition',[1 5 fig_size_x fig_size_y]);
    set(gcf,'PaperOrientation','portrait');
    set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[5 0 0 0]);
    size_x=0.65;
    size_y=0.65;
    
    axes('position',[0.6 0.1 size_x/2 size_y/2])
    [n1,~]=hist(cell_scores,xout);
    n1=n1./sum(n1);
    bar(xout,n1)
    xlim([0 1])
    ylim([0 1])
    xlabel('Overall cell scores','fontsize',16,'fontweight','bold')
    ylabel('Probability','fontsize',16,'fontweight','bold')
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
    
    axes('position',[0.12 0.1 size_x/2 size_y/2])
    [n1,~]=hist(cell_scores_exclusive,xout);
    n1=n1./sum(n1);
    bar(xout,n1)
    xlim([0 1])
    ylim([0 1])
    xlabel('Exclusivity cell scores','fontsize',16,'fontweight','bold')
    ylabel('Probability','fontsize',16,'fontweight','bold')
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
    
    axes('position',[0.6 0.58 size_x/2 size_y/2])
    [n1,~]=hist(cell_scores_positive,xout);
    n1=n1./sum(n1);
    bar(xout,n1)
    xlim([0 1])
    ylim([0 1])
    xlabel('True positive scores','fontsize',16,'fontweight','bold')
    ylabel('Probability','fontsize',16,'fontweight','bold')
    x_label=linspace(0,1,6);
    x=linspace(0,1,6);
    set(gca,'fontsize',14)
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',14,'fontweight','bold')
    text(-0.25,1.2,['Iteration number ' num2str(iteration) ': ' num2str(num_clusters) ' registered cells'],'fontsize',24,'fontweight','bold','HorizontalAlignment','Center')
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
    
    axes('position',[0.12 0.58 size_x/2 size_y/2])
    [n1,~]=hist(cell_scores_negative,xout);
    n1=n1./sum(n1);
    bar(xout,n1)
    xlim([0 1])
    ylim([0 1])
    xlabel('True negative scores','fontsize',16,'fontweight','bold')
    ylabel('Probability','fontsize',16,'fontweight','bold')
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
end

mean_cell_scores=zeros(1,max_iterations);
mean_cell_scores(1)=mean(cell_scores);
mean_cell_scores_positive=zeros(1,max_iterations);
mean_cell_scores_positive(1)=mean(cell_scores_positive(cell_scores_positive>=0));
mean_cell_scores_negative=zeros(1,max_iterations);
mean_cell_scores_negative(1)=mean(cell_scores_negative(cell_scores_negative>=0));
mean_cell_scores_exclusive=zeros(1,max_iterations);
mean_cell_scores_exclusive(1)=mean(cell_scores_exclusive(cell_scores_exclusive>=0));

changes_count=zeros(1,max_iterations);
switch_count=zeros(1,max_iterations);
separation_count=zeros(1,max_iterations);
move_count=zeros(1,max_iterations);
delete_count=zeros(1,max_iterations);
changes_count(1)=-1;

hbar = waitbar(0,'Clustering cells according to the probability model');
while (changes_count(iteration)>num_changes_thresh || changes_count(iteration)==-1) && (iteration)<max_iterations
    waitbar(iteration/max_iterations,hbar ,['Performing cell clustering - iteration number ' num2str(iteration) '/' num2str(max_iterations)])
    
    % Finding the center of mass of all clusters:
    iteration=iteration+1;
    num_clusters=size(cell_to_index_map,1);
    all_clusters_centroids=zeros(2,num_clusters);
    for n=1:num_clusters
        cells_in_cluster=cell_to_index_map(n,:);
        sessions_ind=find(cells_in_cluster>0);
        N=length(sessions_ind);
        centroids_found=zeros(2,N);
        for k=1:N
            centroids_found(:,k)=all_centroids_corrected{sessions_ind(k)}(cells_in_cluster(sessions_ind(k)),:);
        end
        all_clusters_centroids(:,n)=mean(centroids_found,2);
    end
    
    % Checking for each cell which cluster is the most correlated:
    for n=1:num_sessions
        num_cells=sum(cell_to_index_map(:,n)>0);
        cluster_ind=find(cell_to_index_map(:,n)>0);
        for k=1:num_cells
            this_cell=cell_to_index_map(cluster_ind(k),n);
            this_cell_centroid=all_centroids_corrected{n}(this_cell,:);
            centroid=repmat(this_cell_centroid,num_clusters,1)';
            distance_vec=sqrt(sum((centroid-all_clusters_centroids).^2));
            clusters_to_check=find(distance_vec<cluster_distance_thresh);
            num_candidates=length(clusters_to_check);
            total_correlation=zeros(1,num_candidates);
            normalization_factor=zeros(1,num_candidates);
            max_in_cluster_correlation=zeros(1,num_candidates);
            min_in_cluster_correlation=ones(1,num_candidates);
            for m=1:num_candidates
                this_cluster=clusters_to_check(m);
                cells_in_cluster=cell_to_index_map(this_cluster,:);
                % checking if cell from this session is already in the cluster
                sessions_in_cluster=find(cells_in_cluster>0);
                if sum(sessions_in_cluster-n==0)>0 % cell from same session already in cluster
                    sessions_in_cluster(sessions_in_cluster-n==0)=[];
                end
                if ~isempty(sessions_in_cluster)
                    for l=1:length(sessions_in_cluster)
                        in_cluster=cell_to_index_map(this_cluster,sessions_in_cluster(l));
                        in_session=all_to_all_matrix_multi{n}{this_cell,sessions_in_cluster(l)};
                        in_cluster_ind=find(in_session==in_cluster);
                        if ~isempty(all_to_all_p_value_multi{n}{this_cell,sessions_in_cluster(l)}(in_cluster_ind))
                            this_correlation=all_to_all_p_value_multi{n}{this_cell,sessions_in_cluster(l)}(in_cluster_ind);
                            total_correlation(m)=total_correlation(m)+this_correlation;
                            normalization_factor(m)=normalization_factor(m)+1;
                            max_in_cluster_correlation(m)=max(max_in_cluster_correlation(m),this_correlation);
                            min_in_cluster_correlation(m)=min(min_in_cluster_correlation(m),this_correlation);
                        else
                            min_in_cluster_correlation(m)=0;
                        end
                    end
                else
                    min_in_cluster_correlation(m)=0;
                end
            end
            
            if strcmp(decision_type,'Minimal correlation')
                [max_correlation,max_correlation_ind]=max(min_in_cluster_correlation);
            elseif strcmp(decision_type,'Average correlation')
                normalized_correlation=total_correlation./(normalization_factor);
                [max_correlation,max_correlation_ind]=max(normalized_correlation);
                average_correlation=max_correlation;
            elseif strcmp(decision_type,'Maximal correlation')
                [max_correlation,max_correlation_ind]=max(max_in_cluster_correlation);
            end
            
            max_cluster=clusters_to_check(max_correlation_ind);
            cells_in_max_cluster=cell_to_index_map(max_cluster,:);
            sessions_in_max_cluster=find(cells_in_max_cluster>0);
            if sum(sessions_in_max_cluster-n==0)>0 % cell from same session already in cluster
                sessions_in_max_cluster(sessions_in_max_cluster-n==0)=[];
            end
            
            if ~isempty(sessions_in_max_cluster)
                cells_in_original_cluster=cell_to_index_map(cluster_ind(k),:);
                num_cells_in_original_cluster=sum(cells_in_original_cluster>0);
                
                if strcmp(decision_type,'Minimal correlation')
                    average_correlation=max_correlation;
                elseif strcmp(decision_type,'Maximal correlation')
                    average_correlation=max_correlation;
                end
                
                if average_correlation<correlation_thresh && num_cells_in_original_cluster>1 % split the cell to a new cluster
                    num_clusters=num_clusters+1;
                    cell_to_index_map(num_clusters,:)=zeros(1,num_sessions);
                    cell_to_index_map(num_clusters,n)=this_cell;
                    cell_to_index_map(cluster_ind(k),n)=0;
                    all_clusters_centroids(:,num_clusters)=this_cell_centroid;
                    changes_count(iteration)=changes_count(iteration)+1;
                    separation_count(iteration)=separation_count(iteration)+1;
                elseif average_correlation>=correlation_thresh
                    if (clusters_to_check(max_correlation_ind)~=cluster_ind(k)) % need to change clusters
                        this_cluster=clusters_to_check(max_correlation_ind);
                        cells_in_cluster=cell_to_index_map(this_cluster,:);
                        % checking if cell from this session is already in the cluster
                        sessions_in_cluster=find(cells_in_cluster>0);
                        if sum(sessions_in_cluster-n==0)>0 % cell from same session already in cluster - switch
                            temp_cell=cell_to_index_map(this_cluster,n);
                            temp_cell_centroid=all_centroids_corrected{n}(temp_cell,:);
                            temp_correlation=0;
                            % switch only if this cell is more correlated to the cluster
                            cells_in_cluster=cell_to_index_map(this_cluster,:);
                            % checking if cell from this session is already in the cluster
                            sessions_in_cluster=find(cells_in_cluster>0);
                            if sum(sessions_in_cluster-n==0)>0 % cell from same session already in cluster
                                sessions_in_cluster(sessions_in_cluster-n==0)=[];
                            end
                            for l=1:length(sessions_in_cluster)
                                in_cluster=cell_to_index_map(this_cluster,sessions_in_cluster(l));
                                in_session=all_to_all_matrix_multi{n}{temp_cell,sessions_in_cluster(l)};
                                in_cluster_ind=find(in_session==in_cluster);
                                if ~isempty(all_to_all_p_value_multi{n}{temp_cell,sessions_in_cluster(l)}(in_cluster_ind))
                                    temp_correlation=temp_correlation+all_to_all_p_value_multi{n}{temp_cell,sessions_in_cluster(l)}(in_cluster_ind);
                                end
                            end
                            if max_correlation>temp_correlation
                                num_clusters=num_clusters+1;
                                cell_to_index_map(this_cluster,n)=this_cell;
                                cell_to_index_map(num_clusters,:)=zeros(1,num_sessions);
                                cell_to_index_map(num_clusters,n)=temp_cell;
                                changes_count(iteration)=changes_count(iteration)+1;
                                switch_count(iteration)=switch_count(iteration)+1;
                                cell_to_index_map(cluster_ind(k),n)=0;
                                all_clusters_centroids(:,num_clusters)=temp_cell_centroid;
                            end
                        else % just need to add the cell to the cluster
                            cell_to_index_map(cluster_ind(k),n)=0;
                            cell_to_index_map(this_cluster,n)=this_cell;
                            changes_count(iteration)=changes_count(iteration)+1;
                            move_count(iteration)=move_count(iteration)+1;
                        end
                    end
                end
            end
        end
    end
    
    % deleting clusters with no cells left
    cell_to_index_map_temp=cell_to_index_map;
    num_clusters=size(cell_to_index_map_temp,1);
    for n=1:num_clusters
        if sum(cell_to_index_map_temp(n,:))==0
            cell_to_index_map(n-delete_count(iteration),:)=[];
            delete_count(iteration)=delete_count(iteration)+1;
        end
    end
    cell_to_index_map_temp=cell_to_index_map;
    
    % Merging clusters if they cross threshold:
    num_clusters=size(cell_to_index_map_temp,1);
    all_clusters_centroids=zeros(2,num_clusters);
    for n=1:num_clusters
        cells_in_cluster=cell_to_index_map_temp(n,:);
        sessions_ind=find(cells_in_cluster>0);
        N=length(sessions_ind);
        centroids_found=zeros(2,N);
        for k=1:N
            centroids_found(:,k)=all_centroids_corrected{sessions_ind(k)}(cells_in_cluster(sessions_ind(k)),:);
        end
        all_clusters_centroids(:,n)=mean(centroids_found,2);
    end
    
    for n=1:num_clusters
        this_cluster_cells=cell_to_index_map_temp(n,:);
        this_cluster_sessions=find(this_cluster_cells>0);
        centroid=repmat(all_clusters_centroids(:,n),1,num_clusters);
        distance_vec=sqrt(sum((centroid-all_clusters_centroids).^2));
        clusters_to_check=find(distance_vec<cluster_distance_thresh);
        clusters_to_check=setdiff(clusters_to_check,n);
        num_candidates=length(clusters_to_check);
        for k=1:num_candidates
            candidate_cells=cell_to_index_map_temp(clusters_to_check(k),:);
            candidate_sessions=find(candidate_cells>0);
            if isempty(intersect(find(this_cluster_cells>0),find(candidate_cells>0)))
                all_to_all_temp=zeros(sum(this_cluster_cells>0),sum(candidate_cells>0));
                for m=1:sum(this_cluster_cells>0);
                    for l=1:sum(candidate_cells>0)
                        neigbor_cell=all_to_all_matrix_multi{this_cluster_sessions(m)}{this_cluster_cells(this_cluster_sessions(m)),candidate_sessions(l)};
                        this_cell_ind=find(neigbor_cell==candidate_cells(candidate_sessions(l)));
                        if ~isempty(all_to_all_p_value_multi{this_cluster_sessions(m)}{this_cluster_cells(this_cluster_sessions(m)),candidate_sessions(l)}(this_cell_ind))
                            all_to_all_temp(m,l)=all_to_all_p_value_multi{this_cluster_sessions(m)}{this_cluster_cells(this_cluster_sessions(m)),candidate_sessions(l)}(this_cell_ind);
                        end
                    end
                end
                if strcmp(decision_type,'Average correlation')
                    average_all_to_all_temp_1=mean(all_to_all_temp);
                    average_all_to_all_temp_2=mean(all_to_all_temp,2);
                    over_thresh_1=sum(average_all_to_all_temp_1>correlation_thresh);
                    over_thresh_2=sum(average_all_to_all_temp_2>correlation_thresh);
                    if over_thresh_2==sum(this_cluster_cells>0) && over_thresh_1==sum(candidate_cells>0)
                        cell_to_index_map_temp(n,candidate_sessions)=cell_to_index_map_temp(clusters_to_check(k),candidate_sessions);
                        cell_to_index_map_temp(clusters_to_check(k),:)=zeros(1,num_sessions);
                        changes_count(iteration)=changes_count(iteration)+1;
                        all_clusters_centroids(:,clusters_to_check(k))=[1000 1000];
                    end
                elseif strcmp(decision_type,'Maximal correlation')
                    max_all_to_all_temp=max(all_to_all_temp(:));
                    if max_all_to_all_temp>correlation_thresh;
                        cell_to_index_map_temp(n,candidate_sessions)=cell_to_index_map_temp(clusters_to_check(k),candidate_sessions);
                        cell_to_index_map_temp(clusters_to_check(k),:)=zeros(1,num_sessions);
                        changes_count(iteration)=changes_count(iteration)+1;
                        all_clusters_centroids(:,clusters_to_check(k))=[1000 1000];
                    end
                elseif strcmp(decision_type,'Minimal correlation')
                    average_all_to_all_temp_1=min(all_to_all_temp);
                    average_all_to_all_temp_2=min(all_to_all_temp,2);
                    over_thresh_1=sum(average_all_to_all_temp_1>correlation_thresh);
                    over_thresh_2=sum(average_all_to_all_temp_2>correlation_thresh);
                    if over_thresh_2==sum(this_cluster_cells>0) && over_thresh_1==sum(candidate_cells>0)
                        cell_to_index_map_temp(n,candidate_sessions)=cell_to_index_map_temp(clusters_to_check(k),candidate_sessions);
                        cell_to_index_map_temp(clusters_to_check(k),:)=zeros(1,num_sessions);
                        changes_count(iteration)=changes_count(iteration)+1;
                        all_clusters_centroids(:,clusters_to_check(k))=[1000 1000];
                    end
                end
            end
        end
    end
    cell_to_index_map=cell_to_index_map_temp;
    num_clusters=size(cell_to_index_map_temp,1);
    temp_delete_count=0;
    for n=1:num_clusters
        if sum(cell_to_index_map_temp(n,:))==0
            cell_to_index_map(n-temp_delete_count,:)=[];
            temp_delete_count=temp_delete_count+1;
        end
    end
    delete_count(iteration)=delete_count(iteration)+temp_delete_count;
    
    % Calculating the cell scores
    num_clusters=size(cell_to_index_map,1);
    cell_scores=zeros(1,num_clusters);
    cell_scores_positive=zeros(1,num_clusters);
    cell_scores_negative=zeros(1,num_clusters);
    cell_scores_exclusive=zeros(1,num_clusters);
    active_sessions=zeros(1,num_clusters);
    for n=1:num_clusters
        good_pairs=0;
        good_pairs_positive=0;
        good_pairs_negative=0;
        good_pairs_exclusive=0;
        num_comparisons=0;
        num_comparisons_positive=0;
        num_comparisons_negative=0;
        cells_in_cluster=cell_to_index_map(n,:);
        for m=1:num_sessions
            for k=1:num_sessions
                if k~=m && cells_in_cluster(m)>0
                    this_cell=cell_to_index_map(n,m);
                    num_comparisons=num_comparisons+1;
                    cells_to_check=all_to_all_matrix_multi{m}{this_cell,k};
                    if cell_to_index_map(n,k)==0
                        num_comparisons_negative=num_comparisons_negative+1;
                        if isempty(cells_to_check)
                            good_pairs=good_pairs+1;
                            good_pairs_negative=good_pairs_negative+1;
                        else
                            this_p_value=all_to_all_p_value_multi{m}{this_cell,k};
                            if max(this_p_value)<0.05
                                good_pairs=good_pairs+1;
                                good_pairs_negative=good_pairs_negative+1;
                            end
                        end
                    else
                        num_comparisons_positive=num_comparisons_positive+1;
                        this_p_value=all_to_all_p_value_multi{m}{this_cell,k};
                        clustered_cell=cell_to_index_map(n,k);
                        clustered_ind=find(cells_to_check==clustered_cell);
                        if this_p_value(clustered_ind)>0.95
                            this_p_value(clustered_ind)=[];
                            good_pairs_positive=good_pairs_positive+1;
                            if isempty(this_p_value)
                                good_pairs=good_pairs+1;
                                good_pairs_exclusive=good_pairs_exclusive+1;
                            else
                                if max(this_p_value)<0.05
                                    good_pairs=good_pairs+1;
                                    good_pairs_exclusive=good_pairs_exclusive+1;
                                end
                            end
                        else
                            this_p_value(clustered_ind)=[];
                            if isempty(this_p_value)
                                good_pairs_exclusive=good_pairs_exclusive+1;
                            else
                                if max(this_p_value)<0.05
                                    good_pairs_exclusive=good_pairs_exclusive+1;
                                end
                            end
                        end
                    end
                end
            end
        end
        cell_scores_positive(n)=good_pairs_positive/num_comparisons_positive;
        cell_scores_negative(n)=good_pairs_negative/num_comparisons_negative;
        cell_scores_exclusive(n)=good_pairs_exclusive/num_comparisons_positive;
        cell_scores(n)=good_pairs/num_comparisons;
        active_sessions(n)=sum(cells_in_cluster>0);
    end
    
    mean_cell_scores(iteration)=mean(cell_scores(cell_scores>=0));
    mean_cell_scores_positive(iteration)=mean(cell_scores_positive(cell_scores_positive>=0));
    mean_cell_scores_negative(iteration)=mean(cell_scores_negative(cell_scores_negative>=0));
    mean_cell_scores_exclusive(iteration)=mean(cell_scores_exclusive(cell_scores_exclusive>=0));   
end
close(hbar);

if is_figure
    figure
    fig_size_y=20;
    fig_size_x=25;
    set(gcf,'PaperUnits','centimeters','PaperPosition',[1 5 fig_size_x fig_size_y]);
    set(gcf,'PaperOrientation','portrait');
    set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[5 0 0 0]);
    size_x=0.65;
    size_y=0.65;
    
    axes('position',[0.6 0.1 size_x/2 size_y/2])
    [n1,~]=hist(cell_scores,xout);
    n1=n1./sum(n1);
    bar(xout,n1)
    xlim([0 1])
    ylim([0 1])
    xlabel('Overall cell scores','fontsize',16,'fontweight','bold')
    ylabel('Probability','fontsize',16,'fontweight','bold')
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
    
    axes('position',[0.12 0.1 size_x/2 size_y/2])
    [n1,~]=hist(cell_scores_exclusive,xout);
    n1=n1./sum(n1);
    bar(xout,n1)
    xlim([0 1])
    ylim([0 1])
    xlabel('Exclusivity cell scores','fontsize',16,'fontweight','bold')
    ylabel('Probability','fontsize',16,'fontweight','bold')
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
    
    axes('position',[0.6 0.58 size_x/2 size_y/2])
    [n1,~]=hist(cell_scores_positive,xout);
    n1=n1./sum(n1);
    bar(xout,n1)
    xlim([0 1])
    ylim([0 1])
    xlabel('True positive scores','fontsize',16,'fontweight','bold')
    ylabel('Probability','fontsize',16,'fontweight','bold')
    x_label=linspace(0,1,6);
    x=linspace(0,1,6);
    set(gca,'fontsize',14)
    set(gca,'XTick',x)
    set(gca,'XTickLabel',x_label,'fontsize',14,'fontweight','bold')
    text(-0.25,1.2,['Final registration: ' num2str(num_clusters) ' registered cells'],'fontsize',24,'fontweight','bold','HorizontalAlignment','Center')
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
    
    axes('position',[0.12 0.58 size_x/2 size_y/2])
    [n1,~]=hist(cell_scores_negative,xout);
    n1=n1./sum(n1);
    bar(xout,n1)
    xlim([0 1])
    ylim([0 1])
    xlabel('True negative scores','fontsize',16,'fontweight','bold')
    ylabel('Probability','fontsize',16,'fontweight','bold')
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
    cd(figures_dir);
    savefig('Stage 5 - Register scores')
    saveas(gcf,'Stage 5 - Register scores','tif')
end
cd(results_dir);

% Calculating neighbor vs corr/distance curves
num_clusters=size(cell_to_index_map,1);
all_clusters_centroids=zeros(2,num_clusters);
for n=1:num_clusters
    cells_in_cluster=cell_to_index_map(n,:);
    sessions_ind=find(cells_in_cluster>0);
    N=length(sessions_ind);
    centroids_found=zeros(2,N);
    for k=1:N
        centroids_found(:,k)=all_centroids_corrected{sessions_ind(k)}(cells_in_cluster(sessions_ind(k)),:);
    end
    all_clusters_centroids(:,n)=mean(centroids_found,2);
end

end
