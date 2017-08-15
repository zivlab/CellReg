function [cell_scores,cell_scores_positive,cell_scores_negative,cell_scores_exclusive,p_same_registered_pairs]=compute_scores(cell_to_index_map,all_to_all_indexes,all_to_all_p_same,number_of_sessions)
% This function computes the false postive, false negative,
% exclusivity, and cell scores for all registered cells according to the
% clustering procedure

% Inputs:
% 1. cell_to_index_map
% 2. all_to_all_indexes
% 3. all_to_all_p_same
% 4. number_of_sessions

% Outputs:
% 1. cell_scores_positive
% 2. cell_scores_negative
% 3. cell_scores_exclusive
% 4. cell_scores
% 5. p_same_registered_pairs

number_of_clusters=size(cell_to_index_map,1);
cell_scores=zeros(1,number_of_clusters);
cell_scores_positive=zeros(1,number_of_clusters);
cell_scores_negative=zeros(1,number_of_clusters);
cell_scores_exclusive=zeros(1,number_of_clusters);
p_same_registered_pairs=cell(1,number_of_clusters);
for n=1:number_of_clusters
    p_same_registered_pairs{n}=nan(number_of_sessions,number_of_sessions);
    % initialize score counts for each cell:
    good_pairs=0;
    good_pairs_positive=0;
    good_pairs_negative=0;
    good_pairs_exclusive=0;
    num_comparisons=0; % the denominator in the calculation of the score
    num_comparisons_positive=0;
    num_comparisons_negative=0;
    cells_in_cluster=cell_to_index_map(n,:);
    for m=1:number_of_sessions
        for k=1:number_of_sessions
            if k~=m & cells_in_cluster(m)>0 % active cell in this session
                this_cell=cell_to_index_map(n,m);
                num_comparisons=num_comparisons+1; 
                cells_to_check=all_to_all_indexes{m}{this_cell,k};
                if cell_to_index_map(n,k)==0 % active-inactive comparison:
                    num_comparisons_negative=num_comparisons_negative+1; 
                    if isempty(cells_to_check) % no candidate means true negative
                        good_pairs=good_pairs+1;
                        good_pairs_negative=good_pairs_negative+1;
                    else
                        this_p_same=all_to_all_p_same{m}{this_cell,k}; % true negative scores
                        if ~isempty(this_p_same)
                            good_pairs=good_pairs+1-sum(this_p_same);
                            good_pairs_negative=good_pairs_negative+1-sum(this_p_same);
                        else
                            good_pairs=good_pairs+1;
                            good_pairs_negative=good_pairs_negative+1;
                        end
                    end
                else % active-active comparison:
                    num_comparisons_positive=num_comparisons_positive+1;
                    this_p_same=all_to_all_p_same{m}{this_cell,k};
                    clustered_cell=cell_to_index_map(n,k);
                    clustered_ind=find(cells_to_check==clustered_cell);
                    if ~isempty(this_p_same) & ~isempty(clustered_ind)
                        temp_true_positive=this_p_same(clustered_ind);
                        p_same_registered_pairs{n}(m,k)=temp_true_positive;
                        good_pairs_positive=good_pairs_positive+temp_true_positive;
                        this_p_same(clustered_ind)=[];
                    else
                        temp_true_positive=0;
                        good_pairs_positive=good_pairs_positive+temp_true_positive;
                    end
                    if isempty(this_p_same) % cell score and exclusivity score
                        good_pairs=good_pairs+temp_true_positive;
                        good_pairs_exclusive=good_pairs_exclusive+1;
                    else
                        good_pairs=good_pairs+temp_true_positive-sum(this_p_same);
                        good_pairs_exclusive=good_pairs_exclusive+1-sum(this_p_same);
                    end
                end
            end
        end
    end
    % normalizing all the scores by the numbers of comparisons:       
    if num_comparisons_positive>0
        cell_scores_positive(n)=good_pairs_positive/num_comparisons_positive;
        cell_scores_exclusive(n)=good_pairs_exclusive/num_comparisons_positive;
    else
        cell_scores_positive(n)=nan;
        cell_scores_exclusive(n)=nan;
    end
    if num_comparisons_negative>0
        cell_scores_negative(n)=good_pairs_negative/num_comparisons_negative;
    else
        cell_scores_negative(n)=nan;
    end
    if num_comparisons>0
        cell_scores(n)=good_pairs/num_comparisons;
    else
        cell_scores(n)=nan;
    end
end

end

