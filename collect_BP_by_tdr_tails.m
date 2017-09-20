function collect_BP_by_tdr_tails(drug, quantile_used, states)

load('subjects.mat'), load('AP_freqs.mat')
    
state_label = '';

if ~isempty(states)
    
    no_states = length(states);
    
    for state = 1:no_states
        
        state_label = [state_label, '_', states{state}];
        
    end
    
end

name = 'ALL_Frontal';

measure = 'p0.99_IEzs';

BP_drugs = text_read([name,'/',name,'_drugs.txt'],'%s');
BP_subjects = text_read([name,'/',name,'_subjects.txt'],'%s');
BP_states = text_read([name,'/',name,'_states.txt'],'%s');
BP(:, :, 1) = load([name, '/', name, '_BP_pct.txt']);
BP(:, :, 2) = load('ALL_CA1/ALL_CA1_BP_pct.txt');
BP = permute(BP, [1 3 2]);

no_criteria = 4;
no_pairs = nchoosek(no_criteria, 2);
no_figures = no_criteria + no_pairs + 1;

theta_BP = nan(ceil(quantile_used*size(BP, 1)), size(BP, 2), size(BP, 3), no_figures);
non_theta_BP = nan(ceil((1 - quantile_used)*size(BP, 1)), size(BP, 2), size(BP, 3), no_figures);

[mean_subj_tBP, mean_subj_ntBP] = deal(nan(size(BP, 2), size(BP, 3), no_figures, subj_num));

[tBP_marker, ntBP_marker] = deal(zeros(no_figures, 1));

pairs = nchoosek(1:no_criteria, 2);

theta_labels = {'Peak Power \delta/\theta', 'Band Power \delta/\theta', '\delta Power', '\theta Power'};

for s = 1:subj_num
    
    subject = subjects{s};
    
    record_dir = [subject, '_', drug];
    
    subj_BP_index = strcmp(BP_subjects, subject) & strcmp(BP_drugs, drug);
    
    if ~isempty(states)
       
        subj_state_index = zeros(sum(subj_BP_index), 1);
        
        for state = 1:no_states
            
            subj_state_index = subj_state_index | strcmp(BP_states(subj_BP_index), states{state});
            
        end
        
    else
        
        subj_state_index = ones(sum(subj_BP_index), 1);
        
    end
    
    clear max_vals power
    
    load([record_dir, '_tdr.mat'])
    
    if length(max_vals) > length(subj_state_index)
        
        max_vals((length(subj_state_index) + 1):end, :) = [];
        
        power((length(subj_state_index) + 1):end, :) = [];
        
    end
    
    val_ratio = max_vals(:, 1)./max_vals(:, 2);
    
    pow_ratio = power(:, 1)./power(:, 2);
    
    criteria = [val_ratio pow_ratio power];
    
    subj_BP = BP(subj_BP_index, :, :);
    
    [indices, non_indices] = deal(logical(zeros(size(criteria))));
    
    for c = 1:no_criteria
        
        indices(:, c) = logical(criteria(:, c) < quantile(criteria(subj_state_index, c), quantile_used) & subj_state_index);
        
        length_selected_tBP = sum(indices(:, c));
        
        indices(:, c) = logical(indices(:, c));
        
        theta_BP(tBP_marker(c) + (1:length_selected_tBP), :, :, c) = subj_BP(indices(:, c), :, :);
        
        mean_subj_tBP(:, :, c, s) = squeeze(nanmean(subj_BP(indices(:, c), :, :), 1));
        
        tBP_marker(c) = tBP_marker(c) + length_selected_tBP;
        
        non_indices(:, c) = logical(criteria(:, c) > quantile(criteria(subj_state_index, c), 1 - quantile_used) & subj_state_index);
        
        length_selected_ntBP = sum(non_indices(:, c));
        
        non_indices(:, c) = logical(non_indices(:, c));
        
        non_theta_BP(ntBP_marker(c) + (1:length_selected_ntBP), :, :, c) = subj_BP(non_indices(:, c), :, :);
        
        mean_subj_ntBP(:, :, c, s) = squeeze(nanmean(subj_BP(non_indices(:, c), :, :), 1));
        
        ntBP_marker(c) = ntBP_marker(c) + length_selected_ntBP;
        
    end
    
    for p = 1:no_pairs
        
        index = indices(:, pairs(p, 1)) & indices(:, pairs(p, 2));
        
        length_selected_tBP = sum(index);
        
        index = logical(index);
        
        theta_BP(tBP_marker(no_criteria + p) + (1:length_selected_tBP), :, :, no_criteria + p) = subj_BP(index, :, :);
        
        mean_subj_tBP(:, :, no_criteria + p, s) = squeeze(nanmean(subj_BP(index, :, :), 1));
        
        tBP_marker(no_criteria + p) = tBP_marker(no_criteria + p) + length_selected_tBP;
        
        non_index = non_indices(:, pairs(p, 1)) & non_indices(:, pairs(p, 2));
        
        length_selected_ntBP = sum(non_index);
        
        non_index = logical(non_index);
        
        non_theta_BP(ntBP_marker(no_criteria + p) + (1:length_selected_ntBP), :, :, no_criteria + p) = subj_BP(non_index, :, :);
        
        mean_subj_ntBP(:, :, no_criteria + p, s) = squeeze(nanmean(subj_BP(non_index, :, :), 1));
        
        ntBP_marker(no_criteria + p) = ntBP_marker(no_criteria + p) + length_selected_ntBP;
        
    end
    
    index = cumprod(indices, 2); index = logical(index(:, end));
    
    length_selected_tBP = sum(index);
    
    theta_BP(tBP_marker(no_figures) + (1:length_selected_tBP), :, :, no_figures) = subj_BP(index, :, :);
        
    mean_subj_tBP(:, :, no_figures, s) = squeeze(nanmean(subj_BP(index, :, :), 1));
    
    tBP_marker(no_figures) = tBP_marker(no_figures) + length_selected_tBP;
    
    non_index = cumprod(non_indices, 2); non_index = logical(non_index(:, 3));
    
    length_selected_ntBP = sum(non_index);
    
    non_index = logical(non_index);
    
    non_theta_BP(ntBP_marker(no_figures) + (1:length_selected_ntBP), :, :, no_figures) = subj_BP(non_index, :, :);
        
    mean_subj_ntBP(:, :, no_figures, s) = squeeze(nanmean(subj_BP(non_index, :, :), 1));
    
    ntBP_marker(no_figures) = ntBP_marker(no_figures) + length_selected_ntBP;
    
end

% [mean_tBP, mean_ntBP] = deal(nan(size(BP, 2), size(BP, 3), no_figures));

for f = 1:no_figures

    theta_BP(sum(theta_BP(:, :, f), 3) == 0, :, f) = nan;

    non_theta_BP(sum(non_theta_BP(:, :, f), 3) == 0, :, f) = nan;
    
end
    
mean_tBP = squeeze(nanmean(theta_BP, 1));

mean_ntBP = squeeze(nanmean(non_theta_BP, 1));

save([drug, '_theta_BP_q', num2str(quantile_used), state_label, '_tails.mat'], '-v7.3',...
    'drug', 'state', 'quantile_used', 'theta_BP', 'mean_subj_tBP', 'mean_tBP',...
    'non_theta_BP', 'mean_subj_ntBP', 'mean_ntBP')