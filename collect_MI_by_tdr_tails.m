function collect_MI_by_tdr_tails(channel, drug, quantile_used, states, summed_flag)

load('subjects.mat'), load('AP_freqs.mat')
    
state_label = '';

if ~isempty(states)
    
    no_states = length(states);
    
    for state = 1:no_states
        
        state_label = [state_label, '_', states{state}];
        
    end
    
end

name = ['ALL_', channel];

measure = 'p0.99_IEzs';

if summed_flag
    suffix = 'summed.mat';
    summed_tag = 'summed_';
else
    suffix = 'hr_MI.txt';
    summed_tag = '';
end

measure = 'p0.99_IEzs';

MI_drugs = text_read([name,'/',name,'_',measure,'_drugs.txt'],'%s');
MI_subjects = text_read([name,'/',name,'_',measure,'_subjects.txt'],'%s');
MI_states = text_read([name,'/',name,'_',measure,'_states.txt'],'%s');
MI = load([name, '/', name, '_', measure, '_', suffix]);

if summed_flag
    MI = MI.summed_MI;
end

no_criteria = 4;
no_pairs = nchoosek(no_criteria, 2);
no_figures = no_criteria + no_pairs + 1;

theta_MI = nan(ceil(quantile_used*size(MI, 1)), size(MI, 2), no_figures);
non_theta_MI = nan(ceil((1 - quantile_used)*size(MI, 1)), size(MI, 2), no_figures);

[median_subj_tMI, median_subj_ntMI] = deal(nan(size(MI, 2), no_figures, subj_num));

[tMI_marker, ntMI_marker] = deal(zeros(no_figures, 1));

pairs = nchoosek(1:no_criteria, 2);

theta_labels = {'Peak Power \delta/\theta', 'Band Power \delta/\theta', '\delta Power', '\theta Power'};

for s = 1:subj_num
    
    subject = subjects{s};
    
    record_dir = [subject, '_', drug];
    
    subj_MI_index = strcmp(MI_subjects, subject) & strcmp(MI_drugs, drug);
    
    if ~isempty(states)
       
        subj_state_index = zeros(sum(subj_MI_index), 1);
        
        for state = 1:no_states
            
            subj_state_index = subj_state_index | strcmp(MI_states(subj_MI_index), states{state});
            
        end
        
    else
        
        subj_state_index = ones(sum(subj_MI_index), 1);
        
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
    
    subj_MI = MI(subj_MI_index, :);
    
    [indices, non_indices] = deal(logical(zeros(size(criteria))));
    
    for c = 1:no_criteria
        
        indices(:, c) = logical(criteria(:, c) < quantile(criteria(subj_state_index, c), quantile_used) & subj_state_index);
        
        length_selected_tMI = sum(indices(:, c));
        
        indices(:, c) = logical(indices(:, c));
        
        theta_MI(tMI_marker(c) + (1:length_selected_tMI), :, c) = subj_MI(indices(:, c), :);
        
        median_subj_tMI(:, c, s) = nanmedian(subj_MI(indices(:, c), :))';
        
        tMI_marker(c) = tMI_marker(c) + length_selected_tMI;
        
        non_indices(:, c) = logical(criteria(:, c) > quantile(criteria(subj_state_index, c), 1 - quantile_used) & subj_state_index);
        
        length_selected_ntMI = sum(non_indices(:, c));
        
        non_indices(:, c) = logical(non_indices(:, c));
        
        non_theta_MI(ntMI_marker(c) + (1:length_selected_ntMI), :, c) = subj_MI(non_indices(:, c), :);
        
        median_subj_ntMI(:, c, s) = nanmedian(subj_MI(non_indices(:, c), :))';
        
        ntMI_marker(c) = ntMI_marker(c) + length_selected_ntMI;
        
    end
    
    for p = 1:no_pairs
        
        index = indices(:, pairs(p, 1)) & indices(:, pairs(p, 2));
        
        length_selected_tMI = sum(index);
        
        index = logical(index);
        
        theta_MI(tMI_marker(no_criteria + p) + (1:length_selected_tMI), :, no_criteria + p) = subj_MI(index, :);
        
        median_subj_tMI(:, no_criteria + p, s) = nanmedian(subj_MI(index, :))';
        
        tMI_marker(no_criteria + p) = tMI_marker(no_criteria + p) + length_selected_tMI;
        
        non_index = non_indices(:, pairs(p, 1)) & non_indices(:, pairs(p, 2));
        
        length_selected_ntMI = sum(non_index);
        
        non_index = logical(non_index);
        
        non_theta_MI(ntMI_marker(no_criteria + p) + (1:length_selected_ntMI), :, no_criteria + p) = subj_MI(non_index, :);
        
        median_subj_ntMI(:, no_criteria + p, s) = nanmedian(subj_MI(non_index, :))';
        
        ntMI_marker(no_criteria + p) = ntMI_marker(no_criteria + p) + length_selected_ntMI;
        
    end
    
    index = cumprod(indices, 2); index = logical(index(:, end));
    
    length_selected_tMI = sum(index);
    
    theta_MI(tMI_marker(no_figures) + (1:length_selected_tMI), :, no_figures) = subj_MI(index, :);
        
    median_subj_tMI(:, no_figures, s) = nanmedian(subj_MI(index, :))';
    
    tMI_marker(no_figures) = tMI_marker(no_figures) + length_selected_tMI;
    
    non_index = cumprod(non_indices, 2); non_index = logical(non_index(:, 3));
    
    length_selected_ntMI = sum(non_index);
    
    non_index = logical(non_index);
    
    non_theta_MI(ntMI_marker(no_figures) + (1:length_selected_ntMI), :, no_figures) = subj_MI(non_index, :);
        
    median_subj_ntMI(:, no_figures, s) = nanmedian(subj_MI(non_index, :))';
    
    ntMI_marker(no_figures) = ntMI_marker(no_figures) + length_selected_ntMI;
    
end

[median_tMI, median_ntMI] = deal(nan(size(MI, 2), no_figures));

for f = 1:no_figures

    theta_MI(sum(theta_MI(:, :, f), 2) == 0, :, f) = nan;
    
    median_tMI(:, f) = nanmedian(theta_MI(:, :, f))';

    non_theta_MI(sum(non_theta_MI(:, :, f), 2) == 0, :, f) = nan;
    
    median_ntMI(:, f) = nanmedian(non_theta_MI(:, :, f))';
    
end

save([drug, '_theta_', summed_tag, 'MI_q', num2str(quantile_used), state_label, '_tails.mat'], '-v7.3',...
    'drug', 'state', 'quantile_used', 'theta_MI', 'median_subj_tMI', 'median_tMI',...
    'non_theta_MI', 'median_subj_ntMI', 'median_ntMI')