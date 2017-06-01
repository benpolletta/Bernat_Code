function collect_MI_by_whm_tails(drug, quantile_used, states)

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

MI_drugs = text_read([name,'/',name,'_',measure,'_drugs.txt'],'%s');
MI_subjects = text_read([name,'/',name,'_',measure,'_subjects.txt'],'%s');
MI_states = text_read([name,'/',name,'_',measure,'_states.txt'],'%s');
MI = load([name, '/', name, '_', measure, '_hr_MI.txt']);


criteria = {'WHM', 'SHM', 'SHM/WHM', 'Entropy', 'Peak Freq.', 'Peak Pow.', 'Pow.'};
no_criteria = length(criteria);
no_pairs = nchoosek(no_criteria, 2);
no_figures = no_criteria + no_pairs + 1;

delta_MI = nan(ceil(quantile_used*size(MI, 1)), size(MI, 2), no_figures);
non_delta_MI = nan(ceil((1 - quantile_used)*size(MI, 1)), size(MI, 2), no_figures);

[median_subj_dMI, median_subj_ndMI] = deal(nan(size(MI, 2), no_figures, subj_num));

[dMI_marker, ndMI_marker] = deal(zeros(no_figures, 1));

pairs = nchoosek(1:no_criteria, 2);

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

    clear criteria

    whm_struct = load([record_dir, '_chan1_whm.mat']);
    whm_fields = fieldnames(whm_struct);
    
    for f = 1:length(whm_fields), eval(sprintf('%s = whm_struct.%s;', whm_fields{f}, whm_fields{f})), end
    
    entropy = load([record_dir, '_chan1_whm.mat'], 'entropy');
    entropy = entropy.entropy;
    
    power = load([record_dir, '_chan1_whm.mat'], 'power');
    power = power.power;
    
    criteria = [whm shm_sum shm_sum./whm entropy max_freqs max_vals power];

    if length(whm) > length(subj_state_index)

        criteria((length(subj_state_index) + 1):end, :) = [];

    end
    
    [delta_indices, non_delta_indices] = deal(nan(size(criteria)));

    for c = 1:size(criteria, 2)
       
        if c <= 5
            
            delta_indices(:, c) = criteria(:, c) < quantile(criteria(subj_state_index, c), quantile_used)...
                & subj_state_index;
            
            non_delta_indices(:, c) = criteria(:, c) > quantile(criteria(subj_state_index, c), 1 - quantile_used)...
                & subj_state_index;
            
        else
            
            delta_indices(:, c) = criteria(:, c) > quantile(criteria(subj_state_index, c), 1 - quantile_used)...
                & subj_state_index;
            
            non_delta_indices(:, c) = criteria(:, c) > quantile(criteria(subj_state_index, c), quantile_used)...
                & subj_state_index;
            
        end
        
    end
    
    delta_indices = logical(delta_indices); non_delta_indices = logical(non_delta_indices);

    subj_MI = MI(subj_MI_index, :);

    for c = 1:no_criteria

        length_selected_dMI = sum(delta_indices(:, c));

        delta_MI(dMI_marker(c) + (1:length_selected_dMI), :, c) = subj_MI(delta_indices(:, c), :);

        median_subj_dMI(:, c, s) = nanmedian(subj_MI(delta_indices(:, c), :))';

        dMI_marker(c) = dMI_marker(c) + length_selected_dMI;

        length_selected_ndMI = sum(non_delta_indices(:, c));

        non_delta_MI(ndMI_marker(c) + (1:length_selected_ndMI), :, c) = subj_MI(non_delta_indices(:, c), :);

        median_subj_ndMI(:, c, s) = nanmedian(subj_MI(non_delta_indices(:, c), :))';

        ndMI_marker(c) = ndMI_marker(c) + length_selected_ndMI;

    end

    for p = 1:no_pairs

        index = delta_indices(:, pairs(p, 1)) & delta_indices(:, pairs(p, 2));

        length_selected_dMI = sum(index);

        delta_MI(dMI_marker(no_criteria + p) + (1:length_selected_dMI), :, no_criteria + p) = subj_MI(index, :);

        median_subj_dMI(:, no_criteria + p, s) = nanmedian(subj_MI(index, :))';

        dMI_marker(no_criteria + p) = dMI_marker(no_criteria + p) + length_selected_dMI;

        non_index = non_delta_indices(:, pairs(p, 1)) & non_delta_indices(:, pairs(p, 2));

        length_selected_ndMI = sum(non_index);

        non_delta_MI(ndMI_marker(no_criteria + p) + (1:length_selected_ndMI), :, no_criteria + p) = subj_MI(non_index, :);

        median_subj_ndMI(:, no_criteria + p, s) = nanmedian(subj_MI(non_index, :))';

        ndMI_marker(no_criteria + p) = ndMI_marker(no_criteria + p) + length_selected_ndMI;

    end

    index = cumprod(delta_indices, 2); index = logical(index(:, end));

    length_selected_dMI = sum(index);

    delta_MI(dMI_marker(no_figures) + (1:length_selected_dMI), :, no_figures) = subj_MI(index, :);

    median_subj_dMI(:, no_figures, s) = nanmedian(subj_MI(index, :))';

    dMI_marker(no_figures) = dMI_marker(no_figures) + length_selected_dMI;

    non_index = cumprod(non_delta_indices, 2); non_index = logical(non_index(:, 3));

    length_selected_ndMI = sum(non_index);

    non_delta_MI(ndMI_marker(no_figures) + (1:length_selected_ndMI), :, no_figures) = subj_MI(non_index, :);

    median_subj_ndMI(:, no_figures, s) = nanmedian(subj_MI(non_index, :))';

    ndMI_marker(no_figures) = ndMI_marker(no_figures) + length_selected_ndMI;

end

% end_dMI_marker = max(dMI_marker)
%
% end_ndMI_marker = max(ndMI_marker)
%
% delta_MI((end_dMI_marker + 1):end, :, :) = [];
%
% non_delta_MI((end_ndMI_marker + 1):end, :, :) = [];

[median_dMI, median_ndMI] = deal(nan(size(MI, 2), no_figures));

for f = 1:no_figures

    delta_MI(sum(delta_MI(:, :, f), 2) == 0, :, f) = nan;

    median_dMI(:, f) = nanmedian(delta_MI(:, :, f))';

    non_delta_MI(sum(non_delta_MI(:, :, f), 2) == 0, :, f) = nan;

    median_ndMI(:, f) = nanmedian(non_delta_MI(:, :, f))';

end

save([drug, '_delta_MI_q', num2str(quantile_used), state_label, '_tails.mat'], '-v7.3',...
    'drug', 'state', 'quantile_used', 'delta_MI', 'median_subj_dMI', 'median_dMI',...
    'non_delta_MI', 'median_subj_ndMI', 'median_ndMI')
