function plot_MI_by_whm_non_delta(drug, quantile_used, states, matched_flag)

load('subjects.mat'), load('AP_freqs.mat')

state_label = ''; long_state_label = '';

if ~isempty(states)
    
    no_states = length(states);
    
    for s = 1:no_states
        
        state_label = [state_label, '_', states{s}];
        
        long_state_label = [long_state_label, ', ', states{s}];
        
    end
    
end

load([drug, '_delta_MI_q', num2str(quantile_used), state_label, matched_flag, '.mat'])

criteria = {'WHM', 'SHM', 'SHM/WHM', 'Entropy', 'Peak Freq.', 'Peak Pow.', 'Pow.'};

no_criteria = length(criteria);

pairs = nchoosek(1:no_criteria, 2);

no_pairs = size(pairs, 1);

delta_labels = cell(no_criteria + no_pairs + 1, 1);

for c = 1:5 % no_criteria
    
    delta_labels{c} = sprintf('Highest %g \\delta %s', quantile_used, criteria{c});
    
end

for c = 5:no_criteria
    
    delta_labels{c} = sprintf('Lowest %g \\delta %s', quantile_used, criteria{c});
    
end

for p = 1:no_pairs
    
    delta_labels{no_criteria + p} = delta_labels(pairs(p, :));
    
end

delta_labels{end} = 'Intersection of All Criteria';

no_deltas = length(delta_labels);

[no_rows, no_cols] = subplot_size(no_deltas);

for s = 1:subj_num
    
    subject = subjects{s};
    
    figure
    
    for d = 1:no_deltas
        
        subplot(no_rows, no_cols, d)
        
        imagesc(phase_freqs, amp_freqs, reshape(median_subj_ndMI(:, d, s), no_afs, no_pfs))
        
        axis xy
        
        title(delta_labels{d})
        
    end
    
    mtit([subject, ' MI Outside of Narrowband Delta', long_state_label])
    
    save_as_pdf(gcf, [subject, '_', drug, '_non_delta_MI_q', num2str(quantile_used), shm_label, state_label, matched_flag])
    
end

figure

for d = 1:no_deltas
    
    subplot(no_rows, no_cols, d)
    
    imagesc(phase_freqs, amp_freqs, reshape(median_ndMI(:, d), no_afs, no_pfs))
    
    axis xy
    
    colorbar
    
    title(delta_labels{d})
    
end

mtit(['MI Outside of Narrowband Delta', long_state_label])

save_as_pdf(gcf, [drug, '_non_delta_MI_q', num2str(quantile_used), state_label, matched_flag])