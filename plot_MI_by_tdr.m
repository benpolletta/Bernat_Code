function plot_MI_by_tdr(drug, quantile_used, states)

load('subjects.mat'), load('AP_freqs.mat')

state_label = ''; long_state_label = '';

if ~isempty(states)
    
    no_states = length(states);
    
    for s = 1:no_states
        
        state_label = [state_label, '_', states{s}];
        
        long_state_label = [long_state_label, ', ', states{s}];
        
    end
    
end

load([drug, '_theta_MI_q', num2str(quantile_used), state_label, '_tails.mat'])

criteria = {'\delta/\theta (peak)', '\delta/\theta', 'delta', 'theta'};

pairs = nchoosek(1:length(criteria), 2);

for c = 1:length(criteria)
    
    theta_labels{c} = ['Low ', criteria{c}];
    
end

for p = 1:length(pairs)
    
    theta_labels{length(criteria) + p} = ['Low ', criteria{pairs(p, 1)}, ' & ', criteria{pairs(p, 2)}];
    
end

theta_labels{end + 1} = 'Low \theta/\delta (peak) & \theta/\delta & delta & theta';

no_thetas = length(theta_labels);

[no_rows, no_cols] = subplot_size(no_thetas);

for s = 1:subj_num
    
    subject = subjects{s};
    
    figure
    
    for d = 1:length(theta_labels)
        
        subplot(no_rows, no_cols, d)
        
        imagesc(phase_freqs, amp_freqs, reshape(median_subj_tMI(:, d, s), no_afs, no_pfs))
        
        axis xy
        
        title(theta_labels{d})
        
    end
    
    mtit([subject, ' MI by \theta & \delta Power', long_state_label])
    
    save_as_pdf(gcf, [subject, '_', drug, '_theta_MI_q', num2str(quantile_used), state_label])
    
end

figure

for d = 1:no_thetas
    
    subplot(no_rows, no_cols, d)
    
    imagesc(phase_freqs, amp_freqs, reshape(median_tMI(:, d), no_afs, no_pfs))
    
    axis xy
    
    colorbar
    
    title(theta_labels{d})
    
end

mtit(['MI During Narrowband theta', long_state_label])

save_as_pdf(gcf, [drug, '_theta_MI_q', num2str(quantile_used), state_label])