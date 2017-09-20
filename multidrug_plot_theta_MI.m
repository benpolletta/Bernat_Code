function multidrug_plot_theta_MI(quantile_used, states, measure_index)

load('AP_freqs')

load('drugs')

load('subjects')

state_label = '';

if ~isempty(states)
    
    no_states = length(states);
    
    for s = 1:no_states
        
        state_label = [state_label, '_', states{s}];
        
        if s == 1
            
            long_state_label = get_long_state(states{s});
            
        else
        
            long_state_label = [long_state_label, ', ', get_long_state(states{s})];
            
        end
        
    end
    
end

suffix = ['_theta_MI_q', num2str(quantile_used), state_label, '_tails'];

theta_labels = make_theta_labels;

figure

for d = 1:no_drugs
    
    theta_MI_name = [drugs{d}, suffix, '.mat'];
    
    load(theta_MI_name)
    
    cmin = min(min(median_tMI(:, measure_index)), min(median_ntMI(:, measure_index)));
    
    cmax = max(max(median_tMI(:, measure_index)), max(median_ntMI(:, measure_index)));
    
    subplot(no_drugs, 3, 3*(d - 1) + 1)
    
    imagesc(phase_freqs, amp_freqs, reshape(median_tMI(:, measure_index), no_afs, no_pfs))
    
    set(gca, 'FontSize', 16)
    
    axis xy
    
    caxis([cmin cmax])
    
    ylabel({drugs{d}; 'Amp. Freq. (Hz)'})
    
    if d == 1
    
        if measure_index <= 2
            
            title({['Highest ', num2str(quantile_used), ' ', theta_labels{measure_index}]; long_state_label}, 'FontSize', 16)
            
        else
    
            title({['Lowest ', num2str(quantile_used), ' ', theta_labels{measure_index}]; long_state_label}, 'FontSize', 16)
            
        end
        
    elseif d == 4
        
        xlabel('Phase Freq. (Hz)')
        
    end
    
    subplot(no_drugs, 3, 3*(d - 1) + 2)
    
    imagesc(phase_freqs, amp_freqs, reshape(median_ntMI(:, measure_index), no_afs, no_pfs))
    
    set(gca, 'FontSize', 16)
    
    axis xy
    
    caxis([cmin cmax])
    
    if d == 1
    
        if measure_index <= 2
            
            title({['Lowest ', num2str(quantile_used), ' ', theta_labels{measure_index}]; long_state_label}, 'FontSize', 16)
            
        else
            
            title({['Highest ', num2str(quantile_used), ' ', theta_labels{measure_index}]; long_state_label}, 'FontSize', 16)
            
        end
        
    elseif d == 4
        
        xlabel('Phase Freq. (Hz)')
    
    end
    
end

save_as_pdf(gcf, ['multidrug', suffix, '_measure', num2str(measure_index)])

end


function long_state = get_long_state(state)

switch state
    
    case 'W'
        long_state = 'Wake';
    case 'NR'
        long_state = 'QW & nREM';
    case 'R'
        long_state = 'REM';
        
end

end

function theta_labels = make_theta_labels

theta_labels = {'\theta/\delta (peak)', '\theta/\delta', '\delta', '\theta'};

alt_theta_labels = {'\delta/\theta (peak)', '\delta/\theta', '\delta', '\theta'};

pairs = nchoosek(1:length(theta_labels), 2);

for p = 1:length(pairs)
    
    theta_labels{length(theta_labels) + p} = [alt_theta_labels{pairs(p, 1)}, ' & ', alt_theta_labels{pairs(p, 2)}];
    
end

theta_labels{end + 1} = '\delta/\theta (peak) & \delta/\theta & \delta & \theta';

end