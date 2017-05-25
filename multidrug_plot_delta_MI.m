function multidrug_plot_delta_MI(quantile_used, shm_lims, states)

load('AP_freqs')

load('drugs')

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

suffix = ['_delta_MI_q', num2str(quantile_used), '_shm', sprintf('%.03f_%.03f', shm_lims), state_label, '_tails'];

figure

entropy_index = 4;

for d = 1:no_drugs
    
    delta_MI_name = [drugs{d}, suffix, '.mat'];
    
    load(delta_MI_name)
    
    cmin = min(min(median_dMI(:, entropy_index)), min(median_ndMI(:, entropy_index)));
    
    cmax = max(max(median_dMI(:, entropy_index)), max(median_ndMI(:, entropy_index)));
    
    subplot(no_drugs, 2, 2*(d - 1) + 1)
    
    imagesc(phase_freqs, amp_freqs, reshape(median_dMI(:, entropy_index), no_afs, no_pfs))
    
    axis xy
    
    caxis([cmin cmax])
    
    ylabel(drugs{d})
    
    if d == 1
    
        title({'Narrowband Delta'; long_state_label}, 'FontSize', 16)
        
    end
    
    subplot(no_drugs, 2, 2*(d - 1) + 2)
    
    imagesc(phase_freqs, amp_freqs, reshape(median_ndMI(:, entropy_index), no_afs, no_pfs))
    
    axis xy
    
    caxis([cmin cmax])
    
    if d == 1
    
        title({'Broadband Delta'; long_state_label}, 'FontSize', 16)
    
    end
    
end

save_as_pdf(gcf, ['multidrug', suffix])

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