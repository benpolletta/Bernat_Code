function multidrug_plot_delta_MI_w_power(quantile_used, shm_lims, states)

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


suffix = ['_delta_MI_q', num2str(quantile_used), make_label('shm', shm_lims, []), state_label, '_tails'];

figure

entropy_index = 4;

for d = 1:no_drugs
    
    delta_MI_name = [drugs{d}, suffix, '.mat'];
    
    load(delta_MI_name)
    
    cmin = min(min(median_dMI(:, entropy_index)), min(median_ndMI(:, entropy_index)));
    
    cmax = max(max(median_dMI(:, entropy_index)), max(median_ndMI(:, entropy_index)));
    
    subplot(no_drugs + 1, 2, 2*(d - 1) + 1)
    
    imagesc(phase_freqs, amp_freqs, reshape(median_dMI(:, entropy_index), no_afs, no_pfs))
    
    set(gca, 'FontSize', 16)
    
    axis xy
    
    caxis([cmin cmax])
    
    ylabel(drugs{d})
    
    if d == 1
    
        title({'Narrowband Delta'; long_state_label}, 'FontSize', 16)
        
    end
    
    subplot(no_drugs + 1, 2, 2*(d - 1) + 2)
    
    imagesc(phase_freqs, amp_freqs, reshape(median_ndMI(:, entropy_index), no_afs, no_pfs))
    
    set(gca, 'FontSize', 16)
    
    axis xy
    
    caxis([cmin cmax])
    
    if d == 1
    
        title({'Broadband Delta'; long_state_label}, 'FontSize', 16)
    
    end
    
end

shm_label = make_label('shm', [1.025 1.25], []);

load('drugs.mat'), load('subjects.mat'), load('AP_freqs.mat')

state_label = ''; long_state_label = '';

if ~isempty(states)
    
    no_states = length(states);
    
    for s = 1:no_states
        
        state_label = [state_label, '_', states{s}];
        
        long_state_label = [long_state_label, ', ', states{s}];
        
    end
    
end

band_labels = {'Frontal \delta', 'CA1 \theta'}; no_bands = length(band_labels);

power = nan(no_subjects, no_bands, no_drugs, 2);

for b = 1:no_bands
    
    for d = 1:no_drugs
        
        drug = drugs{d};
        
        load([drug, '_delta_BP_q', num2str(quantile_used),  make_label('shm', [1.025 1.25], []), state_label, '_tails.mat'])
        
        temp_power(:, 1, 1, 1) = squeeze(mean_subj_dBP(b, b, entropy_index, :));
        temp_power(:, 1, 1, 2) = squeeze(mean_subj_ndBP(b, b, entropy_index, :));
        
        power(:, b, d, :) = temp_power;
        
    end
    
end

for low_high = 1:2
    
    subplot(no_drugs + 1, 2, no_drugs*2 + low_high)
    
    barwitherr(squeeze(nanstd(power(:, :, :, low_high))), squeeze(nanmean(power(:, :, :, low_high))))
    
    box off
    
    if low_high == 1, legend(drugs, 'FontSize', 10), end
    
    set(gca, 'XTickLabel', band_labels, 'FontSize', 16)
    
    ylabel('Power (%\Delta BL)')
    
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