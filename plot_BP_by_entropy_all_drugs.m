function plot_BP_by_entropy_all_drugs(quantile_used, states)

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

entropy_index = 4;

band_labels = {'Frontal \delta Power', 'CA1 \theta Power', '\delta-\theta PAC'}; no_bands = length(band_labels);
y_labels = {'%\Delta from Baseline'; '%\Delta from Baseline'; 'z-Scored MI'};

power = nan(no_subjects, 2, no_bands, no_drugs);

for d = 1:no_drugs
    
    figure
    
    drug = drugs{d};
    
    load([drug, '_delta_BP_q', num2str(quantile_used), shm_label, state_label, '_tails.mat'])
    
    for b = 1:(no_bands - 1)
        
        temp_power = [squeeze(mean_subj_dBP(b, b, entropy_index, :)) squeeze(mean_subj_ndBP(b, b, entropy_index, :))];
        
        power(:, :, d, b) = temp_power;
        
        subplot(1, no_bands, b)
        
        barwitherr(nanstd(temp_power), nanmean(temp_power))
        
        box off
        
        title(band_labels{b}, 'FontSize', 16)
        
        set(gca, 'XTickLabel', {'Narrowband', 'Broadband'})
        
        ylabel(y_labels{b})
        
    end
    
    for b = no_bands
        
        load(['Occipital_', drug, '_delta_summed_MI_q', num2str(quantile_used), state_label, '_tails.mat'])
        
        temp_power = [squeeze(median_subj_dMI(5, entropy_index, :)) squeeze(median_subj_ndMI(5, entropy_index, :))];
        
        power(:, :, d, b) = temp_power;
        
        subplot(1, no_bands, b)
        
        barwitherr(nanstd(temp_power), nanmean(temp_power))
        
        box off
        
        title(band_labels{b}, 'FontSize', 16)
        
        set(gca, 'XTickLabel', {'Narrowband', 'Broadband'})
        
        ylabel(y_labels{b})
        
    end
    
    mtit([drug, ', BP During Narrowband Delta', long_state_label])
    
    save_as_pdf(gcf, [drug, '_delta_BP_q', num2str(quantile_used), state_label, '_power_comparison'])
    
end

figure

for b = 1:no_bands
    
    subplot(1, no_bands, b)
    
    barwitherr(squeeze(nanstd(power(:, :, :, b))), squeeze(nanmean(power(:, :, :, b))))
    
    box off
    
    legend(drugs)
    
    title(band_labels{b}, 'FontSize', 16)
    
    set(gca, 'XTickLabel', {'Narrowband', 'Broadband'})
    
    ylabel(y_labels{b})
    
end

mtit(['BP During Narrowband Delta', long_state_label])

save_as_pdf(gcf, ['delta_BP_q', num2str(quantile_used), state_label, '_power_comparison'])

power = permute(power, [1 3 2 4]);

figure

p = nan(no_bands, 1);

for b = 1:no_bands
    
    subplot(1, no_bands, b)
    
    power_for_test = reshape(power(:, :, :, b), [no_subjects*no_drugs, 2]);
    
    [~, p(b)] = ttest(power_for_test(:, 1), power_for_test(:, 2));
    
    handle = barwitherr(squeeze(nanstd(power_for_test))/sqrt(no_subjects*no_drugs),...
        squeeze(nanmean(power_for_test)));
        
    bar_pos = get_bar_pos(handle);
    
    if p(b) < .05, sigstar({bar_pos}, p(b)), end
    
    box off
    
    title(band_labels{b}, 'FontSize', 16)
    
    set(gca, 'XTickLabel', {'Narrowband', 'Broadband'})
    
    ylabel(y_labels{b})
    
end

mtit(['BP & PAC During Narrowband Delta', long_state_label])

save_as_pdf(gcf, ['delta_BP_q', num2str(quantile_used), state_label, '_power_comparison_over_drugs'])

end