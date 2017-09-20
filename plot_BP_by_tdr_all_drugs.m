function plot_BP_by_tdr_all_drugs(quantile_used, states, index)

load('drugs.mat'), load('subjects.mat'), load('AP_freqs.mat')

state_label = ''; long_state_label = '';

theta_labels = get_theta_labels;

if ~isempty(states)
    
    no_states = length(states);
    
    for s = 1:no_states
        
        state_label = [state_label, '_', states{s}];
        
        long_state_label = [long_state_label, ', ', states{s}];
        
    end
    
end

band_labels = {'Frontal \delta Power', 'CA1 \theta Power'}; no_bands = length(band_labels);

conditions = {'low', 'high'}; no_conditions = length(conditions);

power = nan(no_subjects, no_conditions, no_drugs, no_bands);

[drug_groups, high_low_groups] = deal(cell(no_subjects*no_drugs*no_bands, 1));

for d = 1:no_drugs
    
    figure
    
    drug = drugs{d};
    
    load([drug, '_theta_BP_q', num2str(quantile_used), state_label, '_tails.mat'])
    
    drug_groups((d - 1)*no_subjects + (1:no_subjects)) = {drug};
    drug_groups(no_drugs*no_subjects + (d - 1)*no_subjects + (1:no_subjects)) = {drug};
    
    high_low_groups((d - 1)*no_subjects + (1:no_subjects)) = {'low'};
    high_low_groups(no_drugs*no_subjects + (d - 1)*no_subjects + (1:no_subjects)) = {'high'};
    
    for b = 1:no_bands
        
        temp_power = [squeeze(mean_subj_tBP(b, b, index, :)) squeeze(mean_subj_ntBP(b, b, index, :))];
        
        power(:, :, d, b) = temp_power;
        
        subplot(1, no_bands, b)
        
        % boxplot(temp_power) 
        barwitherr(nanstd(temp_power), nanmean(temp_power))
        
        box off
        
        title(band_labels{b}, 'FontSize', 16)
        
        set(gca, 'XTickLabel', {['Low ', theta_labels{index}], ['High ', theta_labels{index}]})
        
        ylabel('%\Delta from Baseline')
        
    end
    
    mtit([drug, ', BP During ', theta_labels{index}, long_state_label])
    
    save_as_pdf(gcf, [drug, '_theta_BP_q', num2str(quantile_used), state_label, '_', num2str(index), '_comparison'])
    
end

% power = permute(power, [1 3 2 4]);

figure

for b = 1:no_bands
    
    subplot(1, no_bands, b)
    
    % boxplot(reshape(power(:, :, :, b), [size(power, 1)*size(power, 2)*size(power,3), 1]),...
        % {high_low_groups, drug_groups}, 'colorgroup', drug_groups) 
    barwitherr(squeeze(nanstd(power(:, :, :, b))), squeeze(nanmean(power(:, :, :, b))))
    
    box off
    
    legend(drugs)
    
    title(band_labels{b}, 'FontSize', 16)
    
    set(gca, 'XTickLabel', {['Low ', theta_labels{index}], ['High ', theta_labels{index}]})
    
    ylabel('%\Delta from Baseline')
    
end

mtit(['BP During ', theta_labels{index}, long_state_label])

save_as_pdf(gcf, ['theta_BP_q', num2str(quantile_used), state_label, '_', num2str(index), '_comparison'])

power = permute(power, [1 3 2 4]);

figure

for b = 1:no_bands
    
    subplot(1, no_bands, b)
    
    % boxplot(reshape(power(:, :, :, b), [size(power, 1)*size(power, 2)*size(power,3), 1]),...
        % high_low_groups, 'colorgroup', high_low_groups) 
    barwitherr(squeeze(nanstd(reshape(power(:, :, :, b), [no_subjects*no_drugs, 2]))),...
        squeeze(nanmean(reshape(power(:, :, :, b), [no_subjects*no_drugs, 2]))))
    
    box off
    
    title(band_labels{b}, 'FontSize', 16)
    
    set(gca, 'XTickLabel', {['Low ', theta_labels{index}], ['High ', theta_labels{index}]})
    
    ylabel('%\Delta from Baseline')
    
end

mtit(['BP During ', theta_labels{index}, long_state_label])

save_as_pdf(gcf, ['theta_BP_q', num2str(quantile_used), state_label, '_', num2str(index), '_comparison_over_drugs'])

end


function theta_labels = get_theta_labels

criteria = {'\delta/\theta (peak)', '\delta/\theta', 'delta', 'theta'};

pairs = nchoosek(1:length(criteria), 2);

theta_labels = criteria;

for p = 1:length(pairs)
    
    theta_labels{length(criteria) + p} = [criteria{pairs(p, 1)}, ' & ', criteria{pairs(p, 2)}];
    
end

theta_labels{end + 1} = '\theta/\delta (peak) & \theta/\delta & delta & theta';

end