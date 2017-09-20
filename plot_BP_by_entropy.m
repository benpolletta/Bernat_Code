function plot_BP_by_entropy(drug, quantile_used, states)

shm_label = make_label('shm', [1.025 1.25], []);

load('subjects.mat'), load('AP_freqs.mat')

state_label = ''; long_state_label = '';

if ~isempty(states)
    
    no_states = length(states);
    
    for s = 1:no_states
        
        state_label = [state_label, '_', states{s}];
        
        long_state_label = [long_state_label, ', ', states{s}];
        
    end
    
end

load([drug, '_delta_BP_q', num2str(quantile_used), shm_label, state_label, '_tails.mat'])

entropy_index = 4;

figure

band_labels = {'\delta', '\theta'};

for b = 1:2
    
    power = [squeeze(mean_subj_dBP(b, entropy_index, :)) squeeze(mean_subj_ndBP(b, entropy_index, :))];
    
    subplot(1, 2, b)
    
    barwitherr(nanstd(power), nanmean(power))
    
    title(band_labels{b}, 'FontSize', 16)
    
    set(gca, 'XTickLabel', {'Narrowband', 'Broadband'})
    
    ylabel('%\Delta from Baseline')
    
end

mtit(['BP During Narrowband Delta', long_state_label])

save_as_pdf(gcf, [drug, '_delta_BP_q', num2str(quantile_used), state_label, '_power_comparison'])