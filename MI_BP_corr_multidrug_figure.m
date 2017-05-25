function MI_BP_corr_multidrug_figure(suffix)

figure

load('AP_freqs'), load('subjects'), load('drugs'), load('channels'), load('BP_bands')

%% MI and BP correlations.

pairs = [cumsum(ones(3, 2)); nchoosek(1:3, 2); fliplr(nchoosek(1:3, 2))];

no_pairs = length(pairs);

channel_indices = [1 3];

no_channels = 3;
        
for ch = 1:no_channels
    
    figure
    
    for d = 1:no_drugs
        
        drug = drugs{d};
        
        load([drug, '_MI_BP_pct_0to4hrs', suffix, '.mat'])
        
        for band = 1:2
            
            pair_index = (pairs(:, 1) == ch) & (pairs(:, 2) == channel_indices(band));
            
            if strcmp(suffix, '_6min_by_subject')
                
                pair_corrs = nanmedian(All_corrs(:, band, pair_index, 1:(end - 1)), 4);
                
            else
                
                pair_corrs = All_corrs(:, band, pair_index);
                
            end
            
            subplot(no_drugs, 2, (d - 1)*2 + band)
            
            imagesc(phase_freqs, amp_freqs, reshape(pair_corrs, no_afs, no_pfs))
            
            colorbar
            
            axis xy
            
            if band == 1
                
                ylabel(drug, 'FontSize', 16)
                
            end
            
            if d == 1
                
                title({['Correlation of ', channel_names{ch}]; ['MI & ', channel_names{channel_indices(band)}, ' ', band_labels_long{band}]}, 'FontSize', 16)
                
            elseif d == no_drugs
                
                xlabel('Frequency (Hz)', 'FontSize', 16)
                
            end
            
        end
        
    end

    save_as_pdf(gcf, [channel_names{ch}, 'MI_BP_corr', suffix])

end

end