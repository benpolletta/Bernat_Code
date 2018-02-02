function deltaHFO_thetaHFO_quantiles_by_state(q, state, hours)

HFO_labels = {'\delta-HFO PAC', '\theta-HFO PAC'};

load('subjects.mat'), subject_labels = subjects; clear subjects

load('channels.mat')

load('drugs.mat'), drug_labels = drugs; clear drugs

load('AP_freqs.mat')

measure = 'p0.99_IEzs';
    
median_qPAC = nan(no_afs*no_pfs, 2, no_drugs, no_channels);

overlap = nan(no_subjects, no_drugs, no_channels);

name = sprintf('dHFO_tHFO_quantiles_post%dto%d_q%.3g_%s.mat', min(hours), max(hours), q, state);

for ch = 1:no_channels
    
    channel = channel_names{ch};
    
    ch_dir = ['ALL_', channel];
    
    sMI = load([ch_dir, '/', ch_dir, '_', measure, '_summed.mat']);
    sMI = sMI.summed_MI;
    drugs = text_read([ch_dir,'/',ch_dir, '_', measure,'_drugs.txt'],'%s');
    states = text_read([ch_dir, '/', ch_dir, '_', measure, '_states.txt'], '%s');
    subjects = text_read([ch_dir,'/',ch_dir, '_', measure,'_subjects.txt'],'%s');
    hrs = text_read([ch_dir,'/',ch_dir, '_', measure,'_hrs.txt'],'%s');
    fourhrs = text_read([ch_dir, '/', ch_dir, '_', measure,'_4hrs.txt'],'%s');
    MI = load([ch_dir, '/', ch_dir, '_', measure, '_hr_MI.txt']);
    
    clear index
    
    index = false(size(hrs));
    
    for h = 1:length(hours)
        
        index = index | strcmp(hrs, sprintf('post%d', hours(h)));
        
    end
    
    index = index & strcmp(states, state);
    
    PAC_per_subj = ceil(q*(60*60/(4.096*4)));
    
    selected_PAC = nan(no_subjects*PAC_per_subj, no_afs*no_pfs, 2, no_drugs);
    
    for s = 1:no_subjects
            
        subject = subject_labels{s};
        
        clear subj_index
        
        subj_index = strcmp(subjects, subject) & index;
        
        figure
        
        for d = 1:no_drugs
            
            drug = drug_labels{d};
            
            clear drug_index drug_PAC_index
            
            drug_index = subj_index & strcmp(drugs, drug);
                
            drug_PAC_index = repmat(drug_index, 1, 2);
            
            % Getting high deltaHFO & thetaHFO comodulograms.
            
            for lf = 1:2
                
                high_PAC_index = sMI(drug_index, end - 2 + lf) >=...
                    quantile(sMI(drug_index, end - 2 + lf), 1 - q);
                
                drug_PAC_index(drug_index, lf) = high_PAC_index;
                
                selected_PAC((s - 1)*PAC_per_subj + (1:sum(high_PAC_index)), :, lf, d) = MI(drug_PAC_index(:, lf), :);
            
                subplot(no_drugs, 2, (d - 1)*2 + lf)
                
                imagesc(phase_freqs, amp_freqs, reshape(nanmedian(MI(drug_PAC_index(:, lf), :)), no_afs, no_pfs))
                
                axis xy
                
                colorbar
                
                if lf == 1, ylabel({drug; 'Amp. Freq. (Hz)'}, 'FontSize', 14), end
                
                if d == 1
                    
                    if lf == 1
                        
                        title({[state, ', ', subject, ' ', channel]; HFO_labels{lf}}, 'FontSize', 16)
                    
                    else
                        
                        title(HFO_labels{lf}, 'FontSize', 16)
                    
                    end
                    
                elseif d == no_drugs
                    
                    xlabel('Phase Freq. (Hz)', 'FontSize', 14)
                
                end
            
            end
            
            overlap(s, d, ch) = sum(sum(drug_PAC_index, 2) > 1)/mean(sum(drug_PAC_index));
            
        end
        
        save_as_pdf(gcf, sprintf('%s_%s_%s', subject, channel, name))
        
    end
    
    figure
    
    for d = 1:no_drugs
        
        drug = drug_labels{d};
        
        for lf = 1:2
            
            subplot(no_drugs, 2, (d - 1)*2 + lf)
                
            median_qPAC(:, lf, d, ch) = reshape(nanmedian(selected_PAC(:, :, lf, d)), no_afs*no_pfs, 1);
            
            imagesc(phase_freqs, amp_freqs, reshape(nanmedian(selected_PAC(:, :, lf, d)), no_afs, no_pfs))
            
            axis xy
            
            colorbar
            
            if lf == 1
                
                ylabel({drug; 'Amp. Freq. (Hz)'}, 'FontSize', 14)
            
            end
            
            if d == 1
                
                if lf == 1
                    
                    title({[state, ', ', channel]; HFO_labels{lf}}, 'FontSize', 16)
                
                else
                    
                    title(HFO_labels{lf}, 'FontSize', 16)
                
                end
                
            elseif d == no_drugs
                
                xlabel('Phase Freq. (Hz)', 'FontSize', 14)
                
            end
            
        end
        
    end
    
    save_as_pdf(gcf, sprintf('%s_%s', channel, name))
    
end

save(name, 'median_qPAC', 'overlap')

figure

barwitherr(reshape(nanstd(overlap), no_drugs, no_channels), reshape(nanmean(overlap), no_drugs, no_channels))

set(gca, 'XTickLabel', drug_labels)

title({state; sprintf('Overlap of %.3g^{th} Percentiles', 100*q);'Highest \delta-HFO and \theta-HFO PAC'})

ylabel('Overlap Proportion')

legend({'Fr.', 'Occi.', 'CA1'})

save_as_pdf(gcf, sprintf('%s_overlap', name))