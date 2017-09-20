function deltaHFO_thetaHFO_quantiles_hi_hr_overlap_plot(q)

load('subjects.mat')

load('drugs.mat')

load('channels.mat')

load(sprintf('dHFO_tHFO_quantiles_hihr_q%.3g.mat', q), 'overlap')

% Dimensions: size(overlap) = [no_subjects, no_drugs, no_channels].

figure

barwitherr(reshape(nanstd(overlap), no_drugs, no_channels), reshape(nanmean(overlap), no_drugs, no_channels))

set(gca, 'XTickLabel', drugs)

title({sprintf('Overlap of %.3g^{th} Percentiles', 100*q);'Highest \delta-HFO and \theta-HFO PAC'})

ylabel('Overlap Proportion')

legend({'Fr.', 'Occi.', 'CA1'})

save_as_pdf(gcf, sprintf('dHFO_tHFO_quantiles_hihr_q%.3g_overlap', q))

figure

barwitherr(nanstd(reshape(overlap, no_subjects*no_drugs, no_channels))', nanmean(reshape(overlap, no_drugs*no_subjects, no_channels))')

set(gca, 'XTickLabel', {'Fr.', 'Occi.', 'CA1'})

title({sprintf('Overlap of %.3g^{th} Percentiles', 100*q);'Highest \delta-HFO and \theta-HFO PAC'})

ylabel('Overlap Proportion')

save_as_pdf(gcf, sprintf('dHFO_tHFO_quantiles_hihr_q%.3g_overlap_all_drugs', q))