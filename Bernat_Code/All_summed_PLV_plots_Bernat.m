function All_summed_PLV_plots_Bernat(channel_label1, channel_label2, stat)

state_labels={'W','NR','R'};

drug_labels={'saline','MK801','NVP','Ro25'};

name = sprintf('ALL_%s_by_%s_PLV', channel_label1, channel_label2);
title_name = sprintf('%s by %s, ', channel_label1, channel_label2);

drugs = text_read([name, '/', name(1:end-4), '_drugs.txt'], '%s');
%subjects = text_read([name, '/', name(1:end-4), '_subjects.txt'], '%s');
hrs = text_read([name, '/', name(1:end-4), '_hrs.txt'], '%s');
%fourhrs = text_read([name, '/', name(1:end-4), '_4hr_periods.txt'], '%s');
sixmins = text_read([name, '/', name(1:end-4), '_6mins.txt'], '%s');
states = text_read([name, '/', name(1:end-4), '_states.txt'], '%s');
load([name, '/', name, '_summed.mat'])

c_order = [.5 .5 .5; 1 .75 0; 1 0 1; 0 1 1];%[1 .6 0; .6 .4 .8; .59 .29 0; .21 .37 .23];

%% Barplots by 4 hr period.

% p=1;
% for i=1:1
%     fourhr_labels{p}=['Hours ',num2str(4*i),' to ',num2str(4*(i-1)+1),' Preinjection'];
%     short_fourhr_labels{p}=['pre',num2str(4*i),'to',num2str(4*(i-1)+1)];
%     p=p+1;
% end
% for i=1:4
%     fourhr_labels{p}=['Hours ',num2str(4*(i-1)+1),' to ',num2str(4*i),' Postinjection'];
%     short_fourhr_labels{p}=['post',num2str(4*(i-1)+1),'to',num2str(4*i)];
%     p=p+1;
% end
% 
% cplot_collected_BP_by_3_categories('Modulation Index, z-Scored by Hour'],[name,'/',name(1 : end - 4),'_summed_PLV_hr_by_4hr'],{band_labels, band_labels},{drug_labels, drug_labels},{state_labels, state_labels},{short_fourhr_labels, fourhr_labels},drugs,states,fourhrs,summed_PLV)
% 
% cplot_collected_BP_by_3_categories('Modulation Index, z-Scored by State Every 4 Hrs',[name,'/',name(1 : end - 4),'_summed_MI_4hr_by_4hr'],{band_labels, band_labels},{drug_labels, drug_labels},{state_labels, state_labels},{short_fourhr_labels, fourhr_labels},drugs,states,fourhrs,summed_MI_4hr)

%% Lineplots by hr period.

[hr_labels, short_hr_labels] = deal(cell(24, 1));

p=1;
for i=4:-1:1
    hr_labels{p}=['Hour ',num2str(i),' Preinjection'];
    short_hr_labels{p}=['pre',num2str(i)];
    p=p+1;
end
for i=1:20
    hr_labels{p}=['Hour ',num2str(i),' Postinjection'];
    short_hr_labels{p}=['post',num2str(i)];
    p=p+1;
end

% "Raw" summed MI.

lineplot_collected_BP_by_3_categories([title_name, ' PLV'],[name,'/',name(1 : end - 4),'_summed_PLV_hr_by_state'],{band_labels, band_labels},{drug_labels, drug_labels},{state_labels, state_labels},{short_hr_labels, hr_labels},drugs,states,hrs,summed_PLV,stat)

lineplot_collected_BP_by_categories([title_name, ' PLV'],[name,'/',name(1 : end - 4),'_summed_PLV_hr'],{band_labels, band_labels},{drug_labels, drug_labels},{short_hr_labels, hr_labels},drugs,hrs,summed_PLV,stat,c_order)

% Summed MI as a percentage of baseline.

lineplot_collected_BP_by_3_categories([title_name, ' PLV, Pct. Change from Baseline'],[name,'/',name(1 : end - 4),'_summed_PLV_pct_hr_by_state'],{band_labels, band_labels},{drug_labels, drug_labels},{state_labels, state_labels},{short_hr_labels, hr_labels},drugs,states,hrs,summed_PLV_pct,stat)

lineplot_collected_BP_by_categories([title_name, ' PLV, Pct. Change from Baseline'],[name,'/',name(1 : end - 4),'_summed_PLV_pct_hr'],{band_labels, band_labels},{drug_labels, drug_labels},{short_hr_labels, hr_labels},drugs,hrs,summed_PLV_pct,stat,c_order)

%% Lineplots by 6 min period.

no_pre=2; no_post=8;

[pd_labels, ~]=make_period_labels(no_pre,no_post,'6mins');

% "Raw" summed MI.

lineplot_collected_BP_by_3_categories([title_name, ' PLV'],[name,'/',name(1 : end - 4),'_summed_PLV_6min_by_state'],{band_labels, band_labels},{drug_labels, drug_labels},{state_labels, state_labels},{pd_labels, pd_labels},drugs,states,sixmins,summed_PLV,stat)

lineplot_collected_BP_by_categories([title_name, ' PLV'],[name,'/',name(1 : end - 4),'_summed_PLV_6min'],{band_labels, band_labels},{drug_labels, drug_labels},{pd_labels, pd_labels},drugs,sixmins,summed_PLV,stat,c_order)

% Summed MI as a percentage of baseline.

lineplot_collected_BP_by_3_categories([title_name, ' PLV, Pct. Change from Baseline'],[name,'/',name(1 : end - 4),'_summed_PLV_pct_6min_by_state'],{band_labels, band_labels},{drug_labels, drug_labels},{state_labels, state_labels},{pd_labels, pd_labels},drugs,states,sixmins,summed_PLV_pct,stat)

lineplot_collected_BP_by_categories([title_name, ' PLV, Pct. Change from Baseline'],[name,'/',name(1 : end - 4),'_summed_PLV_pct_6min'],{band_labels, band_labels},{drug_labels, drug_labels},{pd_labels, pd_labels},drugs,sixmins,summed_PLV_pct,stat,c_order)

end
