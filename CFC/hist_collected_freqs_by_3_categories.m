function hist_collected_freqs_by_3_categories(Title,name,band_bins,band_names,c_order,cat3_labels,cat1_labels,cat2_labels,cat3_vec,cat1_vec,cat2_vec,spec)

% cat1_labels is for rows. cat2_labels is for columns. cat3_labels is for
% number of subplots.

no_bands = length(band_bins);

long_cat1_labels=cat1_labels{2};
long_cat2_labels=cat2_labels{2};
% long_cat3_labels=cat3_labels{2};

cat1_labels=cat1_labels{1};
cat2_labels=cat2_labels{1};
cat3_labels=cat3_labels{1};

no_cats1=length(cat1_labels);
no_cats2=length(cat2_labels);
no_cats3=length(cat3_labels);

spec_stats=cell(no_bands);

for b = 1:no_bands
    
    if isempty(band_bins{b})
        
        no_bins = 100;
        
    else
        
        no_bins = length(band_bins{b});
        
    end
    
    spec_stats{b} = nan(no_bins,no_cats1,no_cats2);
   
end


close('all')

for c3=1:no_cats3
    
    cat3=char(cat3_labels{c3});
    
    spec_cat3=spec(strcmp(cat3_vec,cat3),:);
    
    cat1_in_cat3=cat1_vec(strcmp(cat3_vec,cat3));
    
    cat2_in_cat3=cat2_vec(strcmp(cat3_vec,cat3));
    
    for c1=1:no_cats1
        
        cat1=char(cat1_labels{c1});
        
        spec_cat1=spec_cat3(strcmp(cat1_in_cat3,cat1),:);
        
        cat2_in_cat1=cat2_in_cat3(strcmp(cat1_in_cat3,cat1));
        
        for c2=1:no_cats2
            
            cat2=char(cat2_labels{c2});
            
            spec_cat2=spec_cat1(strcmp(cat2_in_cat1,cat2),:);
            
            if ~isempty(spec_cat2) && size(spec_cat2,1)>=5
                
                for b = 1:no_bands
                
                    if isempty(band_bins{b})
                        
                        band_bins{b} = linspace(min(spec_cat2(:,b)),max(spec_cat2(:,b)),100)';
                        
                    end
                    
                    h = hist(spec_cat2(:,b), band_bins{b});
                    
                    spec_stats{b}(:,c1,c2) = reshape(h/sum(h),length(band_bins{b}),1);
                    
                end
                
            else
                
                for  b = 1:no_bands

                    spec_stats{b}(:,c1,c2) = nan;
                
                end
                
            end
            
        end
        
        %% Plots by Band.
        
        for b=1:no_bands
            
            band_freqs=band_bins{b};
            
            band_width=length(band_freqs);

            tick_indices=1:ceil(band_width/10):band_width;
            
            band_ticks=band_freqs(tick_indices);
            
            if mean(band_ticks)<10
                
                band_labels=round(10*band_ticks)/10;
                
            else
                
                band_labels=round(band_ticks);
                
            end
            
            %% Line Plots by Band.
            
            figure(c1)
            
            set(gcf,'DefaultAxesColorOrder',c_order)
            
            subplot(no_cats3,no_bands,(c3-1)*no_bands+b)
            
            plot(band_freqs,reshape(spec_stats{b}(:,c1,:),band_width,no_cats2))
           
            set(gca,'XTick',band_ticks,'XTickLabel',band_labels)
            
            xlim([band_ticks(1) band_ticks(end)])
            
            if b==no_bands && c3==1
                
                legend(long_cat2_labels,'Location','BestOutside')
                
            end
                
            if b==1
                
                ylabel(cat3)
                
            end   
                
            if c3==1
                
                if b==ceil(no_bands/2)
                    
                    if ~isempty(band_names)
                
                        title({Title,['Histogram for ',char(long_cat1_labels{c1})],band_names{b}})
                
                    else
                        
                        title({Title,['Histogram for ',char(long_cat1_labels{c1})]})
                        
                    end
                        
                else
                    
                    title(band_names{b})
                
                end
                
            end
            
        end
        
    end
    
end

for c1=1:no_cats1
    
    saveas(c1,[name,'_',cat1_labels{c1},'_line.fig'])
    
end
