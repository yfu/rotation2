genome_length = 4411532; %MTB genome length.
num_bins = length(bins);

%Get actual binding sites from the database.
query = ['SELECT ce.chipped_gene_locus, cfp.impulse_coordinate, ROUND(cfp.impulse_height/ce.mean_coverage) '...
    'FROM tbnetwork_v5.chip_experiments ce '...
    'JOIN tbnetwork_v5.chip_regions cr ON cr.chip_experiment_id = ce.id '...
    'JOIN tbnetwork_v5.csdeconv_final_peaks cfp ON cfp.chip_region_id = cr.id '...
    'WHERE ce.include_in_network2 = 1 AND cr.passes_filter = 1 AND cfp.impulse_height/ce.mean_coverage >= 1 GROUP BY cfp.id ORDER BY ce.id, cfp.impulse_height'];
[tf, bs_center, bs_height] = mysql(query);

%Computationally predicted binding sites. Software - FIMO.
query = ['SELECT TF, center, nlog10p FROM tbnetwork_v5.motifs_fimo_top_mtb '...
    'WHERE CONCAT(CAST(fc AS CHAR), ''-'', CAST(lane AS CHAR)) IN '...
    '(SELECT CONCAT(CAST(flowcell_number AS CHAR), ''-'', CAST(flowcell_lane AS CHAR)) FROM tbnetwork_v5.chip_experiments WHERE include_in_network2 = 1 AND nlog10p >= 5)'];
[tf_fimo, bs_center_fimo, bs_height_fimo] = mysql(query);

%Generate random binding sites.
bs_center_random = round(random('unif', ones(length(bs_center),1), ones(length(bs_center),1)*genome_length));

bins = 1:10000:genome_length; %Bin size.
num_bs = hist(bs_center, bins); %Actual binding sites.
num_bs_random = hist(bs_center_random,bins); %Uniform distribution.
num_bs_fimo = hist(bs_center_fimo,bins); %Computationally predicted binding sites.

xlab = 0:200000:genome_length;
ylab1 = 0:10:max(max(num_bs),max(num_bs_fimo));

scrsz = get(0,'ScreenSize');
f = figure('OuterPosition', [1 40 scrsz(3) scrsz(4)-40], 'color', 'white');

subplot(3,1,1);
bar(bins, num_bs);
set(gca, 'XTick', xlab, 'XTickLabel', sprintf('%.1f|',xlab./1000000), 'YTick', ylab1, 'YTickLabel', ylab1);
title('number of experimentally found BS');
grid on; axis tight; ylim([0 max(ylab1)]);

subplot(3,1,2);
bar(bins, num_bs_random);
set(gca, 'XTick', xlab, 'XTickLabel', sprintf('%.1f|',xlab./1000000), 'YTick', ylab1, 'YTickLabel', ylab1);
title('number of randomly generated BS');
grid on; axis tight; ylim([0 max(ylab1)]);

subplot(3,1,3);
bar(bins, num_bs_fimo);
set(gca, 'XTick', xlab, 'XTickLabel', sprintf('%.1f|',xlab./1000000), 'YTick', ylab1, 'YTickLabel', ylab1);
title('number of computationally predicted BS');
grid on; axis tight; ylim([0 max(ylab1)]);
