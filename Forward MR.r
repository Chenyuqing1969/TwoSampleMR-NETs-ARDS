library(TwoSampleMR)
library(ieugwasr)
library(gwasvcf)
library(gwasglue)
library(VariantAnnotation)
library(tidyverse)
library(plinkbinr)
library(LDlinkR)
library(data.table)
library(MungeSumstats)

#暴露数据和工具变量筛选

exposure_data <- fread(file = "F:/GCST004448.gz",header = T, sep = '\t')
exposure_data <- as.data.frame(exposure_data)
colnames(exposure_data)
head(exposure_data)

# 筛选EAF > 1%
# 筛选P,一般情况默认为 5e-08，若按该标准筛选出IV少于5个，我们改为 5e-06

exposure_data  <- exposure_data [which(exposure_data$p_value<5e-06, exposure_data $effect_allele_frequency > 0.01),]


# 设置type = "exposure"
exposure_data <- format_data(dat = exposure_data, type = "exposure", 
                             snp_col = 'hm_rsid', chr_col = 'hm_chrom', 
                             pos_col = 'hm_pos',
                             eaf_col = 'effect_allele_frequency',
                             effect_allele_col = 'effect_allele', other_allele_col = 'other_allele',
                             beta_col = 'beta',se_col = 'standard_error',pval_col = 'p_value')


exposure_data <- clump_data(exposure_data, clump_r2=0.001, clump_kb = 10000)
exposure_data$exposure <- 'IL-1β'

#获取结局变量急性呼吸窘迫综合征
outcome_data <- fread("F:/finngen_R10_J10_ARDS.gz",header = T)

head(outcome_data)
outcome_data <- as.data.frame(outcome_data)
colnames(outcome_data)


outcome_data <- format_data(dat = outcome_data, type = "outcome", snps = exposure_data$SNP,
                            snp_col = 'rsids', chr_col = '#chrom', pos_col = 'pos',
                            eaf_col = 'af_alt', effect_allele_col = 'ref', other_allele_col = 'alt',
                            beta_col = 'beta',se_col = 'sebeta',pval_col = 'pval')

ex_name <- 'IL-1β'
exp_dat = exposure_data

##剔除混杂因素

library(FastTraitR) 
res=look_trait(rsids=unique(exp_dat$SNP), pval=1e-5)
confound <- res[which(res$trait %in% unique(res$trait[-grep('IL-1β',res$trait)])),]
fwrite(confound,file = paste0("./result/",ex_name,"_SNP_&_confound.csv"), row.names = FALSE)

exp_dat <- exp_dat[-which(exp_dat$SNP %in% confound$rsid & exp_dat$exposure == ex_name), ] #
dim(exp_dat)


diff <- setdiff(exp_dat$SNP, outcome_data$SNP)
diff


# 计算F值, 公式：F = R²*(SampleSize-2)/(1-R²)
exp_dat$r2 <- (2 * exp_dat$beta.exposure^2 * exp_dat$eaf.exposure * (1 - exp_dat$eaf.exposure)) / 
  (2 * exp_dat$beta.exposure^2 * exp_dat$eaf.exposure * (1 - exp_dat$eaf.exposure) + 
     2 * exp_dat$samplesize.exposure * exp_dat$eaf.exposure * (1 - exp_dat$eaf.exposure) * exp_dat$se.exposure^2)
r2 <- sum(exp_dat$r2)

#公式：F = R²*(SampleSize-2)/(1-R²)

exp_dat$F <- exp_dat$r2 * (exp_dat$samplesize.exposure - 2) / (1 - exp_dat$r2)
exp_dat_meanF <- mean(exp_dat$F)

exp_dat_meanF

# 对暴露和结果数据进行合并 
bind_data <- harmonise_data(
  exposure_dat=exp_dat,
  outcome_dat=outcome_data,
  action= 2
)

#write.csv(bind_data,file = paste0("./data/正向/",ex_name,"_merge_正向.csv"))


## MR分析
res <- mr(dat = bind_data, 
          method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
res_or <- generate_odds_ratios(mr_res = res)

res_or$pval <- round(res_or$pval, 2)
res_or$or <- round(res_or$or, 2)
res_or$or_lci95 <- round(res_or$or_lci95, 2)
res_or$or_uci95 <- round(res_or$or_uci95, 2)
res_or$or_ci <- paste(res_or$or,'(',res_or$or_lci95,'-',res_or$or_uci95,')')

write.csv(res_or,file = "F:/MR分析.csv")




## 异质性检验
res_hetergeneity <- mr_heterogeneity(bind_data)

## 水平多效性检验
res_pleiotropy <- mr_pleiotropy_test(bind_data) #

## MRPRESSO检验
library(MRPRESSO)
mrpresso <- mr_presso(data=bind_data, 
                      BetaOutcome="beta.outcome", BetaExposure="beta.exposure", 
                      SdOutcome="se.outcome",  SdExposure="se.exposure", 
                      OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                      SignifThreshold = 0.05, NbDistribution = 300, seed = NULL)
mrpresso
mrpresso_a <- as.data.frame(mrpresso[["Main MR results"]])    
mrpresso_a$OR <- round(exp(mrpresso_a$`Causal Estimate`), 4)
mrpresso_a$LOR <- round(exp(mrpresso_a$`Causal Estimate` - 1.96*mrpresso_a$Sd),4)
mrpresso_a$UOR <- round(exp(mrpresso_a$`Causal Estimate` + 1.96*mrpresso_a$Sd),4)
mrpresso_a$OR_CI <- paste(mrpresso_a$OR,'(',mrpresso_a$LOR,'-',mrpresso_a$UOR,')')

## 单个SNP效应分析
res_single <- mr_singlesnp(dat = bind_data)
res_single

res_loo <- mr_leaveoneout(bind_data)
p_loo <- mr_leaveoneout_plot(leaveoneout_results = res_loo)
p_loo

# sctter plor
p1 <-mr_scatter_plot(res, bind_data)
p1[[1]]

# forest plot
res_single <- mr_singlesnp(bind_data)#,all_method=c("mr_ivw", "mr_two_sample_ml")) to specify method used
p2 <- mr_forest_plot(res_single) 
p2[[1]]

# Leave-one-out分析
res_loo <- mr_leaveoneout(dat = bind_data)
p3 <- mr_leaveoneout_plot(leaveoneout_results = res_loo)
p3[[1]]

# Funnel plot
res_single <- mr_singlesnp(bind_data)
p4 <- mr_funnel_plot(res_single)
p4[[1]]

save.image(file =paste0("./result/正向MR分析.Rdata"))
