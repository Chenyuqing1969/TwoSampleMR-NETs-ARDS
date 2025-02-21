library(TwoSampleMR)
library(ggplot2)
library(tidyverse)
library(data.table)
library(MRPRESSO)
library(ieugwasr)
library(dplyr)
library(plinkbinr)
library(officer)
library(plyr)

#暴露数据和工具变量筛选
GWAS_data <- fread(paste0("./1.rawdata/finngen_R10_J10_ARDS.gz"))
GWAS_data <- as.data.frame(GWAS_data)
GWAS_data$n <- 357+406536
exp_summary <- format_data(dat = GWAS_data, type = "exposure", snps = NULL,
                           snp_col = 'rsids', chr_col = '#chrom', pos_col = 'pos',
                           eaf_col = 'af_alt', effect_allele_col = 'alt', other_allele_col = 'ref',
                           beta_col = 'beta',se_col = 'sebeta',pval_col = 'pval',samplesize_col = 'n')
exp_summary$exposure <- 'Acute Respiratory Distress Syndrome'  

ex_name <- unique(exp_summary$exposure)


# 筛选EAF > 1%
# 筛选P,一般情况默认为 5e-08，若按该标准筛选出IV少于5个，我们改为 5e-06
exp_summary1 <- filter(exp_summary,pval.exposure < 5e-08 )                 
if(dim(exp_summary1)[1] > 3){
  pval = 5e-08
}  else {
  pval = 5e-06
}              
exp_summary <- filter(exp_summary,pval.exposure < pval )    

snp_iv <- exp_summary %>% 
  select(rsid = SNP, pval = pval.exposure) %>% 
  ld_clump_local(., plink_bin = get_plink_exe(),
                 bfile = '/mnt/bak/yuejiali/SMR/g1000_eur/g1000_eur',
                 clump_r2 = 0.001, clump_kb = 10000, clump_p = 1) %>% 
  {.$rsid}
exp_iv <- filter(exp_summary, SNP %in% snp_iv)
exp_iv$maf.exposure <- ifelse(exp_iv$eaf.exposure > 0.5 ,1 - exp_iv$eaf.exposure, exp_iv$eaf.exposure)
exp_iv <- filter(exp_iv,maf.exposure > 0.01)
exp_dat <- exp_iv

##剔除混杂因素
library(FastTraitR) 
res=look_trait(rsids=unique(exp_dat$SNP), pval=1e-5)
confound <- res[which(res$trait %in% unique(res$trait[-grep('acute respiratory distress syndrome',res$trait)])),]
fwrite(confound,file = paste0("./result/反向/",gsub('\\/','\\.',ex_name),"_SNP_&_confound.csv"), row.names = FALSE)

exp_dat <- exp_dat[-which(exp_dat$SNP %in% confound$rsid & exp_dat$exposure == ex_name), ]

#公式：R² = 2*Beta^2*EAF*(1-EAF)/(2*Beta^2*EAF*(1-EAF)+ 2*SE^2*SampleSize*EAF*(1-EAF))
exp_dat$r2 <- (2 * exp_dat$beta.exposure^2 * exp_dat$eaf.exposure * (1 - exp_dat$eaf.exposure)) / 
  (2 * exp_dat$beta.exposure^2 * exp_dat$eaf.exposure * (1 - exp_dat$eaf.exposure) + 
     2 * exp_dat$samplesize.exposure * exp_dat$eaf.exposure * (1 - exp_dat$eaf.exposure) * exp_dat$se.exposure^2)
r2 <- sum(exp_dat$r2)

#公式：F = R²*(SampleSize-2)/(1-R²)

exp_dat$F <- exp_dat$r2 * (exp_dat$samplesize.exposure - 2) / (1 - exp_dat$r2)
exp_dat_meanF <- mean(exp_dat$F)

exp_dat_meanF
dir.create('./data/')
dir.create('./data/反向/')
write.csv(exp_dat,file = paste0("./data/反向/exposure_data_",ex_name,".csv"))


#获取结局变量：
GWAS_data <- fread(paste0('./1.rawdata/GCST90137414.h.tsv.gz'))
GWAS_data <- as.data.frame(GWAS_data)
GWAS_data$n <- 657
head(GWAS_data)

out_summary <- format_data(dat = GWAS_data, type = "outcome", snps = NULL,
                           snp_col = 'rsid', chr_col = 'chromosome', pos_col = 'base_pair_location',
                             eaf_col = 'effect_allele_frequency', effect_allele_col = 'effect_allele', other_allele_col = 'other_allele',
                             beta_col = 'beta',se_col = 'standard_error',pval_col = 'p_value',samplesize_col = 'n')
out_summary$outcome <- 'NETs'
out_names <- unique(out_summary$outcome)
out_summary = out_summary[which(out_summary$SNP != ''),]
out_summary <- unique(out_summary)
out_summary = out_summary[which(!duplicated(out_summary$SNP)),]

out_ids <- list.files('./1.rawdata/')[grep('fixed',list.files('./1.rawdata/'))][1:8]
out_names <- c('TNF-a','IL-18','IL-13','IL-6','IL-1β','IL-5','IL-4','Neutrophil count')
samplesize <- c(3454,3636,3557,8189,3309,3364,8124,519288)

for(i in 1:7){
    GWAS_data <- fread(paste0('./1.rawdata/',out_ids[i]))
    GWAS_data <- as.data.frame(GWAS_data)
    GWAS_data$n <- samplesize[i]
    out_summary <- format_data(dat = GWAS_data, type = "outcome", snps = NULL,
                           snp_col = 'variant_id', chr_col = 'chromosome', pos_col = 'base_pair_location',
                             eaf_col = 'EAF_rehab', effect_allele_col = 'effect_allele', other_allele_col = 'other_allele',
                             beta_col = 'beta',se_col = 'standard_error',pval_col = 'p_value',samplesize_col = 'n')
    out_summary$outcome <- out_names[i]
    out_summary = out_summary[which(out_summary$SNP != ''),]
    out_summary <- unique(out_summary)
    out_summary = out_summary[which(!duplicated(out_summary$SNP)),]
    out_iv <- dplyr::filter(.data = out_summary, SNP %in% exp_dat$SNP)
    library(LDlinkR)
    out_dat <- MR_proxy(exp = exp_dat, out = out_iv, out_summary = out_summary) 
    out_all <- rbind(out_all,out_dat)
}


GWAS_data <- fread(paste0('./1.rawdata/GCST90013658_fixed.tsv'))
GWAS_data <- as.data.frame(GWAS_data)
GWAS_data$n <- 5590
head(GWAS_data)

out_summary <- format_data(dat = GWAS_data, type = "outcome", snps = NULL,
                           snp_col = 'variant_id', chr_col = 'chromosome', pos_col = 'base_pair_location',
                           eaf_col = 'EAF_rehab', effect_allele_col = 'effect_allele', other_allele_col = 'other_allele',
                           beta_col = 'beta',se_col = 'standard_error',pval_col = 'p_value',samplesize_col = 'n')
out_summary$outcome <- 'MPO-DNA'
out_summary = out_summary[which(out_summary$SNP != ''),]
out_summary <- unique(out_summary)
out_summary = out_summary[which(!duplicated(out_summary$SNP)),]
out_iv <- dplyr::filter(.data = out_summary, SNP %in% exp_dat$SNP)

# 找代理SNP
library(LDlinkR)
MR_proxy <- function(exp, out, out_summary){
  
  missing_iv <- exp$SNP[!(exp$SNP %in% out$SNP)]
  print(paste0("Number of missing iv: ", as.character(length(missing_iv))))
  
  for (i in 1:length(missing_iv)) {
    print(paste0(missing_iv[i]))
  }
  
  if (length(missing_iv) == 0) {
    print("All exposure IVs found in outcome GWAS.")
  } else {
    print("Some exposure IVs missing from outcome GWAS.")
    
    for (i in 1:length(missing_iv)) {
      
      tryCatch({
        proxies <- LDproxy(snp = missing_iv[i], pop = "EUR", r2d = 'r2', token 
                           = "b4344efd0281", file = FALSE)
        
        proxies <- proxies[proxies$R2 > 0.8, ]
        proxy_present <- FALSE
        proxies$Alleles <- gsub("[()]","",proxies$Alleles)
        proxies <- proxies[which(grepl("-",proxies$Alleles) == F),]
        for (j in 1:nrow(proxies)) {
          if(proxies$RS_Number[j] %in% out_summary$SNP){
            proxy_present <- TRUE
            proxy_SNP <- proxies$RS_Number[j]
            alleles <- str_split(proxies$Alleles[j], "/", simplify = TRUE)
            proxy_SNP_allele_1 <- alleles[1, 1]
            proxy_SNP_allele_2 <- alleles[1, 2]
            original_alleles <- str_split(proxies$Alleles[1], "/", simplify = 
                                            TRUE)
            original_SNP_allele_1 <- original_alleles[1, 1]
            original_SNP_allele_2 <- original_alleles[1, 2]
            break
          }
        }
        
        if(proxy_present){
          print(paste0("Proxy SNP found. ", missing_iv[i], " replaced with ", 
                       proxy_SNP))
          proxy_row <- out[1, ]
          proxy_row$SNP <- missing_iv[i]
          matched_out_summary <- out_summary %>% filter(SNP == proxy_SNP)
          proxy_row$beta.outcome <- 
            as.numeric(matched_out_summary$beta.outcome)
          proxy_row$se.outcome <- as.numeric(matched_out_summary$se.outcome)
          proxy_row <- proxy_row %>%
            mutate(effect_allele.outcome = 
                     if_else(matched_out_summary$effect_allele.outcome == proxy_SNP_allele_1, 
                             original_SNP_allele_1, original_SNP_allele_2),
                   other_allele.outcome = 
                     if_else(matched_out_summary$other_allele.outcome == proxy_SNP_allele_1, 
                             original_SNP_allele_1, original_SNP_allele_2),
                   pval.outcome = as.numeric(matched_out_summary$pval.outcome),
                   chr.outcome = as.numeric(exp$chr.exposure[exp$SNP == 
                                                               missing_iv[i]]),
                   pos.outcome = as.numeric(exp$pos.exposure[exp$SNP == 
                                                               missing_iv[i]]),
                   `proxy.outcome` = TRUE,
                   `target_snp.outcome` = missing_iv[i],
                   `proxy_snp.outcome` = proxy_SNP)
          if("eaf.outcome" %in% colnames(out_summary)) {
            proxy_row$eaf.outcome <- 
              as.numeric(matched_out_summary$eaf.outcome)
          }
          if("samplesize.outcome" %in% colnames(out_summary)) {
            samplesize.outcome <- 
              as.numeric(matched_out_summary$samplesize.outcome)
          }
          out <- rbind.fill(out,proxy_row)
          print(paste0("Added proxy for ", missing_iv[i]))
        } else {
          print(paste0("No proxy SNP available for ", missing_iv[i], " in 
outcome GWAS."))
        }
      }, error = function(e){
        cat(e$message, '\n')
      })
    }
  }
  return(out)
}
out_dat <- MR_proxy(exp = exp_dat, out = out_iv, out_summary = out_summary) 
out_all <- rbind(out_all,out_dat)

# 对暴露和结果数据进行合并
bind_data <- harmonise_data(
  exposure_dat=exp_dat,
  outcome_dat=out_all,
  action= 2
)

#write.csv(bind_data,file = paste0("./data/反向/",ex_name,"_merge_反向.csv"))

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


save.image(file =paste0("./result/反向MR分析.Rdata"))
