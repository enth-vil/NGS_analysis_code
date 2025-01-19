library(tidyverse)
library(purrr)
library(caret)

test <- read.csv("~/Downloads/PupilBioTest_PMP_revA.csv", stringsAsFactors=TRUE)

##Two separate data frames were created by filtering the Tissue because of huge size. 


## weights were assigned to each combination to calculate weighted score for each PMP

weights <- c("X.000" = 0, "X.001" = 1/3, "X.010" = 1/3, "X.011" = 2/3, 
             "X.100" = 1/3, "X.101" = 2/3, "X.110" = 2/3, "X.111" = 1)

## Tissue= cfDNA
## The reads counts were aggregated by summing across strands, calculating mean across replicate for each sample 
## read proportions were calculated by dividing read counts for each combination by total reads in that sample
## weighted scores were calculated by multiplying read proportion of each combination with weighted scores 

PMP_wt_scores <- test %>%
  group_by(Tissue, Sample_ID, Replicate, CpG_Coordinates) %>%
  summarise(across(starts_with("X"), ~sum(.)) ) %>%
  group_by(Tissue, Sample_ID, CpG_Coordinates) %>%
  summarise(across(starts_with("X"), ~sum(.)) )
save(PMP_wt_scores, file =  "PMP_wt_scores.RData")
wt_scores <- PMP_wt_scores%>%
  rowwise()%>%
  mutate(Total_Reads = sum(c_across(`X.000`:`X.111`))) %>%
  mutate(Variant_Reads = sum(c_across(`X.001`:`X.111`)))%>%
  mutate(weighted_scores = sum(c_across(`X.000`:`X.111`) * weights))

save(wt_scores, file =  "wt_scores.RData")

rm(test)
## The weighted scores were aggregated by taking calculating mean across samples for each PMP

specificity_sites <- wt_scores %>%
  select(Tissue, Sample_ID, CpG_Coordinates, weighted_scores)%>%
  group_by(CpG_Coordinates) %>%
  group_split() %>%
  purrr::imap_dfr(~ {
    if (n_distinct(.x$Tissue) == 2) {
      tibble(CpG_Coordinates = .x$CpG_Coordinates[1], 
             p_value = wilcox.test(weighted_scores ~ Tissue, exact = FALSE, data = .x)$p.value)
    } else {
      tibble(CpG_Coordinates = .x$CpG_Coordinates[1], p_value = NA)
    }
  })

specificity_sites$padj <- p.adjust(specificity_sites$p_value, method ="BH")
save(specificity_sites, file =  "specificity_sites.RData")

##mean VRF
metrics <- wt_scores %>%
  select(Tissue, Sample_ID, CpG_Coordinates,Total_Reads, Variant_Reads, weighted_scores) %>%
  rowwise()%>%
  mutate(VRF= Variant_Reads / Total_Reads)%>%
  group_by(Tissue, CpG_Coordinates) %>%
  summarise(across("Total_Reads":"VRF", ~mean(.))) 

meanVRF <- metrics %>%
  select(Tissue, CpG_Coordinates,VRF)%>%
  group_by(Tissue, CpG_Coordinates) %>%
  summarise(MVRF= mean(VRF))

metrics_cfDNA <- metrics %>% 
  filter(Tissue=="cfDNA")%>%
  rename(weighted_scores_cf = weighted_scores,Total_Reads_cf = Total_Reads, Variant_Reads_cf = Variant_Reads, VRF_cf = VRF)

metrics_Islet <- metrics %>% 
  filter(Tissue=="Islet")%>%
  rename(weighted_scores_Is = weighted_scores,Total_Reads_Is = Total_Reads, Variant_Reads_Is = Variant_Reads, VRF_Is = VRF)

diff_wt_scores <- inner_join(metrics_cfDNA, metrics_Islet, by="CpG_Coordinates")

difference <- diff_wt_scores %>%
  mutate(wt_score_diff = weighted_scores_cf - weighted_scores_Is)

  
result_cfDNA_specific <- inner_join(specificity_sites, difference , by ="CpG_Coordinates")
result_cfDNA_specific <- result_cfDNA_specific %>%
 select(-Tissue.y, -Tissue.x)
  
save(result_cfDNA_specific, file= "Result_cfDNA_specific.RData")
result_cfDNA_specific_filterd <-  result_cfDNA_specific %>%
    filter(padj<0.05) 
##specificity
significant_sites <- result_cfDNA_specific_filterd %>%
  mutate(Tissue = ifelse(wt_score_diff > 0,"cfDNA", "Islet"))%>%
  mutate(classification = ifelse(wt_score_diff >0.2,"cfDNA", "Islet"))  
 

conf_matrix <- confusionMatrix(as.factor(significant_sites$classification), 
                                 as.factor(significant_sites$Tissue))

result_cfDNA_specific_filterd <- significant_sites %>%
  select(-Tissue)%>%
  rename(Tissue = classification)
save(result_cfDNA_specific_filterd, file= "Result_cfDNA_specific_filtered.RData")

write.csv(result_cfDNA_specific_filterd, file= "result_filtered.csv")
## Here weighted score difference ((mean weighted score of cfDNA) - (mean weighted score of Islet)) is used as a measure of sequencing depth
## effect of sequencing depth on p_value

res <- cor.test(result_cfDNA_specific_filterd$wt_score_diff , result_cfDNA_specific_filterd$padj,  exact = FALSE ,               method = "spearman")

res

## The sites were divided into two groups based on median of weighted score difference and box plots were plotted
median_wt <- median(result_cfDNA_specific_filterd$wt_score_diff)

plot_data <- result_cfDNA_specific_filterd %>%
  mutate(wt_score = case_when(
    wt_score_diff <= median_wt ~ "Low_Difference",
    wt_score_diff > median_wt ~ "High_Difference"
  ))

## Boxplot of p-values by score difference bins
ggplot(plot_data, aes(x = wt_score, y = padj, fill=wt_score)) +scale_y_log10()+scale_fill_brewer(palette="Dark2")+
  labs( title = "sequencing depth and specificity conf",
        x = "Weighted Score Difference",
        y = "log transformed p-value") +geom_boxplot(outlier.stroke = 0.1) +theme(plot.title = element_text(size = 35), axis.text=element_text(size=35, colour = "black"), axis.title=element_text(size=35,face="bold"),legend.text =element_text(size=35), legend.title =element_text(size=35)  )

##calculate the coverage threshold at a depth of 1000000 for the top 10 PMPs to call Tissue 2 
result_Islet <- result_cfDNA_specific_filterd %>%
  filter(Tissue=="Islet")
  
res_ordered_Islet <- result_Islet[order(result_Islet$padj),]

top10 <- res_ordered_Islet[1:10,]
cov_scaled <- top10 %>%
   select(Tissue,CpG_Coordinates,p_value,padj,VRF_Is) %>%
   mutate(scaled_million = VRF_Is * 1000000)
treshold <- min(cov_scaled$scaled_million)
write.csv(cov_scaled, file= "Top10_treshold.csv")
