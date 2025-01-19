library(tidyverse)
test <- read.csv("~/Downloads/PupilBioTest_PMP_revA.csv", stringsAsFactors=TRUE)

test_edit <- test %>%
  separate(CpG_Coordinates, into = c("site1", "site2", "site3"), sep = ":", remove = FALSE, convert = FALSE)
save( test_edit, file="test_edit.RData")       

rm(test)           

site_coverage_summary <- function(test_edit, column_site) {
  if(column_site=="site1"){
    test_edit %>%
      mutate(site_cov = rowSums(test_edit[, c("X.100" , "X.101" , "X.110", "X.111" )])) %>%
      mutate(Total_reads = rowSums(test_edit[, c("X.000", "X.100" , "X.101" , "X.110", "X.111" )])) %>%
      rename(site=site1)%>%
      rename(Unmeth_site_cov=X.000)%>%
      select(-site2, -site3,-X.100, -X.010, -X.001, -X.110, -X.101, -X.011, -X.111 )
  }
  else if(column_site=="site2") {
    test_edit %>% 
      mutate(site_cov = rowSums(test_edit[, c("X.010" , "X.011" , "X.110", "X.111")])) %>%
      mutate(Total_reads = rowSums(test_edit[, c("X.000", "X.010" , "X.011" , "X.110", "X.111")])) %>%
      rename(site=site2)%>%
      rename(Unmeth_site_cov=X.000)%>%
      select(-site1, -site3,-X.100, -X.010, -X.001, -X.110, -X.101, -X.011, -X.111 )
  }
  else{
    test_edit %>%
      mutate(site_cov = rowSums(test_edit[, c("X.001" , "X.011", "X.101",  "X.111")])) %>%
      mutate(Total_reads = rowSums(test_edit[, c("X.000","X.001" , "X.011", "X.101",  "X.111")])) %>%
      rename(site=site3)%>%
      rename(Unmeth_site_cov=X.000)%>%
      select(-site1, -site2,-X.100, -X.010, -X.001, -X.110, -X.101, -X.011, -X.111 )
  }
  
}

site1 <- site_coverage_summary(test_edit , "site1") 

site2 <- site_coverage_summary(test_edit, "site2")


site3 <- site_coverage_summary(test_edit, "site3")
rm(test_edit)
site <- rbind(site1, site2, site3)
save(site, file ="site.RData")
rm(site1_)
rm(site2)
rm(site3)
##
site_proportion <- site %>%
  group_by(Tissue,Sample_ID,Replicate,strand,site) %>%
  summarise(methsite_coverage = sum(site_cov[site_cov != 0], na.rm = TRUE),Total_Reads = sum(Total_reads[Total_reads!= 0], na.rm = TRUE) )%>%
  group_by(Tissue,site,Sample_ID,Replicate) %>% 
  summarise(methsite_coverage_str = sum(methsite_coverage, na.rm = TRUE),Total_Reads_str = sum(Total_Reads, na.rm = TRUE)) %>%
  group_by(Tissue,site, Sample_ID) %>%
  summarise(site_coverage_mean = mean(methsite_coverage_str, na.rm = TRUE),Total_Reads_mean = mean(Total_Reads_str, na.rm = TRUE))%>%
  rowwise()%>%
  mutate(proportion = site_coverage_mean / Total_Reads_mean)
save(site_proportion, file = "site_proportion.RData")

##median coverage and Coefficient of variation
median_cv <- site_proportion %>%
  select(-Total_Reads_mean, -proportion) %>%
  group_by(Tissue,site) %>%
  summarise(site_median= median(site_coverage_mean, na.rm = TRUE),
            site_cv=(sd(site_coverage_mean) / abs(mean(site_coverage_mean)))*100, na.rm = TRUE)%>%
  filter(site_median != 0)
write.csv(median_cv, file="median_cv_CpG.csv")

ggplot(median_cv, aes(x=Tissue, y= site_median, fill=Tissue))+geom_boxplot() +
  labs( title = "median coverage of CpG sites",
        x = "Tissue",
        y = "coverage median")+theme(plot.title = element_text(size = 35), axis.text=element_text(size=35, colour = "black"), axis.title=element_text(size=35,face="bold"),legend.text =element_text(size=35), legend.title =element_text(size=35)  )

ggplot(median_cv, aes(x=Tissue, y= site_cv, fill=Tissue))+geom_boxplot()+
  labs( title = "coefficient of variation of CpG sites",
        x = "Tissue",
        y = "coefficient of variation")+theme(plot.title = element_text(size = 35), axis.text=element_text(size=35, colour = "black"), axis.title=element_text(size=35,face="bold"),legend.text =element_text(size=35), legend.title =element_text(size=35)  )




##significant CpGs and specificity analysis to compare to the top ten PMPs

library(purrr)

wilcox <- site_proportion %>%
  select(Tissue,site, Sample_ID,proportion)%>%
  group_by(site) %>%
  group_split() %>%
  purrr::imap_dfr(~ {
    if (n_distinct(.x$Tissue) == 2) {
      tibble(site = .x$site[1], 
             p_value = wilcox.test(proportion ~Tissue, exact = FALSE, data = .x)$p.value)
    } else {
      tibble(site = .x$site[1], p_value = NA)
    }
  })
  
wilcox$padj <- p.adjust(wilcox$p_value, method ="BH")
  save(wilcox, file =  "specificity_sitesCpG.RData")  

metrics_CpG <- site_proportion %>%
  group_by(Tissue,site) %>%
  summarise(across("site_coverage_mean":"proportion", ~mean(.)))



prop_CpG_cfDNA <- metrics_CpG %>%
  filter(Tissue=="cfDNA") %>%
  rename(proportion_cf = proportion,Total_Reads_cf = Total_Reads_mean, Variant_Reads_cf = site_coverage_mean)
prop_CpG_Islet <- metrics_CpG %>%
  filter(Tissue=="Islet")%>%
  rename(proportion_Is = proportion,Total_Reads_cf = Total_Reads_mean, Variant_Reads_cf = site_coverage_mean)
prop_diff <-  inner_join(prop_CpG_Islet,prop_CpG_cfDNA, by= "site")

Prop_differ <- prop_diff%>%
  mutate(prop_df = proportion_cf - proportion_Is)

final_result <- inner_join(wilcox,Prop_differ, by= "site") 
final_result <- final_result %>%
  select(-Tissue.y, -Tissue.x)

final_result_filtered <-  final_result %>%
  filter(padj<0.05) 

significant_sites <- final_result_filtered %>%
  mutate(Tissue = ifelse(prop_df > 0,"cfDNA", "Islet"))%>%
  mutate(classification = ifelse(prop_df >0.1,"cfDNA", "Islet"))  
library(caret)

conf_matrix <- confusionMatrix(as.factor(significant_sites$classification), 
                               as.factor(significant_sites$Tissue))

final_result_filterd <- significant_sites %>%
  select(-Tissue)%>%
  rename(Tissue = classification)
save(final_result_filterd, file= "CpGfinal_result_filtered.RData")











