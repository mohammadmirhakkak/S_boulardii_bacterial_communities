library(dplyr)
library(scales)
library(ggpubr)
library(stringr)


smetana = read.csv("Documents/S_boulardii_modeling/results/meeting_alex_20230331/smetana_complete/_detailed_sorted.csv",header = TRUE)
smetana_two = read.csv("Documents/S_boulardii_modeling/results/Alex_20230711/smetana/_detailed_compound_name.csv",header = TRUE,row.names = 1)
smetana = rbind(smetana,smetana_two)

smetana_plots = list()

# select two colors for positive and negative communities
#hex <- hue_pal()(2)
#neg_col_code = hex[1]
#pos_col_code = hex[2]


smetana[smetana$donor=="B_longum_subsp_infantis","donor"] = "B. longum"
smetana[smetana$receiver=="B_longum_subsp_infantis","receiver"] = "B. longum"
smetana[smetana$donor=="Lactobacillus_johnsonii_DPC_6026","donor"] = "L. johnsonii"
smetana[smetana$receiver=="Lactobacillus_johnsonii_DPC_6026","receiver"] = "L. johnsonii"
smetana[smetana$donor=="yeast_smetana","donor"] = "Yeast"
smetana[smetana$receiver=="yeast_smetana","receiver"] = "Yeast"
smetana[smetana$donor=="Lactobacillus_buchneri_CD034","donor"] = "L. buchneri"
smetana[smetana$receiver=="Lactobacillus_buchneri_CD034","receiver"] = "L. buchneri"
smetana[smetana$donor=="Lactobacillus_brevis_ATCC_367","donor"] = "L. brevis"
smetana[smetana$receiver=="Lactobacillus_brevis_ATCC_367","receiver"] = "L. brevis"
smetana[smetana$donor=="Lactobacillus_crispatus_125_2_CHN","donor"] = "L. crispatus"
smetana[smetana$receiver=="Lactobacillus_crispatus_125_2_CHN","receiver"] = "L. crispatus"
smetana[smetana$donor=="Lactobacillus_crispatus_ST1","donor"] = "L. crispatus"
smetana[smetana$receiver=="Lactobacillus_crispatus_ST1","receiver"] = "L. crispatus"
smetana[smetana$donor=="Lactobacillus_reuteri_DSM_20016","donor"] = "L. reuteri"
smetana[smetana$receiver=="Lactobacillus_reuteri_DSM_20016","receiver"] = "L. reuteri"
smetana[smetana$donor=="Lactobacillus_jensenii_SNUV360","donor"] = "L. jensenii"
smetana[smetana$receiver=="Lactobacillus_jensenii_SNUV360","receiver"] = "L. jensenii"
smetana[smetana$donor=="L_acidophilus","donor"] = "L. acidophilus"
smetana[smetana$receiver=="L_acidophilus","receiver"] = "L. acidophilus"
smetana[smetana$donor=="L_gasseri","donor"] = "L. gasseri"
smetana[smetana$receiver=="L_gasseri","receiver"] = "L. gasseri"
smetana[smetana$donor=="Lactobacillus_rhamnosus_GG_GG_ATCC_53103","donor"] = "L. rhamnosus"
smetana[smetana$receiver=="Lactobacillus_rhamnosus_GG_GG_ATCC_53103","receiver"] = "L. rhamnosus"
smetana[smetana$donor=="L_paracasei","donor"] = "L. paracasei"
smetana[smetana$receiver=="L_paracasei","receiver"] = "L. paracasei"
smetana[smetana$donor=="L_delbrueckii_subsp_delbrueckii","donor"] = "L. delbrueckii"
smetana[smetana$receiver=="L_delbrueckii_subsp_delbrueckii","receiver"] = "L. delbrueckii"

smetana[smetana$compound_name == "D-Lactate","compound_name"] = "Lactate"
smetana[smetana$compound_name == "L-Lactate","compound_name"] = "Lactate"
  
comm_ids <- list()
comm_ids[[1]] <- c("quad-green-643","quad-green-yeast-643")
comm_ids[[2]] <- c("tri-green-154","tri-green-yeast-154")
comm_ids[[3]] <- c("tri-green-315","tri-green-yeast-315")
comm_ids[[4]] <- c("tri-green-322","tri-green-yeast-322")
comm_ids[[5]] <- c("two-green-74","two-green-yeast-74")
comm_ids[[6]] <- c("tri-green-325","tri-green-yeast-325")
comm_ids[[7]] <- c("tri-green-309","tri-green-yeast-309")
comm_ids[[8]] <- c("tri-green-321","tri-green-yeast-321")
comm_ids[[9]] <- c("tri-green-119","tri-green-yeast-119")
comm_ids[[10]] <- c("tri-red-green-41","tri-red-green-yeast-41")
comm_ids[[13]] <- c("tri-green-130","tri-green-yeast-130")
comm_ids[[14]] <- c("quad-green-446","quad-green-yeast-446")
comm_ids[[15]] <- c("tri-red-green-31","tri-red-green-yeast-31")
comm_ids[[16]] <- c("tri-green-353","tri-green-yeast-353")
comm_ids[[17]] <- c("quad-green-1070","quad-green-yeast-1070")
comm_ids[[18]] <- c("tri-red-green-50","tri-red-green-yeast-50")
comm_ids[[20]] <- c("one-green-yeast-7")
comm_ids[[21]] <- c("one-green-yeast-2")
comm_ids[[22]] <- c("one-green-yeast-6")
comm_ids[[23]] <- c("one-green-yeast-5")
comm_ids[[24]] <- c("one-green-yeast-8")
comm_ids_num = c(c(1:10),c(13:18),c(20:24))
comm_sizes_wo_yeast = c(4,3,3,3,2,3,3,3,3,3,3,4,3,3,4,3,1,1,1,1,1)

#####################
####  ####
#####################
unq_met = unique(smetana$compound_name)
smetana_per_met_y = list()
smetana_per_met = list()
for (m in 1:124){# length(unq_met)){
  
  sum_smetana = vector()
  sum_smetana_y = vector()
  
  for (i in 1:length(comm_ids_num)){
    
    community <- smetana %>%
      filter(community %in% comm_ids[[comm_ids_num[i]]])
    
    if (nrow(community)==0){
      next
    }
    
    community_reshaped <- rbind(community,community)
    community_reshaped <- rbind(community_reshaped,community)
    community_reshaped$category <- c(community$receiver,community$compound_name,community$donor)
    community_reshaped$collection <- c(rep("Receiver",nrow(community)),rep("Metabolite",nrow(community)),rep("Donor",nrow(community)))
    community_reshaped$subject <- rep(seq(nrow(community)),3)
    yeast_presence <- ifelse(str_detect(community_reshaped$community, "yeast"), "yes", "no")
    community_reshaped$Yeast <- yeast_presence
    community_reshaped_yless <- community_reshaped[community_reshaped$Yeast=='no',]
    community_reshaped_y <- community_reshaped[community_reshaped$Yeast=='yes',]
    
    sum_smetana[i] <- sum(community_reshaped_yless[community_reshaped_yless$compound_name==unq_met[m],'smetana']) / comm_sizes_wo_yeast[i]
    sum_smetana_y[i] <- sum(community_reshaped_y[community_reshaped_y$compound_name==unq_met[m],'smetana']) / (comm_sizes_wo_yeast[i]+1)
  }
  
  smetana_per_met[[m]] <- sum_smetana
  smetana_per_met_y[[m]] <- sum_smetana_y
  
}


mean_pos_y = vector()
mean_neg_y = vector()
mean_pos = vector()
mean_neg = vector()
p_neg_y_pos_y = vector()
p_neg_pos = vector()
p_pos_pos_y = vector()
p_neg_neg_y = vector()
smetana_score = vector()
community_type = vector()
metabolites = vector()
for (i in 1:124){
  #if (length(smetana_per_met[[i]])!=21){
  #  print(i)
  #}
  pos = smetana_per_met[[i]][1:10]
  neg = smetana_per_met[[i]][11:16]
  pos_y = smetana_per_met_y[[i]][1:10]
  neg_y = smetana_per_met_y[[i]][11:16]
  
  smetana_score = c(smetana_score,pos,neg,pos_y,neg_y)
  community_type = c(community_type,
                     rep("pos_no_yeast",length(pos)),
                     rep("neg_no_yeast",length(neg)),
                     rep("pos_with_yeast",length(pos_y)),
                     rep("neg_with_yeast",length(neg_y)))
  metabolites = c(metabolites,rep(unq_met[i],length(c(pos,neg,pos_y,neg_y))))
  
  test = wilcox.test(neg_y,pos_y)
  p_neg_y_pos_y[i] = test$p.value
  test = wilcox.test(neg,pos)
  p_neg_pos[i] = test$p.value
  test = wilcox.test(pos,pos_y)
  p_pos_pos_y[i] = test$p.value
  test = wilcox.test(neg,neg_y)
  p_neg_neg_y[i] = test$p.value
  
  mean_pos_y[i] = mean(pos_y)
  mean_neg_y[i] = mean(neg_y)
  mean_pos[i] = mean(pos)
  mean_neg[i] = mean(neg)
}

df_stat_mets = tibble(unq_met,mean_neg,mean_pos,mean_neg_y,mean_pos_y,p_neg_pos,p_neg_y_pos_y,p_neg_neg_y,p_pos_pos_y)
df_stat_mets$q_neg_pos = p.adjust(df_stat_mets$p_neg_pos)
df_stat_mets$q_neg_y_pos_y = p.adjust(df_stat_mets$p_neg_y_pos_y)
df_stat_mets$q_neg_neg_y = p.adjust(df_stat_mets$p_neg_neg_y)
df_stat_mets$q_pos_pos_y = p.adjust(df_stat_mets$p_pos_pos_y)

df_stat_mets %>% filter(p_neg_y_pos_y < 0.05) %>% filter(p_neg_pos >= 0.05) -> selected_stat_mets
df_stat_mets %>% filter(q_neg_y_pos_y < 0.2) %>% filter(q_neg_pos >= 0.2)


df_stat_mets %>% filter(unq_met == "Acetaldehyde")
df_stat_mets %>% filter(unq_met == "Acetate")
df_stat_mets %>% filter(unq_met == "Lactate")



df_mets = tibble(community_type,smetana_score,metabolites)
list_plots <- list()
selected_mets <- c(selected_stat_mets$unq_met,"Acetaldehyde","Acetate","Lactate")
for (i in 1:length(selected_mets)){ 
  df_single_met <- df_mets %>% filter(metabolites == selected_mets[[i]])
  my_comparisons <- list(c("pos_no_yeast","pos_with_yeast"),
                         c("neg_no_yeast","neg_with_yeast"),
                         c("pos_no_yeast","neg_no_yeast"),
                         c("pos_with_yeast","neg_with_yeast"))
  p <- ggboxplot(df_single_met, x = "community_type", y = "smetana_score",
                 color = "community_type", palette = "jco", facet.by = "metabolites",
                 add = "jitter")
  #  Add p-value
  list_plots[[i]] <- p + stat_compare_means(comparisons = my_comparisons)
}

ggexport(plotlist = list_plots,
         filename = "Documents/S_boulardii_modeling/results/Alex_20230711/ribbon_plot_alex_list/met_transfer_difference.pdf")





mean_pos_y = vector()
mean_neg_y = vector()
mean_pos = vector()
mean_neg = vector()
p_pos_pos_y = vector()
p_neg_neg_y = vector()
smetana_score = vector()
community_type = vector()
metabolites = vector()
for (i in 1:124){
  #if (length(smetana_per_met[[i]])!=21){
  #  print(i)
  #}
  pos = smetana_per_met[[i]][1:10]
  neg = smetana_per_met[[i]][11:16]
  pos_y = smetana_per_met_y[[i]][1:10]
  neg_y = smetana_per_met_y[[i]][11:16]
  
  smetana_score = c(smetana_score,pos,neg,pos_y,neg_y)
  community_type = c(community_type,
                     rep("pos_no_yeast",length(pos)),
                     rep("neg_no_yeast",length(neg)),
                     rep("pos_with_yeast",length(pos_y)),
                     rep("neg_with_yeast",length(neg_y)))
  metabolites = c(metabolites,rep(unq_met[i],length(c(pos,neg,pos_y,neg_y))))
  
  
  test = wilcox.test(pos,pos_y,paired = TRUE)
  p_pos_pos_y[i] = test$p.value
  test = wilcox.test(neg,neg_y,paired = TRUE)
  p_neg_neg_y[i] = test$p.value
  
  mean_pos_y[i] = mean(pos_y)
  mean_neg_y[i] = mean(neg_y)
  mean_pos[i] = mean(pos)
  mean_neg[i] = mean(neg)
}

df_stat_mets = tibble(unq_met,mean_neg,mean_pos,mean_neg_y,mean_pos_y,p_neg_neg_y,p_pos_pos_y)
df_stat_mets$q_neg_neg_y = p.adjust(df_stat_mets$p_neg_neg_y)
df_stat_mets$q_pos_pos_y = p.adjust(df_stat_mets$p_pos_pos_y)

df_stat_mets %>% filter(p_pos_pos_y < 0.05) %>% filter(p_neg_neg_y >= 0.05) -> selected_stat_mets
df_stat_mets %>% filter(q_pos_pos_y < 0.2) %>% filter(q_neg_neg_y >= 0.2)

df_stat_mets %>% filter(unq_met == "Acetaldehyde")
df_stat_mets %>% filter(unq_met == "Acetate")
df_stat_mets %>% filter(unq_met == "Lactate")

### boxplot paired test ###
df_mets = tibble(community_type,smetana_score,metabolites)
impact = ifelse(str_detect(df_mets$community_type, "pos"), "positive", "negative")
yeast = ifelse(str_detect(df_mets$community_type, "with"), "yes", "no")
df_mets$impact = impact
df_mets$yeast = yeast
list_plots <- list()
selected_mets <- c(selected_stat_mets$unq_met,"Acetaldehyde","Lactate")
for (i in 1:length(selected_mets)){ 
  df_single_met <- df_mets %>% filter(metabolites == selected_mets[[i]])
  #my_comparisons <- list(c("pos_no_yeast","pos_with_yeast"),
  #                       c("neg_no_yeast","neg_with_yeast"))
  p <- ggpaired(df_single_met, x = "yeast", y = "smetana_score",
                 line.color = "gray", line.size = 0.4,
                 color = "yeast", palette = "jco", facet.by = "impact",short.panel.labs = FALSE) + 
    ggtitle(selected_mets[[i]]) + ylab("SMETANA score") +
    theme(axis.ticks.x = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank())
  #  Add p-value
  list_plots[[i]] <- p + stat_compare_means(label = "p.format",paired = TRUE)
}

ggexport(plotlist = list_plots,
         filename = "Documents/S_boulardii_modeling/results/Alex_20230711/ribbon_plot_alex_list/met_transfer_difference_paired.pdf")


## acetate and lactate for publication ##
box_annot <- ifelse(str_detect(df_mets$community_type, "with_yeast"), "With yeast", "Without yeast")
df_mets$box_annot = box_annot
df_mets$box_annot <- factor(df_mets$box_annot, levels = c("Without yeast","With yeast"))
line_annot <- ifelse(str_detect(df_mets$community_type, "pos_"), "positive", "negative")
df_mets$line_annot = line_annot
point_annot = rep("no_yeast",nrow(df_mets))
point_annot[df_mets$community_type=="pos_with_yeast"] = "pos_with_yeast"
point_annot[df_mets$community_type=="neg_with_yeast"] = "neg_with_yeast"
df_mets$point_annot = point_annot

filtered_df <- df_mets %>% filter(metabolites == "Acetate")
line_group = rep(seq(16),2)
filtered_df$line_group = line_group
p1 <- ggplot(filtered_df,aes(x = box_annot, y = smetana_score)) + 
  geom_boxplot(aes(fill = box_annot), col = "black", show.legend = FALSE) +
  geom_point(aes(fill = point_annot), size = 2, shape = 21, color = "black", show.legend = FALSE) +
  geom_line(aes(group = line_group, col = line_annot),show.legend = FALSE) +
  ylab(paste("Acetate","exchange SMETANA score")) + 
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size = 8)) +
  stat_compare_means(label = "p.format",paired = TRUE, label.x = 1.4) +
  scale_colour_manual(values = c("Without yeast" = "#6e829133",
                                 "With yeast" = "#00afff33",
                                 "positive" = "gray",
                                 "negative" = "#cd534c")) +
  scale_fill_manual(values = c("Without yeast" = "#6e829133",
                               "With yeast" = "#00afff33",
                               "pos_with_yeast" = "#00b0ffff",
                               "no_yeast" = "#6c818fff",
                               "neg_with_yeast" = "#cd534c"))


filtered_df <- df_mets %>% filter(metabolites == "Lactate")
line_group = rep(seq(16),2)
filtered_df$line_group = line_group
p2 <- ggplot(filtered_df,aes(x = box_annot, y = smetana_score)) + 
  geom_boxplot(aes(fill = box_annot), col = "black", show.legend = FALSE) +
  geom_point(aes(fill = point_annot), size = 2, shape = 21, color = "black", show.legend = FALSE) +
  geom_line(aes(group = line_group, col = line_annot),show.legend = FALSE) +
  ylab(paste("Lactate","exchange SMETANA score")) + 
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size = 8)) +
  stat_compare_means(label = "p.format",paired = TRUE, label.x = 1.4) +
  scale_colour_manual(values = c("Without yeast" = "#6e829133",
                                 "With yeast" = "#00afff33",
                                 "positive" = "gray",
                                 "negative" = "#cd534c")) +
  scale_fill_manual(values = c("Without yeast" = "#6e829133",
                               "With yeast" = "#00afff33",
                               "pos_with_yeast" = "#00b0ffff",
                               "no_yeast" = "#6c818fff",
                               "neg_with_yeast" = "#cd534c"))

multi_plot <- ggarrange(p1,p2,ncol = 2, nrow = 1)
list_plots <- list()  
list_plots[[1]] <- multi_plot


ggexport(plotlist = list_plots,
         filename = "Documents/S_boulardii_modeling/results/pairwise_carveme_231005/ac_lac_smetana.pdf",width = 4, height = 3.5)


# save df_mets
community_id <- vector()
for (i in c(1:16)){
  community_id[i] <- comm_ids[[comm_ids_num[i]]][1]
}
for (i in c(1:16)){
  community_id[i+16] <- comm_ids[[comm_ids_num[i]]][2]
}
df_mets$community_id = rep(community_id,length(unq_met))
write.csv(df_mets,"Documents/S_boulardii_modeling/results/Alex_20230711/ribbon_plot_alex_list/met_transfer_difference_paired.csv")


### fisher's test pos vs neg ###
fisher_pvalue <- vector()
for (i in 1:length(unq_met)){ 
  df_single_met <- df_mets %>% filter(metabolites == unq_met[i])
  quan_thr <- quantile(df_single_met$smetana_score,0.25)
  pos_with_yeast <- df_single_met %>% filter(community_type=="pos_with_yeast") %>% filter(smetana_score > quan_thr) %>% nrow()
  pos_no_yeast <- df_single_met %>% filter(community_type=="pos_no_yeast") %>% filter(smetana_score > quan_thr) %>% nrow()
  neg_with_yeast <- df_single_met %>% filter(community_type=="neg_with_yeast") %>% filter(smetana_score > quan_thr) %>% nrow()
  neg_no_yeast <- df_single_met %>% filter(community_type=="neg_no_yeast") %>% filter(smetana_score > quan_thr) %>% nrow()
  possible_exchange <-
    matrix(c(pos_with_yeast, pos_no_yeast, neg_with_yeast, neg_no_yeast),
           nrow = 2,
           dimnames = list(positive = c("with yeast", "without yeast"),
                           negative = c("with yeast", "without yeast")))
  fisher_pvalue[i] <- fisher.test(possible_exchange)$p.value
}
fisher_results = data_frame(unq_met,fisher_pvalue)
fisher_results %>% filter(fisher_pvalue < 0.05)

### fisher's test with yeast vs without yeast ###
fisher_pvalue <- vector()
for (i in 1:length(unq_met)){ 
  df_single_met <- df_mets %>% filter(metabolites == unq_met[i])
  quan_thr <- quantile(df_single_met$smetana_score,0.25)
  pos_with_yeast <- df_single_met %>% filter(community_type=="pos_with_yeast") %>% filter(smetana_score > quan_thr) %>% nrow()
  pos_no_yeast <- df_single_met %>% filter(community_type=="pos_no_yeast") %>% filter(smetana_score > quan_thr) %>% nrow()
  neg_with_yeast <- df_single_met %>% filter(community_type=="neg_with_yeast") %>% filter(smetana_score > quan_thr) %>% nrow()
  neg_no_yeast <- df_single_met %>% filter(community_type=="neg_no_yeast") %>% filter(smetana_score > quan_thr) %>% nrow()
  possible_exchange <-
    matrix(c(pos_with_yeast, neg_with_yeast, pos_no_yeast, neg_no_yeast),
           nrow = 2,
           dimnames = list(with_yeast = c("positive", "negative"),
                           without_yeast = c("positive", "negative")))
  fisher_pvalue[i] <- fisher.test(possible_exchange)$p.value
}
fisher_results = data_frame(unq_met,fisher_pvalue)
fisher_results %>% filter(fisher_pvalue < 0.05)







