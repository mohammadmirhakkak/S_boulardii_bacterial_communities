library(dplyr)
library(scales)
library(ggpubr)
library(RColorBrewer)

single_vs_pair_fva_sc <- read.csv("Documents/S_boulardii_modeling/results/pairwise_carveme_20230328/yeast_bacteria_pairwise_fva_single_pair_sc_bacSingleFVA_exchange.csv")
single_vs_pair_fva_mrs <- read.csv("Documents/S_boulardii_modeling/results/pairwise_carveme_231005/yeast_bacteria_fva_clean_mrs.csv")

single_vs_pair_fva_sc <- as_tibble(single_vs_pair_fva_sc)
single_vs_pair_fva_mrs <- as_tibble(single_vs_pair_fva_mrs)

# merge two datasets
single_vs_pair_fva_sc$media = 'SCmod'
single_vs_pair_fva_mrs$media = 'MRS'
single_vs_pair_fva <- rbind(single_vs_pair_fva_sc,single_vs_pair_fva_mrs)

unq_mets <- unique(single_vs_pair_fva_sc$metabolite_name)

single_vs_pair_fva_sc[single_vs_pair_fva_sc$simulation == 'single','model'] = 'bac_single'
single_vs_pair_fva_mrs[single_vs_pair_fva_mrs$simulation == 'single','model'] = 'bac_single'

single_vs_pair_fva_sc[is.na(single_vs_pair_fva_sc$value),'value'] = 0
single_vs_pair_fva_mrs[is.na(single_vs_pair_fva_mrs$value),'value'] = 0

pairs = list()
pairs[[1]] = c("bac_single","bac")
pairs[[2]] = c("bac_single","yeast")
pairs[[3]] = c("bac_single","joined")

bounds = c("lb","ub")

list_plots <- list()
for (i in 1:length(unq_mets)){
  
  filtered_df <- single_vs_pair_fva_sc %>% filter(metabolite_name == unq_mets[i]) %>% filter(model %in% pairs[[1]]) %>% filter(bound %in% bounds[1])
  p1 <- ggpaired(filtered_df, x = "model", y = "value",
                 line.color = "gray", line.size = 0.4,
                 color = "model", palette = "jco" ,short.panel.labs = FALSE) + 
    ggtitle(paste(unq_mets[i],"SC media",sep = '; ')) + ylab("LB") +
    theme(axis.ticks.x = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),text = element_text(size = 8)) +
    stat_compare_means(label = "p.format",paired = TRUE,size = 1.5)
  
  filtered_df <- single_vs_pair_fva_sc %>% filter(metabolite_name == unq_mets[i]) %>% filter(model %in% pairs[[2]]) %>% filter(bound %in% bounds[1])
  p2 <- ggpaired(filtered_df, x = "model", y = "value",
                 line.color = "gray", line.size = 0.4,
                 color = "model", palette = "jco" ,short.panel.labs = FALSE) + 
    ggtitle(paste(unq_mets[i],"SC media",sep = '; ')) + ylab("LB") +
    theme(axis.ticks.x = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),text = element_text(size = 8)) +
    stat_compare_means(label = "p.format",paired = TRUE,size = 1.5)
  
  filtered_df <- single_vs_pair_fva_sc %>% filter(metabolite_name == unq_mets[i]) %>% filter(model %in% pairs[[3]]) %>% filter(bound %in% bounds[1])
  p3 <- ggpaired(filtered_df, x = "model", y = "value",
                 line.color = "gray", line.size = 0.4,
                 color = "model", palette = "jco" ,short.panel.labs = FALSE) + 
    ggtitle(paste(unq_mets[i],"SC media",sep = '; ')) + ylab("LB") +
    theme(axis.ticks.x = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),text = element_text(size = 8)) +
    stat_compare_means(label = "p.format",paired = TRUE,size = 1.5)
  
  filtered_df <- single_vs_pair_fva_sc %>% filter(metabolite_name == unq_mets[i]) %>% filter(model %in% pairs[[1]]) %>% filter(bound %in% bounds[2])
  p4 <- ggpaired(filtered_df, x = "model", y = "value",
                 line.color = "gray", line.size = 0.4,
                 color = "model", palette = "jco" ,short.panel.labs = FALSE) + 
    ggtitle(paste(unq_mets[i],"SC media",sep = '; ')) + ylab("UB") +
    theme(axis.ticks.x = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),text = element_text(size = 8)) +
    stat_compare_means(label = "p.format",paired = TRUE,size = 1.5)
  
  filtered_df <- single_vs_pair_fva_sc %>% filter(metabolite_name == unq_mets[i]) %>% filter(model %in% pairs[[2]]) %>% filter(bound %in% bounds[2])
  p5 <- ggpaired(filtered_df, x = "model", y = "value",
                 line.color = "gray", line.size = 0.4,
                 color = "model", palette = "jco" ,short.panel.labs = FALSE) + 
    ggtitle(paste(unq_mets[i],"SC media",sep = '; ')) + ylab("UB") +
    theme(axis.ticks.x = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),text = element_text(size = 8)) +
    stat_compare_means(label = "p.format",paired = TRUE,size = 1.5)
  
  filtered_df <- single_vs_pair_fva_sc %>% filter(metabolite_name == unq_mets[i]) %>% filter(model %in% pairs[[3]]) %>% filter(bound %in% bounds[2])
  p6 <- ggpaired(filtered_df, x = "model", y = "value",
                 line.color = "gray", line.size = 0.4,
                 color = "model", palette = "jco" ,short.panel.labs = FALSE) + 
    ggtitle(paste(unq_mets[i],"SC media",sep = '; ')) + ylab("UB") +
    theme(axis.ticks.x = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),text = element_text(size = 8)) +
    stat_compare_means(label = "p.format",paired = TRUE,size = 1.5)
  
  filtered_df <- single_vs_pair_fva_mrs %>% filter(metabolite_name == unq_mets[i]) %>% filter(model %in% pairs[[1]]) %>% filter(bound %in% bounds[1])
  p7 <- ggpaired(filtered_df, x = "model", y = "value",
                 line.color = "gray", line.size = 0.4,
                 color = "model", palette = "jco" ,short.panel.labs = FALSE) + 
    ggtitle(paste(unq_mets[i],"MRS media",sep = '; ')) + ylab("LB") +
    theme(axis.ticks.x = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),text = element_text(size = 8)) +
    stat_compare_means(label = "p.format",paired = TRUE,size = 1.5)
  
  filtered_df <- single_vs_pair_fva_mrs %>% filter(metabolite_name == unq_mets[i]) %>% filter(model %in% pairs[[2]]) %>% filter(bound %in% bounds[1])
  p8 <- ggpaired(filtered_df, x = "model", y = "value",
                 line.color = "gray", line.size = 0.4,
                 color = "model", palette = "jco" ,short.panel.labs = FALSE) + 
    ggtitle(paste(unq_mets[i],"MRS media",sep = '; ')) + ylab("LB") +
    theme(axis.ticks.x = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),text = element_text(size = 8)) +
    stat_compare_means(label = "p.format",paired = TRUE,size = 1.5)
  
  filtered_df <- single_vs_pair_fva_mrs %>% filter(metabolite_name == unq_mets[i]) %>% filter(model %in% pairs[[3]]) %>% filter(bound %in% bounds[1])
  p9 <- ggpaired(filtered_df, x = "model", y = "value",
                 line.color = "gray", line.size = 0.4,
                 color = "model", palette = "jco" ,short.panel.labs = FALSE) + 
    ggtitle(paste(unq_mets[i],"MRS media",sep = '; ')) + ylab("LB") +
    theme(axis.ticks.x = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),text = element_text(size = 8)) +
    stat_compare_means(label = "p.format",paired = TRUE,size = 1.5)
  
  filtered_df <- single_vs_pair_fva_mrs %>% filter(metabolite_name == unq_mets[i]) %>% filter(model %in% pairs[[1]]) %>% filter(bound %in% bounds[2])
  p10 <- ggpaired(filtered_df, x = "model", y = "value",
                 line.color = "gray", line.size = 0.4,
                 color = "model", palette = "jco" ,short.panel.labs = FALSE) + 
    ggtitle(paste(unq_mets[i],"MRS media",sep = '; ')) + ylab("UB") +
    theme(axis.ticks.x = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),text = element_text(size = 8)) +
    stat_compare_means(label = "p.format",paired = TRUE,size = 1.5)
  
  filtered_df <- single_vs_pair_fva_mrs %>% filter(metabolite_name == unq_mets[i]) %>% filter(model %in% pairs[[2]]) %>% filter(bound %in% bounds[2])
  p11 <- ggpaired(filtered_df, x = "model", y = "value",
                 line.color = "gray", line.size = 0.4,
                 color = "model", palette = "jco" ,short.panel.labs = FALSE) + 
    ggtitle(paste(unq_mets[i],"MRS media",sep = '; ')) + ylab("UB") +
    theme(axis.ticks.x = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),text = element_text(size = 8)) +
    stat_compare_means(label = "p.format",paired = TRUE,size = 1.5)
  
  filtered_df <- single_vs_pair_fva_mrs %>% filter(metabolite_name == unq_mets[i]) %>% filter(model %in% pairs[[3]]) %>% filter(bound %in% bounds[2])
  p12 <- ggpaired(filtered_df, x = "model", y = "value",
                 line.color = "gray", line.size = 0.4,
                 color = "model", palette = "jco" ,short.panel.labs = FALSE) + 
    ggtitle(paste(unq_mets[i],"MRS media",sep = '; ')) + ylab("UB") +
    theme(axis.ticks.x = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),text = element_text(size = 8)) +
    stat_compare_means(label = "p.format",paired = TRUE,size = 1.5)
  
  multi_plot <- ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12, 
            #labels = c("A", "B", "C","D","E","F"),
            ncol = 3, nrow = 4)
  
  list_plots[[i]] <- multi_plot
  
}

ggexport(plotlist = list_plots,
         filename = "Documents/S_boulardii_modeling/results/Alex_20230711/ribbon_plot_alex_list/met_exchange_single_pair.pdf")







bounds = c("lb","ub")

# Box plots on S04 and S99 (best ones), S108, S07, S12 (worst ones)
# S04 -> L. brevis
# S99 -> L. crispatus
# S108 -> L. salivarus
# S07 -> L. gasseri
# S12 -> L. delbrueckii
single_vs_pair_fva_sc <- read.csv("Documents/S_boulardii_modeling/results/pairwise_carveme_20230328/yeast_bacteria_pairwise_fva_single_pair_sc_bacSingleFVA_exchange.csv")
single_vs_pair_fva_mrs <- read.csv("Documents/S_boulardii_modeling/results/pairwise_carveme_231005/yeast_bacteria_fva_clean_mrs.csv")

single_vs_pair_fva_sc <- as_tibble(single_vs_pair_fva_sc)
single_vs_pair_fva_mrs <- as_tibble(single_vs_pair_fva_mrs)

# merge two datasets
single_vs_pair_fva_sc$media = 'SCmod'
single_vs_pair_fva_mrs$media = 'MRS'
single_vs_pair_fva <- rbind(single_vs_pair_fva_sc,single_vs_pair_fva_mrs)

# make joined as the overall production of the metabolite and not the net exchange between the two microbes
single_vs_pair_fva[is.na(single_vs_pair_fva$value),'value'] = 0
single_vs_pair_fva[single_vs_pair_fva$value < 0,'value'] = 0
single_vs_pair_fva_1 <- single_vs_pair_fva %>% filter(simulation == 'pairwise')
single_vs_pair_fva_2 <- single_vs_pair_fva %>% filter(simulation == 'single')
single_vs_pair_fva = tibble()
for (i in 1:(nrow(single_vs_pair_fva_1)/6)){
  df <- single_vs_pair_fva_1[(i*6-5):(i*6),]
  df[5,'value'] = df[1,'value'] + df[3,'value']
  df[6,'value'] = df[2,'value'] + df[4,'value']
  single_vs_pair_fva = rbind(single_vs_pair_fva,df)
}
single_vs_pair_fva = rbind(single_vs_pair_fva,single_vs_pair_fva_2)


desired_mets <- c("Lactate","Acetate")

single_vs_pair_fva[single_vs_pair_fva$simulation == 'single','model'] = 'Mono-culture'
single_vs_pair_fva[single_vs_pair_fva$model == 'joined','model'] = 'Pairwise co-culture'

#pairs = list()
#pairs[[1]] = c("bac_single","bac")
#pairs[[2]] = c("bac_single","yeast")
#pairs[[3]] = c("bac_single","joined")

bounds = c("lb","ub")

desired_interactions <- c("yeast-Lactobacillus_brevis_ATCC_367",
                          "yeast-Lactobacillus_crispatus_ST1",
                          "yeast-L_salivarius",
                          "yeast-L_gasseri",
                          "yeast-L_delbrueckii_subsp_delbrueckii")

single_vs_pair_fva <- single_vs_pair_fva %>% filter(interaction %in% desired_interactions)

single_vs_pair_fva <- single_vs_pair_fva %>% filter(metabolite_name %in% desired_mets)

# add labels for good and bad species
single_vs_pair_fva[single_vs_pair_fva$interaction=='yeast-Lactobacillus_brevis_ATCC_367','interaction'] = "S04"
single_vs_pair_fva[single_vs_pair_fva$interaction=='yeast-Lactobacillus_crispatus_ST1','interaction'] = "S99"
single_vs_pair_fva[single_vs_pair_fva$interaction=='yeast-L_salivarius','interaction'] = "S108"
single_vs_pair_fva[single_vs_pair_fva$interaction=='yeast-L_gasseri','interaction'] = "S07"
single_vs_pair_fva[single_vs_pair_fva$interaction=='yeast-L_delbrueckii_subsp_delbrueckii','interaction'] = "S12"
simulation <- single_vs_pair_fva$simulation
point_annot <- ifelse(simulation == "pairwise", "pair_point", "single_point")
single_vs_pair_fva$point_annot = point_annot
single_vs_pair_fva$model <- factor(single_vs_pair_fva$model, levels = c("Mono-culture","Pairwise co-culture"))

colnames(single_vs_pair_fva)[4] = "Species"

dark_col_codes = brewer.pal(n=5,"Dark2")
# 00b0ffff -> dark blue
# 00afff33 -> light blue
# 6c818fff -> dark gray
# 6e829133 -> light gray



list_plots <- list()

filtered_df <- single_vs_pair_fva %>% filter(metabolite_name == "Lactate") %>% filter(model %in% c("Mono-culture","Pairwise co-culture")) %>% filter(bound %in% "lb")
p1 <- ggplot(filtered_df,aes(x = model, y = value)) + 
  geom_boxplot(aes(fill = model), col = "black", show.legend = FALSE) +
  geom_point(aes(fill = point_annot, shape = media), size = 2, color = "black", show.legend = FALSE) +
  geom_line(aes(group = Species, col = Species)) +
  ylab(paste("Minimum","Lactate","production","[mmol/grDW/hr]")) + 
  theme_classic() +
  facet_wrap(~media, ncol=2) +
  theme(strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size = 8)) +
  #stat_compare_means(label = "p.format",paired = TRUE,size = 1.5) +
  scale_colour_manual(values = c("S04" = dark_col_codes[1],
                                 "S99" = dark_col_codes[5],
                                 "S108" = dark_col_codes[2],
                                 "S07" = dark_col_codes[3],
                                 "S12" = dark_col_codes[4],
                                 "Mono-culture" = "#6e829133",
                                 "Pairwise co-culture" = "#00afff33")) +
  scale_fill_manual(values = c("Mono-culture" = "#6e829133",
                               "Pairwise co-culture" = "#00afff33",
                               "pair_point" = "#00b0ffff",
                               "single_point" = "#6c818fff")) +
  scale_shape_manual(values = c("SCmod" = 23,
                                "MRS" = 21))


filtered_df <- single_vs_pair_fva %>% filter(metabolite_name == "Lactate") %>% filter(model %in% c("Mono-culture","Pairwise co-culture")) %>% filter(bound %in% "ub")
p2 <- ggplot(filtered_df,aes(x = model, y = value)) + 
  geom_boxplot(aes(fill = model), col = "black", show.legend = FALSE) +
  geom_point(aes(fill = point_annot, shape = media), size = 2, color = "black", show.legend = FALSE) +
  geom_line(aes(group = Species, col = Species)) +
  ylab(paste("Maximum","Lactate","production","[mmol/grDW/hr]")) + 
  theme_classic() +
  facet_wrap(~media, ncol=2) +
  theme(strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size = 8)) +
  #stat_compare_means(label = "p.format",paired = TRUE,size = 1.5) +
  scale_colour_manual(values = c("S04" = dark_col_codes[1],
                                 "S99" = dark_col_codes[5],
                                 "S108" = dark_col_codes[2],
                                 "S07" = dark_col_codes[3],
                                 "S12" = dark_col_codes[4],
                                 "Mono-culture" = "#6e829133",
                                 "Pairwise co-culture" = "#00afff33")) +
  scale_fill_manual(values = c("Mono-culture" = "#6e829133",
                               "Pairwise co-culture" = "#00afff33",
                               "pair_point" = "#00b0ffff",
                               "single_point" = "#6c818fff")) +
  scale_shape_manual(values = c("SCmod" = 23,
                                "MRS" = 21))


filtered_df <- single_vs_pair_fva %>% filter(metabolite_name == "Acetate") %>% filter(model %in% c("Mono-culture","Pairwise co-culture")) %>% filter(bound %in% "lb")
p3 <- ggplot(filtered_df,aes(x = model, y = value)) + 
  geom_boxplot(aes(fill = model), col = "black", show.legend = FALSE) +
  geom_point(aes(fill = point_annot, shape = media), size = 2, color = "black", show.legend = FALSE) +
  geom_line(aes(group = Species, col = Species)) +
  ylab(paste("Minimum","Acetate","production","[mmol/grDW/hr]")) + 
  theme_classic() +
  facet_wrap(~media, ncol=2) +
  theme(strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size = 8)) +
  #stat_compare_means(label = "p.format",paired = TRUE,size = 1.5) +
  scale_colour_manual(values = c("S04" = dark_col_codes[1],
                                 "S99" = dark_col_codes[5],
                                 "S108" = dark_col_codes[2],
                                 "S07" = dark_col_codes[3],
                                 "S12" = dark_col_codes[4],
                                 "Mono-culture" = "#6e829133",
                                 "Pairwise co-culture" = "#00afff33")) +
  scale_fill_manual(values = c("Mono-culture" = "#6e829133",
                               "Pairwise co-culture" = "#00afff33",
                               "pair_point" = "#00b0ffff",
                               "single_point" = "#6c818fff")) +
  scale_shape_manual(values = c("SCmod" = 23,
                                "MRS" = 21))


filtered_df <- single_vs_pair_fva %>% filter(metabolite_name == "Acetate") %>% filter(model %in% c("Mono-culture","Pairwise co-culture")) %>% filter(bound %in% "ub")
p4 <- ggplot(filtered_df,aes(x = model, y = value)) + 
  geom_boxplot(aes(fill = model), col = "black", show.legend = FALSE) +
  geom_point(aes(fill = point_annot, shape = media), size = 2, color = "black", show.legend = FALSE) +
  geom_line(aes(group = Species, col = Species)) +
  ylab(paste("Maximum","Acetate","production","[mmol/grDW/hr]")) + 
  theme_classic() +
  facet_wrap(~media, ncol=2) +
  theme(strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size = 8)) +
  #stat_compare_means(label = "p.format",paired = TRUE,size = 1.5) +
  scale_colour_manual(values = c("S04" = dark_col_codes[1],
                                 "S99" = dark_col_codes[5],
                                 "S108" = dark_col_codes[2],
                                 "S07" = dark_col_codes[3],
                                 "S12" = dark_col_codes[4],
                                 "Mono-culture" = "#6e829133",
                                 "Pairwise co-culture" = "#00afff33")) +
  scale_fill_manual(values = c("Mono-culture" = "#6e829133",
                               "Pairwise co-culture" = "#00afff33",
                               "pair_point" = "#00b0ffff",
                               "single_point" = "#6c818fff")) +
  scale_shape_manual(values = c("SCmod" = 23,
                                "MRS" = 21))
  
  
  
multi_plot <- ggarrange(p1,p2,p3,p4,ncol = 2, nrow = 2,common.legend = T)
list_plots <- list()  
list_plots[[1]] <- multi_plot
  

ggexport(plotlist = list_plots,
         filename = "Documents/S_boulardii_modeling/results/pairwise_carveme_231005/ac_lac_production.pdf",width = 5, height = 6.7)






#################
# Box plots ALL #
#################
single_vs_pair_fva_sc <- read.csv("Documents/S_boulardii_modeling/results/pairwise_carveme_20230328/yeast_bacteria_pairwise_fva_single_pair_sc_bacSingleFVA_exchange.csv")
single_vs_pair_fva_mrs <- read.csv("Documents/S_boulardii_modeling/results/pairwise_carveme_231005/yeast_bacteria_fva_clean_mrs.csv")

single_vs_pair_fva_sc <- as_tibble(single_vs_pair_fva_sc)
single_vs_pair_fva_mrs <- as_tibble(single_vs_pair_fva_mrs)

# merge two datasets
single_vs_pair_fva_sc$media = 'SCmod'
single_vs_pair_fva_mrs$media = 'MRS'
single_vs_pair_fva <- rbind(single_vs_pair_fva_sc,single_vs_pair_fva_mrs)

# make joined as the overall production of the metabolite and not the net exchange between the two microbes
single_vs_pair_fva[is.na(single_vs_pair_fva$value),'value'] = 0
single_vs_pair_fva[single_vs_pair_fva$value < 0,'value'] = 0
single_vs_pair_fva_1 <- single_vs_pair_fva %>% filter(simulation == 'pairwise')
single_vs_pair_fva_2 <- single_vs_pair_fva %>% filter(simulation == 'single')
single_vs_pair_fva = tibble()
for (i in 1:(nrow(single_vs_pair_fva_1)/6)){
  df <- single_vs_pair_fva_1[(i*6-5):(i*6),]
  df[5,'value'] = df[1,'value'] + df[3,'value']
  df[6,'value'] = df[2,'value'] + df[4,'value']
  single_vs_pair_fva = rbind(single_vs_pair_fva,df)
}
single_vs_pair_fva = rbind(single_vs_pair_fva,single_vs_pair_fva_2)


desired_mets <- c("Lactate","Acetate")

single_vs_pair_fva[single_vs_pair_fva$simulation == 'single','model'] = 'Mono-culture'
single_vs_pair_fva[single_vs_pair_fva$model == 'joined','model'] = 'Pairwise co-culture'

#pairs = list()
#pairs[[1]] = c("bac_single","bac")
#pairs[[2]] = c("bac_single","yeast")
#pairs[[3]] = c("bac_single","joined")

bounds = c("lb","ub")

single_vs_pair_fva <- single_vs_pair_fva %>% filter(metabolite_name %in% desired_mets)


simulation <- single_vs_pair_fva$simulation
point_annot <- ifelse(simulation == "pairwise", "pair_point", "single_point")
single_vs_pair_fva$point_annot = point_annot
single_vs_pair_fva$model <- factor(single_vs_pair_fva$model, levels = c("Mono-culture","Pairwise co-culture"))

# 00b0ffff -> dark blue
# 00afff33 -> light blue
# 6c818fff -> dark gray
# 6e829133 -> light gray
list_plots <- list()

filtered_df <- single_vs_pair_fva %>% filter(metabolite_name == "Lactate") %>% filter(model %in% c("Mono-culture","Pairwise co-culture")) %>% filter(bound %in% "lb")
line_group = rep(seq(nrow(filtered_df)/2),2)
filtered_df$line_group = line_group
p1 <- ggplot(filtered_df,aes(x = model, y = value)) + 
  geom_boxplot(aes(fill = model), col = "black", show.legend = FALSE) +
  geom_point(aes(fill = point_annot, shape = media), size = 2, color = "black", show.legend = FALSE) +
  geom_line(aes(group = line_group), col = 'gray') +
  ylab(paste("Minimum","Lactate","production","[mmol/grDW/hr]")) + 
  theme_classic() +
  facet_wrap(~media, ncol=2) +
  theme(strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size = 8)) +
  stat_compare_means(label = "p.format",paired = TRUE,label.x = 1.4) +
  scale_colour_manual(values = c("S04" = dark_col_codes[1],
                                 "S99" = dark_col_codes[5],
                                 "S108" = dark_col_codes[2],
                                 "S07" = dark_col_codes[3],
                                 "S12" = dark_col_codes[4],
                                 "Mono-culture" = "#6e829133",
                                 "Pairwise co-culture" = "#00afff33")) +
  scale_fill_manual(values = c("Mono-culture" = "#6e829133",
                               "Pairwise co-culture" = "#00afff33",
                               "pair_point" = "#00b0ffff",
                               "single_point" = "#6c818fff")) +
  scale_shape_manual(values = c("SCmod" = 23,
                                "MRS" = 21))


filtered_df <- single_vs_pair_fva %>% filter(metabolite_name == "Lactate") %>% filter(model %in% c("Mono-culture","Pairwise co-culture")) %>% filter(bound %in% "ub")
line_group = rep(seq(nrow(filtered_df)/2),2)
filtered_df$line_group = line_group
p2 <- ggplot(filtered_df,aes(x = model, y = value)) + 
  geom_boxplot(aes(fill = model), col = "black", show.legend = FALSE) +
  geom_point(aes(fill = point_annot, shape = media), size = 2, color = "black", show.legend = FALSE) +
  geom_line(aes(group = line_group), col = 'gray') +
  ylab(paste("Maximum","Lactate","production","[mmol/grDW/hr]")) + 
  theme_classic() +
  facet_wrap(~media, ncol=2) +
  theme(strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size = 8)) +
  stat_compare_means(label = "p.format",paired = TRUE,label.x = 1.4) +
  scale_colour_manual(values = c("S04" = dark_col_codes[1],
                                 "S99" = dark_col_codes[5],
                                 "S108" = dark_col_codes[2],
                                 "S07" = dark_col_codes[3],
                                 "S12" = dark_col_codes[4],
                                 "Mono-culture" = "#6e829133",
                                 "Pairwise co-culture" = "#00afff33")) +
  scale_fill_manual(values = c("Mono-culture" = "#6e829133",
                               "Pairwise co-culture" = "#00afff33",
                               "pair_point" = "#00b0ffff",
                               "single_point" = "#6c818fff")) +
  scale_shape_manual(values = c("SCmod" = 23,
                                "MRS" = 21))


filtered_df <- single_vs_pair_fva %>% filter(metabolite_name == "Acetate") %>% filter(model %in% c("Mono-culture","Pairwise co-culture")) %>% filter(bound %in% "lb")
line_group = rep(seq(nrow(filtered_df)/2),2)
filtered_df$line_group = line_group
p3 <- ggplot(filtered_df,aes(x = model, y = value)) + 
  geom_boxplot(aes(fill = model), col = "black", show.legend = FALSE) +
  geom_point(aes(fill = point_annot, shape = media), size = 2, color = "black", show.legend = FALSE) +
  geom_line(aes(group = line_group), col = 'gray') +
  ylab(paste("Minimum","Acetate","production","[mmol/grDW/hr]")) + 
  theme_classic() +
  facet_wrap(~media, ncol=2) +
  theme(strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size = 8)) +
  stat_compare_means(label = "p.format",paired = TRUE,label.x = 1.4) +
  scale_colour_manual(values = c("S04" = dark_col_codes[1],
                                 "S99" = dark_col_codes[5],
                                 "S108" = dark_col_codes[2],
                                 "S07" = dark_col_codes[3],
                                 "S12" = dark_col_codes[4],
                                 "Mono-culture" = "#6e829133",
                                 "Pairwise co-culture" = "#00afff33")) +
  scale_fill_manual(values = c("Mono-culture" = "#6e829133",
                               "Pairwise co-culture" = "#00afff33",
                               "pair_point" = "#00b0ffff",
                               "single_point" = "#6c818fff")) +
  scale_shape_manual(values = c("SCmod" = 23,
                                "MRS" = 21))


filtered_df <- single_vs_pair_fva %>% filter(metabolite_name == "Acetate") %>% filter(model %in% c("Mono-culture","Pairwise co-culture")) %>% filter(bound %in% "ub")
line_group = rep(seq(nrow(filtered_df)/2),2)
filtered_df$line_group = line_group
p4 <- ggplot(filtered_df,aes(x = model, y = value)) + 
  geom_boxplot(aes(fill = model), col = "black", show.legend = FALSE) +
  geom_point(aes(fill = point_annot, shape = media), size = 2, color = "black", show.legend = FALSE) +
  geom_line(aes(group = line_group), col = 'gray') +
  ylab(paste("Maximum","Acetate","production","[mmol/grDW/hr]")) + 
  theme_classic() +
  facet_wrap(~media, ncol=2) +
  theme(strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size = 8)) +
  stat_compare_means(label = "p.format",paired = TRUE,label.x = 1.4) +
  scale_colour_manual(values = c("S04" = dark_col_codes[1],
                                 "S99" = dark_col_codes[5],
                                 "S108" = dark_col_codes[2],
                                 "S07" = dark_col_codes[3],
                                 "S12" = dark_col_codes[4],
                                 "Mono-culture" = "#6e829133",
                                 "Pairwise co-culture" = "#00afff33")) +
  scale_fill_manual(values = c("Mono-culture" = "#6e829133",
                               "Pairwise co-culture" = "#00afff33",
                               "pair_point" = "#00b0ffff",
                               "single_point" = "#6c818fff")) +
  scale_shape_manual(values = c("SCmod" = 23,
                                "MRS" = 21))



multi_plot <- ggarrange(p1,p2,p3,p4,ncol = 2, nrow = 2,common.legend = T)
list_plots <- list()  
list_plots[[1]] <- multi_plot


ggexport(plotlist = list_plots,
         filename = "Documents/S_boulardii_modeling/results/pairwise_carveme_231005/ac_lac_production_all.pdf",width = 5, height = 6.7)
