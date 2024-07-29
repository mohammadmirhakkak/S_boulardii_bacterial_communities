library(ggplot2)
library(ggalluvial)
library(dplyr)
library(ggpubr)
library(scales)
library(stringr)



smetana = read.csv("mohammadmirhakkak/S_boulardii_bacterial_communities/res/_detailed_sorted.csv",header = TRUE)

# yeast model ID
smetana[smetana$receiver=='yeast_smetana','receiver'] = "yeastGEM_v8.6.2"
smetana[smetana$donor=='yeast_smetana','donor'] = "yeastGEM_v8.6.2"

smetana_plots = list()

comm_ids <- list()
comm_ids[[1]] <- c("quad-green-643","quad-green-yeast-643")
comm_ids[[2]] <- c("tri-green-322","tri-green-yeast-322")
titles = c("Quad-643","Tri-322")

reshaped_df <- list()

for (i in 1:length(comm_ids)){
  community <- smetana %>%
    filter(community %in% comm_ids[[i]]) %>%
    filter(compound_name %in% c("L-Tryptophan",
                                "L-Arginine",
                                "L-Isoleucine",
                                "L-Alanine",
                                "L-Lysine",
                                "L-Phenylalanine",
                                "L-Leucine",
                                "L-Asparagine",
                                "L-Cysteine",
                                "L-Threonine",
                                "L-Histidine",
                                "L-Tyrosine",
                                "L-Methionine",
                                "L-Proline",
                                "L-Glutamate",
                                "L-Aspartate",
                                "Glycine",
                                "L-Glutamine",
                                "L-Serine",
                                "L-Valine"))
  
  if (nrow(community)==0){
    next
  }
  
  # filter for top-10 metabiltes
  #aa = unique(community$compound_name)
  #sum_smetana = vector()
  #for (a in 1:length(aa)){
  #  sum_smetana[i] <- sum(community[community$compound_name==aa[a],"smetana"])
  #}
  #aa_ordered = aa[order(sum_smetana, decreasing = T)]
  #aa_ordered = aa_ordered[1:10]
  #community <- community %>% filter(compound_name %in% aa_ordered)
  
  
  community_reshaped <- rbind(community,community)
  community_reshaped <- rbind(community_reshaped,community)
  community_reshaped$category <- c(community$receiver,community$compound_name,community$donor)
  community_reshaped$collection <- c(rep("Receiver",nrow(community)),rep("Metabolite",nrow(community)),rep("Donor",nrow(community)))
  community_reshaped$subject <- rep(seq(nrow(community)),3)
  yeast_presence <- ifelse(str_detect(community_reshaped$community, "yeast"), "yes", "no")
  community_reshaped$Yeast <- yeast_presence
  #community_reshaped$collection <- factor(community_reshaped$collection, levels = c("Receiver","Metabolite","Donor"))
  reshaped_df[[i]] <- community_reshaped
  p <- ggplot(community_reshaped, aes(y = smetana, x = collection, stratum = category, alluvium = subject)) +
    geom_alluvium(width = 0.45,aes(fill = Yeast, color = after_scale(fill)),alpha = 1) + 
    geom_stratum(width = 0.45) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          axis.title.y = element_blank(),
          axis.text.x = element_text(color = "black"),axis.title.x = element_blank()) +
    ylab("Smetana score") + ggtitle(titles[i]) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_manual(values = c("no" = "#6e829133",
                                   "yes" = "#00afff33"))
  p_built <- ggplot_build(p)
  strata <- p_built$data[[2]]$stratum
  italic_labels <- ifelse(str_detect(strata, " "), "italic", "plain")
  p <- p + ggfittext::geom_fit_text(stat = "stratum", min.size = 2,aes(label = after_stat(stratum)), size = 6, fontface = italic_labels)
  
  smetana_plots[[i]] <- p
  
}

ggexport(plotlist = smetana_plots,
         filename = "mohammadmirhakkak/S_boulardii_bacterial_communities/res/suppl_smetana_ribbon.pdf",width = 8.7,height = 5)

# calculate percentages
reshaped_df[[1]] <- as_tibble(reshaped_df[[1]])
sum_smetana <- sum(reshaped_df[[1]]$smetana)/3*2 # summation of donation and receive
reshaped_df[[1]] %>% filter(collection == "Donor") %>% filter(category == "B_longum_subsp_infantis") %>% filter(Yeast == "yes") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[1]] %>% filter(collection == "Donor") %>% filter(category == "B_longum_subsp_infantis") %>% filter(Yeast == "no") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[1]] %>% filter(collection == "Donor") %>% filter(category == "Lactobacillus_brevis_ATCC_367") %>% filter(Yeast == "yes") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[1]] %>% filter(collection == "Donor") %>% filter(category == "Lactobacillus_brevis_ATCC_367") %>% filter(Yeast == "no") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[1]] %>% filter(collection == "Donor") %>% filter(category == "Lactobacillus_buchneri_CD034") %>% filter(Yeast == "yes") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[1]] %>% filter(collection == "Donor") %>% filter(category == "Lactobacillus_buchneri_CD034") %>% filter(Yeast == "no") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[1]] %>% filter(collection == "Donor") %>% filter(category == "Lactobacillus_johnsonii_DPC_6026") %>% filter(Yeast == "yes") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[1]] %>% filter(collection == "Donor") %>% filter(category == "Lactobacillus_johnsonii_DPC_6026") %>% filter(Yeast == "no") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[1]] %>% filter(collection == "Donor") %>% filter(category == "yeastGEM_v8.6.2") %>% filter(Yeast == "yes") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[1]] %>% filter(collection == "Donor") %>% filter(category == "yeastGEM_v8.6.2") %>% filter(Yeast == "no") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[1]] %>% filter(collection == "Receiver") %>% filter(category == "B_longum_subsp_infantis") %>% filter(Yeast == "yes") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[1]] %>% filter(collection == "Receiver") %>% filter(category == "B_longum_subsp_infantis") %>% filter(Yeast == "no") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[1]] %>% filter(collection == "Receiver") %>% filter(category == "Lactobacillus_brevis_ATCC_367") %>% filter(Yeast == "yes") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[1]] %>% filter(collection == "Receiver") %>% filter(category == "Lactobacillus_brevis_ATCC_367") %>% filter(Yeast == "no") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[1]] %>% filter(collection == "Receiver") %>% filter(category == "Lactobacillus_buchneri_CD034") %>% filter(Yeast == "yes") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[1]] %>% filter(collection == "Receiver") %>% filter(category == "Lactobacillus_buchneri_CD034") %>% filter(Yeast == "no") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[1]] %>% filter(collection == "Receiver") %>% filter(category == "Lactobacillus_johnsonii_DPC_6026") %>% filter(Yeast == "yes") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[1]] %>% filter(collection == "Receiver") %>% filter(category == "Lactobacillus_johnsonii_DPC_6026") %>% filter(Yeast == "no") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[1]] %>% filter(collection == "Receiver") %>% filter(category == "yeastGEM_v8.6.2") %>% filter(Yeast == "yes") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[1]] %>% filter(collection == "Receiver") %>% filter(category == "yeastGEM_v8.6.2") %>% filter(Yeast == "no") %>% select(smetana) %>% sum()/sum_smetana

reshaped_df[[2]] <- as_tibble(reshaped_df[[2]])
sum_smetana <- sum(reshaped_df[[2]]$smetana)/3*2 # summation of donation and receive
reshaped_df[[2]] %>% filter(collection == "Donor") %>% filter(category == "Lactobacillus_buchneri_CD034") %>% filter(Yeast == "yes") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[2]] %>% filter(collection == "Donor") %>% filter(category == "Lactobacillus_buchneri_CD034") %>% filter(Yeast == "no") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[2]] %>% filter(collection == "Donor") %>% filter(category == "Lactobacillus_jensenii_SNUV360") %>% filter(Yeast == "yes") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[2]] %>% filter(collection == "Donor") %>% filter(category == "Lactobacillus_jensenii_SNUV360") %>% filter(Yeast == "no") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[2]] %>% filter(collection == "Donor") %>% filter(category == "Lactobacillus_reuteri_DSM_20016") %>% filter(Yeast == "yes") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[2]] %>% filter(collection == "Donor") %>% filter(category == "Lactobacillus_reuteri_DSM_20016") %>% filter(Yeast == "no") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[2]] %>% filter(collection == "Donor") %>% filter(category == "yeastGEM_v8.6.2") %>% filter(Yeast == "yes") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[2]] %>% filter(collection == "Donor") %>% filter(category == "yeastGEM_v8.6.2") %>% filter(Yeast == "no") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[2]] %>% filter(collection == "Receiver") %>% filter(category == "Lactobacillus_buchneri_CD034") %>% filter(Yeast == "yes") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[2]] %>% filter(collection == "Receiver") %>% filter(category == "Lactobacillus_buchneri_CD034") %>% filter(Yeast == "no") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[2]] %>% filter(collection == "Receiver") %>% filter(category == "Lactobacillus_jensenii_SNUV360") %>% filter(Yeast == "yes") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[2]] %>% filter(collection == "Receiver") %>% filter(category == "Lactobacillus_jensenii_SNUV360") %>% filter(Yeast == "no") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[2]] %>% filter(collection == "Receiver") %>% filter(category == "Lactobacillus_reuteri_DSM_20016") %>% filter(Yeast == "yes") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[2]] %>% filter(collection == "Receiver") %>% filter(category == "Lactobacillus_reuteri_DSM_20016") %>% filter(Yeast == "no") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[2]] %>% filter(collection == "Receiver") %>% filter(category == "yeastGEM_v8.6.2") %>% filter(Yeast == "yes") %>% select(smetana) %>% sum()/sum_smetana
reshaped_df[[2]] %>% filter(collection == "Receiver") %>% filter(category == "yeastGEM_v8.6.2") %>% filter(Yeast == "no") %>% select(smetana) %>% sum()/sum_smetana


reshaped_df[[1]] %>% filter(collection == "Donor") %>% filter(category == "yeastGEM_v8.6.2") %>% filter(Yeast == "yes") %>% select(smetana) %>% sum()/sum_smetana*2
reshaped_df[[1]] %>% filter(collection == "Receiver") %>% select(smetana) %>% sum()


############################
#### PUBLICATION FORMAT ####
############################
smetana = read.csv("mohammadmirhakkak/S_boulardii_bacterial_communities/res/_detailed_sorted.csv",header = TRUE)

smetana <- as_tibble(smetana)

smetana[smetana$receiver=='yeast_smetana','receiver'] = "yeastGEM_v8.6.2"
smetana[smetana$donor=='yeast_smetana','donor'] = "yeastGEM_v8.6.2"

# Amino acid name change
smetana[smetana$compound_name == "L-Tryptophan","compound_name"] = "Trp"
smetana[smetana$compound_name == "L-Arginine","compound_name"] = "Arg"
smetana[smetana$compound_name == "L-Isoleucine","compound_name"] = "Ile"
smetana[smetana$compound_name == "L-Alanine","compound_name"] = "Ala"
smetana[smetana$compound_name == "L-Lysine","compound_name"] = "Lys"
smetana[smetana$compound_name == "L-Phenylalanine","compound_name"] = "Phe"
smetana[smetana$compound_name == "L-Leucine","compound_name"] = "Leu"
smetana[smetana$compound_name == "L-Asparagine","compound_name"] = "Asn"
smetana[smetana$compound_name == "L-Cysteine","compound_name"] = "Cys"
smetana[smetana$compound_name == "L-Threonine","compound_name"] = "Thr"
smetana[smetana$compound_name == "L-Histidine","compound_name"] = "His"
smetana[smetana$compound_name == "L-Proline","compound_name"] = "Pro"
smetana[smetana$compound_name == "L-Glutamate","compound_name"] = "Glu"
smetana[smetana$compound_name == "L-Aspartate","compound_name"] = "Asp"
smetana[smetana$compound_name == "Glycine","compound_name"] = "Gly"
smetana[smetana$compound_name == "L-Glutamine","compound_name"] = "Gln"
smetana[smetana$compound_name == "L-Serine","compound_name"] = "Ser"
smetana[smetana$compound_name == "L-Valine","compound_name"] = "Val"
smetana[smetana$compound_name == "L-Methionine","compound_name"] = "Met"
smetana[smetana$compound_name == "L-Tyrosine","compound_name"] = "Tyr"

# based on labels in fig4a
comm_ids[[1]] <- c("quad-green-yeast-643")
comm_ids[[2]] <- c("tri-green-yeast-154")
comm_ids[[3]] <- c("tri-green-yeast-315")
comm_ids[[4]] <- c("tri-green-yeast-322")
comm_ids[[5]] <- c("tri-red-green-yeast-31")
comm_ids[[6]] <- c("tri-green-yeast-325")
comm_ids[[7]] <- c("tri-green-yeast-309")
comm_ids[[8]] <- c("tri-green-yeast-321")
comm_ids[[9]] <- c("tri-green-yeast-119")
comm_ids[[10]] <- c("tri-red-green-yeast-41")
comm_ids[[11]] <- c("tri-red-green-yeast-50")
comm_ids[[12]] <- c("tri-green-yeast-193")
comm_ids[[13]] <- c("tri-green-yeast-294")
comm_ids[[14]] <- c("tri-green-yeast-323")
comm_ids[[15]] <- c("quad-green-yeast-560")
comm_ids[[16]] <- c("tri-green-yeast-229")
comm_ids[[17]] <- c("tri-red-green-yeast-3")
comm_ids[[18]] <- c("tri-red-green-yeast-39")
comm_ids[[19]] <- c("tri-red-green-yeast-43")
comm_ids[[20]] <- c("tri-red-green-yeast-86")
comm_ids[[21]] <- c("tri-red-green-yeast-67")
comm_ids[[22]] <- c("tri-green-yeast-163")
comm_ids[[23]] <- c("tri-green-yeast-451")
comm_ids[[24]] <- c("tri-green-yeast-402")
comm_ids[[25]] <- c("tri-green-yeast-430")
comm_ids[[26]] <- c("tri-green-yeast-130")
comm_ids[[27]] <- c("quad-green-yeast-446")
comm_ids[[28]] <- c("quad-green-yeast-447")
comm_ids[[29]] <- c("quad-green-yeast-525")
comm_ids[[30]] <- c("tri-green-yeast-105")
comm_ids[[31]] <- c("tri-green-yeast-112")
comm_ids[[32]] <- c("tri-red-green-yeast-73")
comm_ids[[33]] <- c("tri-green-yeast-164")
comm_ids_num = c(c(1:23),c(24:33))


AA <- c("Trp","Arg","Ile","Ala","Lys","Phe","Leu","Asn","Cys","Thr",
        "His","Tyr","Met","Pro","Glu","Asp","Gly","Gln","Ser","Val")

SMETANA <- vector()
"Amino acid" <- vector()
Donor <- vector()
ID <- vector()
Type <- vector()
for (i in 1:length(comm_ids_num)){
  community <- smetana %>%
    filter(community %in% comm_ids[[comm_ids_num[i]]]) %>%
    filter(compound_name %in% AA)
  
  if (nrow(community)==0){
    next
  }
  
  # filter for top-10 metabiltes
  #aa = unique(community$compound_name)
  #sum_smetana = vector()
  #for (a in 1:length(aa)){
  #  sum_smetana[i] <- sum(community[community$compound_name==aa[a],"smetana"])
  #}
  #aa_ordered = aa[order(sum_smetana, decreasing = T)]
  #aa_ordered = aa_ordered[1:10]
  #community <- community %>% filter(compound_name %in% aa_ordered)
  
  for (j in 1:length(AA)){
    score_y <- community %>% filter(compound_name == AA[j]) %>% filter(donor == "yeastGEM_v8.6.2") %>% select(smetana) %>% sum()
    score_b <- community %>% filter(compound_name == AA[j]) %>% filter(donor != "yeastGEM_v8.6.2") %>% select(smetana) %>% sum()
    SMETANA <- c(SMETANA,score_y,score_b)
    Donor <- c(Donor,"Yeast","Bacteria")
    `Amino acid` <- c(`Amino acid` , rep(AA[j],2))
    ID <- c(ID, rep(comm_ids[[comm_ids_num[i]]],2))
    if (comm_ids_num[i] %in% 1:21){
      Type <- c(Type,rep("Cooperative",2))
    }
    if (comm_ids_num[i] %in% 22:25){
      Type <- c(Type,rep("Competitive",2))
    }
  }
  
}

donation <- tibble(ID,Type,`Amino acid`,Donor,SMETANA)
donation$media <- "complete"

line_group = rep(1:(nrow(donation)/2),each = 2)
donation$line_group = line_group
donor <- donation$Donor
point_annot <- ifelse(donor == "Yeast", "y", "b")
donation$point_annot = point_annot
p1 <- ggplot(donation,aes(x = Donor, y = SMETANA)) + 
  geom_boxplot(aes(fill = Donor), col = "black", show.legend = FALSE) +
  geom_line(aes(group = line_group, col = Type)) +
  geom_point(aes(fill = point_annot, shape = media), size = 1, color = "black", show.legend = FALSE) +
  ylab("SMETANA score") + 
  theme_classic() +
  facet_wrap(~`Amino acid`, ncol=length(AA)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7),
        #axis.ticks.x = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1,size = 7),
        text = element_text(size = 8)) +
  stat_compare_means(label = "p.signif",paired = TRUE,label.x = 1.4, size = 2) +
  scale_colour_manual(values = c("Bacteria" = "#6e829133",
                                 "Yeast" = "#00afff33",
                                 "Cooperative" = "gray",
                                 "Competitive" = "#cd534c")) +
  scale_fill_manual(values = c("Bacteria" = "#6e829133",
                               "Yeast" = "#00afff33",
                               "y" = "#B6F1FF",
                               "b" = "#6c818fff")) +
  scale_shape_manual(values = c("complete" = 21))


pdf("mohammadmirhakkak/S_boulardii_bacterial_communities/res/smetana_figS13_all_AA.pdf",width = 8.5,height = 2)
p1
dev.off()










############################
## SUPPL BOXPLOTs SMETANA ##
############################
smetana = read.csv("mohammadmirhakkak/S_boulardii_bacterial_communities/res/_detailed_sorted.csv",header = TRUE)

smetana <- as_tibble(smetana)

# metabolite name correction
smetana[smetana$compound_name == "Fe2+ mitochondria","compound_name"] = "Fe2+"
smetana[smetana$compound_name == "Riboflavin C17H20N4O6","compound_name"] = "Riboflavin"
smetana[smetana$compound_name == "Maltotriose C18H32O16","compound_name"] = "Maltotriose"
smetana[smetana$compound_name == "Indole 3 acetaldehyde C10H9NO","compound_name"] = "Indole 3 acetaldehyde"
smetana[smetana$compound_name == "Thymine C5H6N2O2","compound_name"] = "Thymine"
smetana[smetana$compound_name == "Thymidine C10H14N2O5","compound_name"] = "Thymidine"
smetana[smetana$compound_name == "N-Acetyl-D-glucosamine(anhydrous)N-Acetylmuramic acid","compound_name"] = "*anhgm"
smetana[smetana$compound_name == "Maltose C12H22O11","compound_name"] = "Maltose"
smetana[smetana$compound_name == "L Arabinitol C5H12O5","compound_name"] = "L Arabinitol"
smetana[smetana$compound_name == "Glycolate C2H3O3","compound_name"] = "Glycolate"
smetana[smetana$compound_name == "UMP C9H11N2O9P","compound_name"] = "UMP"
smetana[smetana$compound_name == "O2 O2","compound_name"] = "Oxygen"
smetana[smetana$compound_name == "Undecaprenyl-diphospho-N-acetylmuramoyl-(N-acetylglucosamine)-L-ala-D-glu-meso-2,6-diaminopimeloyl-D-ala-D-ala","compound_name"] = "*uaagmda"
smetana[smetana$compound_name == "Choline C5H14NO","compound_name"] = "Choline"
smetana[smetana$compound_name == "Three linked disacharide pentapeptide murein units (uncrosslinked, middle of chain)","compound_name"] = "*murein5p5p5p"
smetana[smetana$compound_name == "Taurine C2H7NO3S","compound_name"] = "Taurine"
smetana[smetana$compound_name == "Two linked disacharide pentapeptide murein units (uncrosslinked, middle of chain)","compound_name"] = "*murein5p5p"
smetana[smetana$compound_name == "Melibiose C12H22O11","compound_name"] = "Melibiose"
smetana[smetana$compound_name == "L-alanine-D-glutamate-meso-2,6-diaminoheptanedioate","compound_name"] = "*LalaDgluMdap"
smetana[smetana$compound_name == "Urea CH4N2O","compound_name"] = "Urea"
smetana[smetana$compound_name == "CO2 CO2","compound_name"] = "CO2"
smetana[smetana$compound_name == "H2O H2O","compound_name"] = "H2O"
smetana[smetana$compound_name == "Three disacharide linked murein units (pentapeptide crosslinked tetrapeptide (A2pm->D-ala) tetrapeptide corsslinked tetrapeptide (A2pm->D-ala)) (middle of chain)","compound_name"] = "*murein5px4px4p"
smetana[smetana$compound_name == "Two disacharide linked murein units, pentapeptide crosslinked tetrapeptide (A2pm->D-ala) (middle of chain)","compound_name"] = "*murein5px4p"
smetana[smetana$compound_name == "Two linked disacharide pentapeptide and tetrapeptide murein units (uncrosslinked, middle of chain)","compound_name"] = "*murein5p4p"
smetana[smetana$compound_name == " R R  2 3 Butanediol C4H10O2","compound_name"] = "(R,R)-2,3-Butanediol"
smetana[smetana$compound_name == "Sn-Glycero-3-phosphoethanolamine","compound_name"] = "*g3pe"
smetana[smetana$compound_name == "2-Dehydro-3-deoxy-D-gluconate","compound_name"] = "*2ddglcn"



# based on labels in fig4a
comm_ids <- list()
comm_ids[[1]] <- c("quad-green-643","quad-green-yeast-643")
comm_ids[[2]] <- c("tri-green-154","tri-green-yeast-154")
comm_ids[[3]] <- c("tri-green-315","tri-green-yeast-315")
comm_ids[[4]] <- c("tri-green-322","tri-green-yeast-322")
comm_ids[[5]] <- c("tri-red-green-31","tri-red-green-yeast-31")
comm_ids[[6]] <- c("tri-green-325","tri-green-yeast-325")
comm_ids[[7]] <- c("tri-green-309","tri-green-yeast-309")
comm_ids[[8]] <- c("tri-green-321","tri-green-yeast-321")
comm_ids[[9]] <- c("tri-green-119","tri-green-yeast-119")
comm_ids[[10]] <- c("tri-red-green-41","tri-red-green-yeast-41")
comm_ids[[11]] <- c("tri-red-green-50","tri-red-green-yeast-50")
comm_ids[[12]] <- c("tri-green-193","tri-green-yeast-193")
comm_ids[[13]] <- c("tri-green-294","tri-green-yeast-294")
comm_ids[[14]] <- c("tri-green-323","tri-green-yeast-323")
comm_ids[[15]] <- c("quad-green-560","quad-green-yeast-560")
comm_ids[[16]] <- c("tri-green-229","tri-green-yeast-229")
comm_ids[[17]] <- c("tri-red-green-3","tri-red-green-yeast-3")
comm_ids[[18]] <- c("tri-red-green-39","tri-red-green-yeast-39")
comm_ids[[19]] <- c("tri-red-green-43","tri-red-green-yeast-43")
comm_ids[[20]] <- c("tri-red-green-86","tri-red-green-yeast-86")
comm_ids[[21]] <- c("tri-red-green-67","tri-red-green-yeast-67")
comm_ids[[22]] <- c("tri-green-163","tri-green-yeast-163")
comm_ids[[23]] <- c("tri-green-451","tri-green-yeast-451")
comm_ids[[24]] <- c("tri-green-402","tri-green-yeast-402")
comm_ids[[25]] <- c("tri-green-430","tri-green-yeast-430")
comm_ids[[26]] <- c("tri-green-130","tri-green-yeast-130")
comm_ids[[27]] <- c("quad-green-446","quad-green-yeast-446")
comm_ids[[28]] <- c("quad-green-447","quad-green-yeast-447")
comm_ids[[29]] <- c("quad-green-525","quad-green-yeast-525")
comm_ids[[30]] <- c("tri-green-105","tri-green-yeast-105")
comm_ids[[31]] <- c("tri-green-112","tri-green-yeast-112")
comm_ids[[32]] <- c("tri-red-green-73","tri-red-green-yeast-73")
comm_ids[[33]] <- c("tri-green-164","tri-green-yeast-164")
comm_ids_num = c(c(1:23),c(24:33))
comm_sizes_wo_yeast = c(4,rep(3,13),4,rep(3,11),rep(4,3),rep(3,4))

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
    
    yeast_presence <- ifelse(str_detect(community$community, "yeast"), "yes", "no")
    community$Yeast <- yeast_presence
    community_yless <- community[community$Yeast=='no',]
    community_y <- community[community$Yeast=='yes',]
    
    sum_smetana[i] <- sum(community_yless[community_yless$compound_name==unq_met[m],'smetana']) / comm_sizes_wo_yeast[i]
    sum_smetana_y[i] <- sum(community_y[community_y$compound_name==unq_met[m],'smetana']) / (comm_sizes_wo_yeast[i]+1)
  }
  
  smetana_per_met[[m]] <- sum_smetana
  smetana_per_met_y[[m]] <- sum_smetana_y
  
}


SMETANA = vector()
Type = vector()
Yeast <- vector()
Metabolite = vector()
for (i in 1:length(unq_met)){
  #if (length(smetana_per_met[[i]])!=21){
  #  print(i)
  #}
  pos = smetana_per_met[[i]][1:21]
  neg = smetana_per_met[[i]][22:25]
  pos_y = smetana_per_met_y[[i]][1:21]
  neg_y = smetana_per_met_y[[i]][22:25]
  
  SMETANA = c(SMETANA,pos,neg,pos_y,neg_y)
  
  Type = c(Type,
           rep("Cooperative",length(pos)),
           rep("Competitive",length(neg)),
           rep("Cooperative",length(pos_y)),
           rep("Competitive",length(neg_y)))
  
  Yeast = c(Yeast,rep("Without yeast",length(c(pos,neg))),rep("With yeast",length(c(pos_y,neg_y))))
                    
  Metabolite = c(Metabolite,rep(unq_met[i],length(c(pos,neg,pos_y,neg_y))))

}

df_smetana_norm = tibble(Metabolite,Type,Yeast,SMETANA)
df_smetana_norm$media = "complete"
df_smetana_norm$Yeast <- factor(df_smetana_norm$Yeast, levels = c("Without yeast","With yeast"))

ggexport_plots <- list()
#len(unq_met) = 124; divide the number by 12 (No. plots per page)
for (i in 0:10){
  ggarrange_plots <- list()
  for (j in 1:12){
    if (i*12+j > length(unq_met)){
      break
    }
    df_filtered <- df_smetana_norm %>% filter(Metabolite == unq_met[i*12+j])
    line_group = rep(1:(nrow(df_filtered)/2),2)
    df_filtered$line_group = line_group
    yeast <- df_filtered$Yeast
    point_annot <- ifelse(yeast == "With yeast", "y", "n")
    df_filtered$point_annot = point_annot
    p <- ggplot(df_filtered,aes(x = Yeast, y = SMETANA)) + 
      geom_boxplot(aes(fill = Yeast), col = "black", show.legend = FALSE) +
      geom_line(aes(group = line_group, col = Type)) +
      geom_point(aes(fill = point_annot, shape = media), size = 1, color = "black", show.legend = FALSE) +
      ylab("Normalized SMETANA score") + 
      ggtitle(unq_met[i*12+j]) +
      theme_classic() +
      facet_wrap(~Type, ncol=2) +
      theme(strip.background = element_blank(),
            strip.text = element_text(size = 7),
            title = element_text(size = 7),
            axis.title.y = element_text(size = 7),
            #axis.ticks.x = element_blank(),
            legend.position = "none",
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust=1,size = 7),
            text = element_text(size = 8)) +
      stat_compare_means(label = "p.signif",paired = TRUE,label.x = 1.4, size = 2) +
      scale_colour_manual(values = c("Without yeast" = "#6e829133",
                                     "With yeast" = "#00afff33",
                                     "Cooperative" = "gray",
                                     "Competitive" = "#cd534c")) +
      scale_fill_manual(values = c("Without yeast" = "#6e829133",
                                   "With yeast" = "#00afff33",
                                   "y" = "#B6F1FF",
                                   "n" = "#6c818fff")) +
      scale_shape_manual(values = c("complete" = 21))
    
    ggarrange_plots[[j]] <- p
  }
  ggexport_plots[[i+1]] <- ggarrange(plotlist = ggarrange_plots,ncol = 4, nrow = 3)
}

ggexport(plotlist = ggexport_plots,
         filename = "mohammadmirhakkak/S_boulardii_bacterial_communities/res/Suppl_smetana_all_compounds.pdf",width = 8, height = 8)



#################################
######## Random vs Study ########
#################################

final_res_combined <- data.frame()
ggexport_plots <- list()
#20 randomized community sets; divide the number by 9 (No. plots per page)
for (i in 0:2){ # two pages
  ggarrange_plots <- list()
  for (j in 1:9){
    if (i*9+j > 20){
      break
    }
    
    org_res = read.csv("mohammadmirhakkak/S_boulardii_bacterial_communities/res/_global_sorted.csv")
    random_res_directory = paste0("mohammadmirhakkak/S_boulardii_bacterial_communities/res/random_communities_",j,"_global.tsv")
    random_res = read.table(random_res_directory,sep = '\t',header = TRUE)
    
    # numerize random_res
    random_res$size <- as.numeric(random_res$size)
    random_res$mip <- as.numeric(random_res$mip)
    random_res$mro <- as.numeric(random_res$mro)
    
    # add a column indicating presence/absence of yeast/fungi
    yeast_presence <- ifelse(str_detect(org_res$community, "yeast"), "yes", "no")
    org_res$with_yeast <- yeast_presence
    fungi_presence <- ifelse(str_detect(random_res$community, "_with_fungi"), "yes", "no")
    random_res$with_fungi <- fungi_presence
    
    # remove "-yeast"/"_with_fungi" from community name
    for (k in 1:nrow(org_res)){
      org_res$community[k] <- str_remove(org_res$community[k],'-yeast')
    }
    for (k in 1:nrow(random_res)){
      random_res$community[k] <- str_remove(random_res$community[k],'_with_fungi')
    }
    
    # split tables into with and without yeast/fungi
    org_res_with <- org_res[org_res$with_yeast=='yes',]
    org_res_with <- org_res_with[order(org_res_with$community),]
    org_res_without <- org_res[org_res$with_yeast=='no',]
    org_res_without <- org_res_without[order(org_res_without$community),]
    random_res_with <- random_res[random_res$with_fungi=='yes',]
    random_res_with <- random_res_with[order(random_res_with$community),]
    random_res_without <- random_res[random_res$with_fungi=='no',]
    random_res_without <- random_res_without[order(random_res_without$community),]
    
    # calculate log2(with/without)
    fc_mip_org <- log2(org_res_with$mip/org_res_without$mip)
    fc_mro_org <- log2(org_res_with$mro/org_res_without$mro)
    fc_mip_random <- log2(random_res_with$mip/random_res_without$mip)
    fc_mro_random <- log2(random_res_with$mro/random_res_without$mro)
    fc_org <- log2((org_res_with$mip/org_res_with$mro)/(org_res_without$mip/org_res_without$mro))
    fc_random <- log2((random_res_with$mip/random_res_with$mro)/(random_res_without$mip/random_res_without$mro))
    
    # make data table for visualization
    t1 <- org_res_with[,c("community","medium")]
    t1$score <- fc_mip_org
    t1$score_type <- "MIP"
    t1$group <- "Study"
    t1$modified_name <- c(1:nrow(t1))
    
    t2 <- org_res_with[,c("community","medium")]
    t2$score <- fc_mro_org
    t2$score_type <- "MRO"
    t2$group <- "Study"
    t2$modified_name <- c(1:nrow(t2))
    
    t3 <- org_res_with[,c("community","medium")]
    t3$score <- fc_org
    t3$score_type <- "MIP/MRO"
    t3$group <- "Study"
    t3$modified_name <- c(1:nrow(t3))
    
    t4 <- random_res_with[,c("community","medium")]
    t4$score <- fc_mip_random
    t4$score_type <- "MIP"
    t4$group <- "Random"
    t4$modified_name <- c((nrow(t1)+1):(nrow(t1)+nrow(t4)))
    
    t5 <- random_res_with[,c("community","medium")]
    t5$score <- fc_mro_random
    t5$score_type <- "MRO"
    t5$group <- "Random"
    t5$modified_name <- c((nrow(t2)+1):(nrow(t2)+nrow(t5)))
    
    t6 <- random_res_with[,c("community","medium")]
    t6$score <- fc_random
    t6$score_type <- "MIP/MRO"
    t6$group <- "Random"
    t6$modified_name <- c((nrow(t3)+1):(nrow(t3)+nrow(t6)))
    
    final_res <- rbind(t1,t2,t3,t4,t5,t6)
    final_res$group <- as.factor(final_res$group)
    final_res$random_community <- paste("Random_communities",i*9+j)
    
    final_res_combined <- rbind(final_res_combined,final_res)
    
    # save the results table
    #save_directory = paste0("Documents/S_boulardii_modeling/results/ISME_revision/random",j,"_vs_study.csv")
    #write.csv(final_res,save_directory)
    
    
    # box plots
    p <- ggplot(final_res,aes(x = group, y = score)) + 
      geom_boxplot(aes(fill = group), col = "black", show.legend = FALSE) +
      ylab("log2FC") + 
      ggtitle(paste("Random communities",9*i+j)) +
      theme_classic() +
      facet_wrap(~score_type, ncol=3) +
      theme(strip.background = element_blank(),
            strip.text = element_text(size = 7),
            title = element_text(size = 7),
            axis.title.y = element_text(size = 7),
            #axis.ticks.x = element_blank(),
            legend.position = "none",
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust=1,size = 7),
            text = element_text(size = 8)) +
      stat_compare_means(label = "p.signif",paired = FALSE,label.x = 1.4, size = 2)
    
    ggarrange_plots[[j]] <- p
  }
  ggexport_plots[[i+1]] <- ggarrange(plotlist = ggarrange_plots,ncol = 3, nrow = 3)
}

ggexport(plotlist = ggexport_plots,
         filename = "mohammadmirhakkak/S_boulardii_bacterial_communities/res/random_vs_study.pdf",width = 8, height = 8)

write.csv(final_res_combined,"mohammadmirhakkak/S_boulardii_bacterial_communities/res/random_vs_study.csv")
