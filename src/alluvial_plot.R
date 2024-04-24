library(ggplot2)
library(ggalluvial)
library(dplyr)



smetana = read.csv("Documents/S_boulardii_modeling/results/meeting_alex_20230331/smetana_complete/_detailed_sorted.csv",header = TRUE)

tri <- smetana %>%
  filter(community %in% c("tri-green-229","tri-green-yeast-229")) %>%
  filter(compound_name %in% c("D-Lactate","Acetate","Acetaldehyde"))

tri[tri$donor=="Lactobacillus_farciminis_KCTC_3681_DSM_20184","donor"] = "L. farciminis"
tri[tri$receiver=="Lactobacillus_farciminis_KCTC_3681_DSM_20184","receiver"] = "L. farciminis"
tri[tri$donor=="Lactobacillus_johnsonii_DPC_6026","donor"] = "L. johnsonii"
tri[tri$receiver=="Lactobacillus_johnsonii_DPC_6026","receiver"] = "L. johnsonii"
tri[tri$donor=="yeast_smetana","donor"] = "Yeast"
tri[tri$receiver=="yeast_smetana","receiver"] = "Yeast"
tri[tri$donor=="Lactobacillus_brevis_ATCC_367","donor"] = "L. brevis"
tri[tri$receiver=="Lactobacillus_brevis_ATCC_367","receiver"] = "L. brevis"

pdf("Documents/S_boulardii_modeling/results/Alex_20230711/ac_transfer_v2.pdf")

ggplot(as.data.frame(tri),
       aes(y = smetana,
           axis1 = donor, axis2 = compound_name, axis3 = receiver)) +
  geom_alluvium(aes(fill = community),
                width = 1/8, knot.pos = 0, reverse = FALSE) +
  #scale_fill_manual(values = c(Brown = "#70493D", Hazel = "#E2AC76",
  #                             Green = "#3F752B", Blue = "#81B0E4")) +
  guides(fill = "none") +
  geom_stratum(alpha = .25, width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            reverse = FALSE) +
  scale_x_continuous(breaks = 1:3, labels = c("Donor", "Metabolite", "Receiver")) +
  ggtitle("Metabolite tranfer") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))

dev.off()



# top different transfers
tri <- smetana %>% filter(community %in% c("tri-green-229","tri-green-yeast-229"))

tri[tri$donor=="Lactobacillus_farciminis_KCTC_3681_DSM_20184","donor"] = "L. farciminis"
tri[tri$receiver=="Lactobacillus_farciminis_KCTC_3681_DSM_20184","receiver"] = "L. farciminis"
tri[tri$donor=="Lactobacillus_johnsonii_DPC_6026","donor"] = "L. johnsonii"
tri[tri$receiver=="Lactobacillus_johnsonii_DPC_6026","receiver"] = "L. johnsonii"
tri[tri$donor=="yeast_smetana","donor"] = "Yeast"
tri[tri$receiver=="yeast_smetana","receiver"] = "Yeast"
tri[tri$donor=="Lactobacillus_brevis_ATCC_367","donor"] = "L. brevis"
tri[tri$receiver=="Lactobacillus_brevis_ATCC_367","receiver"] = "L. brevis"

diff_connection = vector()
unq_mets = unique(tri$compound_name)
for (i in 1:length(unq_mets)){
  met = unique(tri$compound_name)[i]
  tri_wo_yeast <- tri %>% filter(compound_name == met) %>% filter(community == "tri-green-229")
  tri_with_yeast <- tri %>% filter(compound_name == met) %>% filter(community == "tri-green-yeast-229")
  a <- nrow(tri_wo_yeast) / 3
  b <- nrow(tri_with_yeast) / 6
  diff_connection[i] = a-b
}

unq_met_stat = data.frame(unq_mets,diff_connection)
unq_met_stat = unq_met_stat[order(unq_met_stat$diff_connection),]

chosen_mets = c(head(unq_met_stat$unq_mets,5),tail(unq_met_stat$unq_mets,5))

tri <- tri %>% filter(compound_name %in% chosen_mets)


pdf("Documents/S_boulardii_modeling/results/Alex_20230711/top_transfer_v2.pdf")

ggplot(as.data.frame(tri),
       aes(y = smetana,
           axis1 = donor, axis2 = compound_name, axis3 = receiver)) +
  geom_alluvium(aes(fill = community),
                width = 1/8, knot.pos = 0, reverse = FALSE) +
  #scale_fill_manual(values = c(Brown = "#70493D", Hazel = "#E2AC76",
  #                             Green = "#3F752B", Blue = "#81B0E4")) +
  guides(fill = "none") +
  geom_stratum(alpha = .25, width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            reverse = FALSE) +
  scale_x_continuous(breaks = 1:3, labels = c("Donor", "Metabolite", "Receiver")) +
  ggtitle("Top metabolite tranfer") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))

dev.off()