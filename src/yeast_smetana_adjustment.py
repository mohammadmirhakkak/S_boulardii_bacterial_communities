#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 09:59:40 2023

@author: mohammadmirhakkak

This code changes the exchange reaction IDs of the Yeast model to BiGG IDs
for SMATANA simulations
"""

import sys
sys.path
sys.path.append('/path/to/cplex/solver/')

import cobra
import pandas as pd
import numpy as np

yeast = cobra.io.read_sbml_model("mohammadmirhakkak/S_boulardii_bacterial_communities/GEMs/yeast-GEM.xml")

id_mapping = pd.read_csv('mohammadmirhakkak/S_boulardii_bacterial_communities/dat/mets_map_v3.csv', index_col=0)

bigg = pd.read_table("mohammadmirhakkak/S_boulardii_bacterial_communities/dat/bigg_models_reactions.txt")

bigg_rxn_id = []
for i in id_mapping.AGORA_abbr.values:
    if not pd.isna(i):
        created_rxn_id = 'EX_' + i.replace('_','__') + '_e'
        if created_rxn_id in bigg.bigg_id.values:
            bigg_rxn_id.append(created_rxn_id)
        else:
            bigg_rxn_id.append(np.nan)
    else:
        bigg_rxn_id.append(np.nan)
id_mapping['bigg_reaction'] = bigg_rxn_id

# check if the remainder can be mapped by metabolite names
bigg_ex = bigg.loc[bigg.bigg_id.str.startswith('EX_'),]
ex_met_names = []
for i in bigg_ex.name.values:
    if not pd.isna(i):
        ex_met_names.append(i.lower().replace(' exchange',''))
    else:
        ex_met_names.append(np.nan)
bigg_ex['metabolite_names_lower'] = ex_met_names

for i in range(id_mapping.shape[0]):
    if pd.isna(id_mapping.bigg_reaction[i]):
        met_name = id_mapping.index[i].lower()
        if met_name in bigg_ex.metabolite_names_lower.values:
            print(bigg_ex.bigg_id[bigg_ex.metabolite_names_lower==met_name].values[0])
            print(i)
            id_mapping.iloc[i,-1] = bigg_ex.bigg_id[bigg_ex.metabolite_names_lower==met_name].values[0]

# manual mapping
id_mapping.loc['D-tagatose','bigg_reaction'] = 'EX_tag__D_e'
id_mapping.loc['(R,R)-2,3-butanediol','bigg_reaction'] = 'EX_btd_RR_e'
id_mapping.loc['1-(sn-glycero-3-phospho)-1D-myo-inositol','bigg_reaction'] = 'EX_g3pi_e'
id_mapping.loc['Gly-Asn','bigg_reaction'] = 'EX_gly_asn__L_e'
id_mapping.loc['Ala-Glu','bigg_reaction'] = 'EX_ala_L_glu__L_e'
id_mapping.loc['Ala-Asp','bigg_reaction'] = 'EX_ala_L_asp__L_e'
id_mapping.loc['Met-Ala','bigg_reaction'] = 'EX_met_L_ala__L_e'

# for the exchanged metabolites not present in the table, check if they can be 
# mapped to BiGG IDs
met_names_add = []
bigg_reactions_add = []
met_ids_add = []
for i in yeast.reactions:
    if len(i.metabolites)==1:
        met = list(i.metabolites)[0]
        if met.id not in id_mapping.yeast_REPLACEMENT_ID.values:
            met_name = met.name.lower()
            if met_name in bigg_ex.metabolite_names_lower.values:
                #print(bigg_ex.bigg_id[bigg_ex.metabolite_names_lower==met_name].values[0])
                #print(i)
                met_names_add.append(met_name)
                met_ids_add.append(met.id)
                bigg_reactions_add.append(bigg_ex.bigg_id[bigg_ex.metabolite_names_lower==met_name].values[0])

table_add = pd.DataFrame({'yeast_ID':[np.nan]*len(met_names_add), \
                          'yeast_REPLACEMENT_ID':met_ids_add,\
                          'AGORA_abbr':[np.nan]*len(met_names_add), \
                          'bigg_reaction':bigg_reactions_add},index = met_names_add)

id_mapping = pd.concat([id_mapping,table_add])      


for i in yeast.reactions:
    if len(i.metabolites)==1:
        met = list(i.metabolites)[0]
        if met.id not in id_mapping.yeast_REPLACEMENT_ID.values:
            print(met.id)
            print(met.name)
            
# manual
met_names_add = []
bigg_reactions_add = []
met_ids_add = []

met_names_add.append('(1->3)-beta-D-glucan')
bigg_reactions_add.append('EX_13BDglcn_e')
met_ids_add.append('s_0003')

met_names_add.append('(R)-mevalonate')
bigg_reactions_add.append('EX_mev__R_e')
met_ids_add.append('s_0029')

met_names_add.append('(S)-3-methyl-2-oxopentanoate')
bigg_reactions_add.append('EX_3mob_e')
met_ids_add.append('s_0058')

met_names_add.append("2'-deoxyadenosine")
bigg_reactions_add.append('EX_dad_2_e')
met_ids_add.append('s_0133')

met_names_add.append("2'-deoxyuridine")
bigg_reactions_add.append('EX_duri_e')
met_ids_add.append('s_0139')

met_names_add.append('2-methylbutanal')
bigg_reactions_add.append('EX_2mbald_e')
met_ids_add.append('s_0167')

met_names_add.append('2-methylbutanol')
bigg_reactions_add.append('EX_2mbtoh_e')
met_ids_add.append('s_0170')

met_names_add.append('2-methylbutyl acetate')
bigg_reactions_add.append('EX_2mbac_e')
met_ids_add.append('s_0173')

met_names_add.append('3-methylbutanal')
bigg_reactions_add.append('EX_3mbald_e')
met_ids_add.append('s_0235')

met_names_add.append('5-formyltetrahydrofolic acid')
bigg_reactions_add.append('EX_5fthf_e')
met_ids_add.append('s_0320')

met_names_add.append("adenosine 3',5'-bismonophosphate")
bigg_reactions_add.append('EX_pap_e')
met_ids_add.append('s_0391')

met_names_add.append('D-glucitol')
bigg_reactions_add.append('EX_sbt__D_e')
met_ids_add.append('s_0562')

met_names_add.append('decanoate')
bigg_reactions_add.append('EX_dca_e')
met_ids_add.append('s_0597')

met_names_add.append('isoamylol')
bigg_reactions_add.append('EX_iamoh_e')
met_ids_add.append('s_0930')

met_names_add.append('isobutanol')
bigg_reactions_add.append('EX_ibtol_e')
met_ids_add.append('s_0933')

met_names_add.append("N,N'-diformyldityrosine")
bigg_reactions_add.append('EX_Nbfortyr_e')
met_ids_add.append('s_1186')

met_names_add.append('tryptophol')
bigg_reactions_add.append('EX_ind3eth_e')
met_ids_add.append('s_1530')

met_names_add.append('hexanoate')
bigg_reactions_add.append('EX_hxcoa_e')
met_ids_add.append('s_2824')

met_names_add.append('alpha-maltotriose')
bigg_reactions_add.append('EX_malttr_e')
met_ids_add.append('s_4140')

met_names_add.append('methyl alpha-D-glucopyranoside')
bigg_reactions_add.append('EX_madg_e')
met_ids_add.append('s_4125')

met_names_add.append('glycerone')
bigg_reactions_add.append('EX_dha_e')
met_ids_add.append('s_4167')

met_names_add.append('N(omega)-phospho-L-arginine')
bigg_reactions_add.append('EX_argp_e')
met_ids_add.append('s_4093')

met_names_add.append('2-hydroxyethane-1-sulfonate')
bigg_reactions_add.append('EX_isetac_e')
met_ids_add.append('s_4116')

met_names_add.append('N-acetyl-L-glutamate')
bigg_reactions_add.append('EX_acglu_e')
met_ids_add.append('s_4172')


table_add = pd.DataFrame({'yeast_ID':[np.nan]*len(met_names_add), \
                          'yeast_REPLACEMENT_ID':met_ids_add,\
                          'AGORA_abbr':[np.nan]*len(met_names_add), \
                          'bigg_reaction':bigg_reactions_add},index = met_names_add)

id_mapping = pd.concat([id_mapping,table_add])

# The id_mapping table is the most complete map of yeast metabolites IDs to
# corresponding exchange reactions IDs of BiGG
# make a yeast model for SMETANA

for i in range(len(yeast.reactions)):
    r = yeast.reactions[i]
    if len(r.metabolites)==1:
        met = list(r.metabolites)[0]
        if met.id in id_mapping.yeast_REPLACEMENT_ID.values and not pd.isna(id_mapping.loc[id_mapping.yeast_REPLACEMENT_ID == met.id,'bigg_reaction'].values[0]):
            mapped_id = id_mapping.loc[id_mapping.yeast_REPLACEMENT_ID == met.id,'bigg_reaction'].values[0]
            r.id = mapped_id
            met.id = mapped_id[3:-2]
            
r = yeast.reactions.r_2111
r.id = 'biomass'

bigg_metabolite = []
for i in id_mapping.bigg_reaction.values:
    if not pd.isna(i):
        bigg_metabolite.append(i[3:-2])
    else:
        bigg_metabolite.append(i)

id_mapping['CARVEME_abbr'] = bigg_metabolite


# final check to see if there are metabolites having yeast metabolite ID, but
# not having bigg IDs
for i in range(id_mapping.shape[0]):
    met_id = id_mapping.iloc[i,1]
    bigg_met_id = id_mapping.iloc[i,4]
    if not pd.isna(met_id) and pd.isna(bigg_met_id):
        print(id_mapping.index[i])
        print(id_mapping.AGORA_abbr[i])
        print(i)
        print()

id_mapping.iloc[154,3] = 'EX_lpchol_hs_e'
id_mapping.iloc[154,4] = 'lpchol_hs'

id_mapping.iloc[185,3] = 'EX_gly_asn__L_e'
id_mapping.iloc[185,4] = 'gly_asn__L'

id_mapping = id_mapping.drop('triphosphate',axis=0)
id_mapping = id_mapping.drop('carbamoyl phosphate',axis=0)
id_mapping = id_mapping.drop('hexadecanal',axis=0)
id_mapping = id_mapping.drop('3-sulfino-L-alanine',axis=0)

id_mapping.to_csv("S_boulardii_modeling/data/20230206/mets_map_v3_carveme.csv")


# switch the yeast GEM to anaerobic mode #

# change the coefficients of ATP, ADP, phosphate, H+ and water in pseudo
# biomass reaction to anaerobic values (taken from the github of the 
# original Yeast-GEM nature paper)
yeast.reactions.r_4041.reaction = '30.49 s_0434 + 30.49 s_0803 + s_1096 + \
    s_3717 + s_3718 + s_3719 + s_3720 + s_4205 + s_4206 --> 30.49 s_0394 \
    + s_0450 + 30.49 s_0794 + 30.49 s_1322'

# force NGAM (maintenance) reaction to not carry flux
yeast.reactions.r_4046.bounds = (0,0)

# heme a, NAD, NADH, NADP, NADPH, and CoA must be deleted from the
# Cofactor pseudo reaction. They are not used under anaerobic condition
# Only FAD, riboflavin, TDP, and THF are kept.
yeast.reactions.r_4598.reaction = '9.99999974737875e-06 s_0687 + \
    0.000989999971352518 s_1405 + 1.20000004244503e-06 s_1475 + \
    6.34000025456771e-05 s_1487 --> s_4205'

# Change media to anaerobic (no O2 uptake and allows sterol
# and fatty acid exchanges)
yeast.reactions.EX_o2_e.lower_bound = -0.2        #O2
yeast.reactions.EX_ergst_e.lower_bound = -1000    #ergosterol
yeast.reactions.EX_lanost_e.lower_bound = -1000    #lanosterol
yeast.reactions.EX_hdcea_e.lower_bound = -1000    #palmitoleate
yeast.reactions.EX_zymst_e.lower_bound = -1000    #zymosterol
yeast.reactions.r_2134.lower_bound = -1000    #14-demethyllanosterol
yeast.reactions.r_2137.lower_bound = -1000    #ergosta-5,7,22,24(28)-tetraen-3beta-ol
yeast.reactions.EX_ocdcea_e.lower_bound = -1000    #oleate

# Block pathways for proper glycerol production
# Block oxaloacetate-malate shuttle (not present in anaerobic conditions)
yeast.reactions.r_0713.lower_bound = 0 #Mithocondria
yeast.reactions.r_0714.lower_bound = 0 #Cytoplasm
#Block glycerol dehydroginase (only acts in microaerobic conditions)
yeast.reactions.r_0487.upper_bound = 0
#Block 2-oxoglutarate + L-glutamine -> 2 L-glutamate (alternative pathway)
yeast.reactions.r_0472.upper_bound = 0

# after saving the model, kineticLaw must be added to the biomass reaction in
# the sbml file
cobra.io.write_sbml_model(yeast,"mohammadmirhakkak/S_boulardii_bacterial_communities/GEMs/yeast_smetana.xml")
