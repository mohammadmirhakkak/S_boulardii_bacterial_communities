#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 10:21:41 2024

@author: mohammadmirhakkak
"""
import cobra
import pandas as pd
import glob
import numpy as np

fungal_model_dir = glob.glob("Desktop/candida albicans coreco refinement/GSMMs used as db/*.xml")

# read the models one by one and adjust the exchange reaction and metabolite IDs



# Penicillium chrysogenum genome-scale model
iAL1006 = cobra.io.read_sbml_model("Desktop/candida albicans coreco refinement/GSMMs used as db/iAL1006.xml")

id_mapping = pd.read_csv('Documents/S_boulardii_modeling/data/20230206/mets_map_v3.csv', index_col=0)

bigg = pd.read_table("Documents/S_boulardii_modeling/data/20230206/bigg_models_reactions.txt")

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

id_mapping.index = [i.lower() for i in id_mapping.index]
metname2biggex = id_mapping.to_dict()['bigg_reaction']


for ex in iAL1006.exchanges:
    met_name = list(ex.metabolites)[0].name.lower()
    if met_name in metname2biggex.keys():
        mapped_ex_id = metname2biggex[met_name]
        if type(mapped_ex_id) is not str:
            continue
        ex.id = mapped_ex_id
        ex_met = list(ex.metabolites)[0]
        ex_met_id = mapped_ex_id[3:-2] + '_e'
        ex_met.id = ex_met_id

# check the unmapped exchanges for futher investigation and manual curation
for ex in iAL1006.exchanges:
    met_name = list(ex.metabolites)[0].name.lower()
    if met_name not in metname2biggex.keys():
        print(ex)
        print(met_name)
        print()

metname2biggex['(s)-lactate'] = 'EX_lac__L_e'
metname2biggex['1,3-beta-d-glucan'] = 'EX_13BDglcn_e'
metname2biggex['alpha,alpha-trehalose'] = 'EX_tre_e'
metname2biggex['alpha-l-arabinan'] = 'EX_araban__L_e'
metname2biggex['aminoacetaldehyde'] = 'EX_aacald_e'
metname2biggex['anthranilate'] = 'EX_anth_e'
metname2biggex['beta-alanine'] = 'EX_ala_B_e'
metname2biggex['cellobiose'] = 'EX_cellb_e'
metname2biggex['chitin'] = 'EX_chtn_e'
metname2biggex['cellulose'] = 'EX_cell4_e'
metname2biggex['chitobiose'] = 'EX_chitob_e'
metname2biggex['chitosan'] = 'EX_chitos_e'
metname2biggex['co2'] = 'EX_co2_e'
metname2biggex['cyanate'] = 'EX_cynt_e'
metname2biggex['d-arabinitol'] = 'EX_abt__D_e'
metname2biggex['decanoate'] = 'EX_dca_e'
metname2biggex['d-gluconate'] = 'EX_glcn_e'
metname2biggex['d-mannitol'] = 'EX_mnl_e'
metname2biggex['glutathione'] = 'EX_gthox_e'
metname2biggex['glycine betaine'] = 'EX_glyb_e'
metname2biggex['glycogen'] = 'EX_glygn4_e'
metname2biggex['h2o'] = 'EX_h2o_e'
metname2biggex['h2o2'] = 'EX_h2o2_e'
metname2biggex['h2s'] = 'EX_h2s_e'
metname2biggex['heptadecanoate'] = 'EX_hpdca_e'
metname2biggex['heptadecenoate'] = 'EX_M00003_e'
metname2biggex['heptanoate'] = 'EX_hpta_e'
metname2biggex['hexadecenoate'] = 'EX_hdcea_e'
metname2biggex['hexanoate'] = 'EX_hxa_e'
metname2biggex['icosanoate'] = 'EX_arach_e'
metname2biggex['isocitrate'] = 'EX_icit_e'
metname2biggex['l-2-aminoadipate'] = 'EX_L2aadp_e'
metname2biggex['lactose'] = 'EX_lcts_e'
metname2biggex['laurate'] = 'EX_ddca_e'
metname2biggex['l-cystine'] = 'EX_cysi__L_e'
metname2biggex['l-homocysteine'] = 'EX_hcys__L_e'
metname2biggex['l-ornithine'] = 'EX_orn_e'
metname2biggex['l-ribulose'] = 'EX_rbl__L_e'
metname2biggex['maltotriose'] = 'EX_malttr_e'
metname2biggex['mannans'] = 'EX_mannan_e'
metname2biggex['n-acetyl-d-glucosamine'] = 'EX_acgam_e'
metname2biggex['nh3'] = 'EX_nh3_e'
metname2biggex['nicotinamide'] = 'EX_ncam_e'
metname2biggex['nitrate'] = 'EX_no3_e'
metname2biggex['nitrite'] = 'EX_no2_e'
metname2biggex['nonanoate'] = 'EX_nona_e'
metname2biggex['o2'] = 'EX_o2_e'
metname2biggex['octadecadienoate'] = 'EX_ocdcya_e'
metname2biggex['octadecenoate'] = 'EX_ocdcea_e'
metname2biggex['oxalate'] = 'EX_oxa_e'
metname2biggex['pentadecanoate'] = 'EX_ptdca_e'
metname2biggex['phenylacetate'] = 'EX_pac_e'
metname2biggex['pimelate'] = 'EX_pime_e'
metname2biggex['quinolinate'] = 'EX_quln_e'
metname2biggex['starch'] = 'EX_strch1_e'
metname2biggex['stearate'] = 'EX_ocdca_e'
metname2biggex['sulfite'] = 'EX_so3_e'
metname2biggex['sulfur'] = 'EX_so2_e'
metname2biggex['xylans'] = 'EX_xylan4_e'


for ex in iAL1006.exchanges:
    met_name = list(ex.metabolites)[0].name.lower()
    if met_name in metname2biggex.keys():
        mapped_ex_id = metname2biggex[met_name]
        if type(mapped_ex_id) is not str:
            continue
        ex.id = mapped_ex_id
        ex_met = list(ex.metabolites)[0]
        ex_met_id = mapped_ex_id[3:-2] + '_e'
        ex_met.id = ex_met_id

iAL1006_unmapped = []
for ex in iAL1006.exchanges:
    met_name = list(ex.metabolites)[0].name.lower()
    if met_name not in metname2biggex.keys():
        iAL1006_unmapped.append(met_name)





# Aspergillus niger genome-scale model
iMA871 = cobra.io.read_sbml_model("Desktop/candida albicans coreco refinement/GSMMs used as db/iMA871.xml")

for ex in iMA871.exchanges:
    met_name = list(ex.metabolites)[0].name.lower()
    if met_name in metname2biggex.keys():
        mapped_ex_id = metname2biggex[met_name]
        if type(mapped_ex_id) is not str:
            continue
        if mapped_ex_id in iMA871.reactions:
            continue
        ex.id = mapped_ex_id
        ex_met = list(ex.metabolites)[0]
        ex_met_id = mapped_ex_id[3:-2] + '_e'
        ex_met.id = ex_met_id

# check the unmapped exchanges for futher investigation and manual curation
for ex in iMA871.exchanges:
    met_name = list(ex.metabolites)[0].name.lower()
    if met_name not in metname2biggex.keys() and met_name not in iAL1006_unmapped:
        print(ex)
        print(met_name)
        print()

metname2biggex['orthophosphate'] = 'EX_pi_e'
metname2biggex['d-ribulose'] = 'EX_rbl__D_e'
metname2biggex['d-xylulose'] = 'EX_xylu__D_e'
metname2biggex['l-xylulose'] = 'EX_xylu__L_e'
metname2biggex['mannan'] = 'EX_mannan_e'
metname2biggex['arabinan'] = 'EX_araban__L_e'
metname2biggex['d-arabitol'] = 'EX_abt__D_e'
metname2biggex['l-arabitol'] = 'EX_abt__L_e'
metname2biggex['ribitol'] = 'EX_rbt_e'
metname2biggex['galactitol'] = 'EX_galt_e'
metname2biggex['d-sorbitol'] = 'EX_sbt__D_e'
metname2biggex['cis-aconitate'] = 'EX_acon_C_e'
metname2biggex['d-galactonate'] = 'EX_galctn__D_e'
metname2biggex['d-glucuronate'] = 'EX_glcur_e'
metname2biggex['d-lactate'] = 'EX_lac__D_e'
metname2biggex['tartrate'] = 'EX_tartr__D_e'
metname2biggex['4-aminobutanoate'] = 'EX_4abut_e'
metname2biggex['quinate'] = 'EX_quin_e'
metname2biggex['shikimate'] = 'EX_skm_e'
metname2biggex['butanoate'] = 'EX_but_e'
metname2biggex['dodecanoate'] = 'EX_ddca_e'
metname2biggex['tetradecanoate'] = 'EX_ttdca_e'
metname2biggex['tetradecenoate'] = 'EX_ttdcea_e'
metname2biggex['hexadecanoate'] = 'EX_hdca_e'
metname2biggex['octadecanoate'] = 'EX_ocdca_e'
metname2biggex['nonadecanoate'] = 'EX_M02613_e'
metname2biggex['eicosanoate'] = 'EX_arach_e'
metname2biggex['benzoic acid'] = 'EX_bz_e'
metname2biggex['coumarate'] = 'EX_2coum_e'
metname2biggex['ferulic acid'] = 'EX_fer_e'
metname2biggex['vanillate'] = 'EX_vanlt_e'
metname2biggex['indole'] = 'EX_indole_e'
metname2biggex['hydroxybenzoate'] = 'EX_3hbz_e'
metname2biggex['salicylate'] = 'EX_salc_e'
metname2biggex['gallic acid'] = 'EX_ga_e'
metname2biggex['phenylpyruvate'] = 'EX_phpyr_e'
metname2biggex['atp'] = 'EX_atp_e'
metname2biggex['l-kynurenine'] = 'EX_Lkynr_e'
metname2biggex['dethiobiotin'] = 'EX_dtbt_e'
metname2biggex['glycerone'] = 'EX_dha_e'


for ex in iMA871.exchanges:
    met_name = list(ex.metabolites)[0].name.lower()
    if met_name in metname2biggex.keys():
        mapped_ex_id = metname2biggex[met_name]
        if type(mapped_ex_id) is not str:
            continue
        if mapped_ex_id in iMA871.reactions:
            continue
        ex.id = mapped_ex_id
        ex_met = list(ex.metabolites)[0]
        ex_met_id = mapped_ex_id[3:-2] + '_e'
        ex_met.id = ex_met_id

iMA871_unmapped = []
for ex in iMA871.exchanges:
    met_name = list(ex.metabolites)[0].name.lower()
    if met_name not in metname2biggex.keys():
        iMA871_unmapped.append(met_name)




# Aspergillus oryzae genome-scale model
iWV1314 = cobra.io.read_sbml_model("Desktop/candida albicans coreco refinement/GSMMs used as db/iWV1314.xml")

# remove formula from met name
for ex in iWV1314.exchanges:
    met_name = list(ex.metabolites)[0].name.lower()
    list(ex.metabolites)[0].name = met_name.split('_')[0]

for ex in iWV1314.exchanges:
    met_name = list(ex.metabolites)[0].name.lower()
    if met_name in metname2biggex.keys():
        mapped_ex_id = metname2biggex[met_name]
        if type(mapped_ex_id) is not str:
            continue
        if mapped_ex_id in iWV1314.reactions:
            continue
        ex.id = mapped_ex_id
        ex_met = list(ex.metabolites)[0]
        ex_met_id = mapped_ex_id[3:-2] + '_e'
        ex_met.id = ex_met_id
    
# check the unmapped exchanges for futher investigation and manual curation
for ex in iWV1314.exchanges:
    met_name = list(ex.metabolites)[0].name.lower()
    if met_name not in metname2biggex.keys() and met_name not in iAL1006_unmapped and met_name not in iMA871_unmapped:
        print(ex)
        print(met_name)
        print()


iWV1314_unmapped = []
for ex in iWV1314.exchanges:
    met_name = list(ex.metabolites)[0].name.lower()
    if met_name not in metname2biggex.keys():
        iWV1314_unmapped.append(met_name)
        

# make common exchange IDs for common exchanges among fungal models
common1 = set(iWV1314_unmapped) & set(iMA871_unmapped)
for ex in iWV1314.exchanges:
    met_name = list(ex.metabolites)[0].name.lower()
    if met_name in common1:
        print(ex)
ex = iWV1314.reactions.r2217
ex.id = 'EX_Ee_e'
ex = iWV1314.reactions.r2220
ex.id = 'EX_PTATEe_e'
ex = iWV1314.reactions.r2221
ex.id = 'EX_CELLOTe_e'
ex = iWV1314.reactions.r2224
ex.id = 'EX_AMYLPe_e'
ex = iWV1314.reactions.r2225
ex.id = 'EX_AMYLSe_e'
ex = iWV1314.reactions.r2235
ex.id = 'EX_EOLe_e'
for ex in iMA871.exchanges:
    met_name = list(ex.metabolites)[0].name.lower()
    if met_name in common1:
        print(ex)
ex = iMA871.reactions.get_by_id('1236')
ex.id = 'EX_Ee_e'
ex = iMA871.reactions.get_by_id('1245')
ex.id = 'EX_PTATEe_e'
ex_met = iMA871.metabolites.get_by_id('PECTATEe')
ex_met.id = 'PTATEe'
ex = iMA871.reactions.get_by_id('1246')
ex.id = 'EX_CELLOTe_e'
ex = iMA871.reactions.get_by_id('1249')
ex.id = 'EX_AMYLPe_e'
ex = iMA871.reactions.get_by_id('1250')
ex.id = 'EX_AMYLSe_e'
ex = iMA871.reactions.get_by_id('1261')
ex.id = 'EX_EOLe_e'

common2 = set(iWV1314_unmapped) & set(iAL1006_unmapped)

common3 = set(iMA871_unmapped) & set(iAL1006_unmapped)
for ex in iMA871.exchanges:
    met_name = list(ex.metabolites)[0].name.lower()
    if met_name in common3:
        print(ex)
ex = iMA871.reactions.get_by_id('1227')
ex.id = 'EX_IDOLe_e'
ex = iMA871.reactions.get_by_id('1257')
ex.id = 'EX_GLCN15LACe_e'
ex = iMA871.reactions.get_by_id('1303')
ex.id = 'EX_C162e_e'
ex = iMA871.reactions.get_by_id('1309')
ex.id = 'EX_C183e_e'
for ex in iAL1006.exchanges:
    met_name = list(ex.metabolites)[0].name.lower()
    if met_name in common3:
        print(ex)
ex = iAL1006.reactions.get_by_id('EX_GLCN15LACb')
ex.id = 'EX_GLCN15LACe_e'
ex = iAL1006.reactions.get_by_id('EX_C162b')
ex.id = 'EX_C162e_e'
ex = iAL1006.reactions.get_by_id('EX_IDOLb')
ex.id = 'EX_IDOLe_e'
ex = iAL1006.reactions.get_by_id('EX_C183b')
ex.id = 'EX_C183e_e'
iAL1006.reactions.EX_GLCN15LACe_e.reaction = 'GLCN15LACe <=>'
iAL1006.reactions.EX_C162e_e.reaction = 'C162e <=>'
iAL1006.reactions.EX_IDOLe_e.reaction = 'IDOLe <=>'
iAL1006.reactions.EX_C183e_e.reaction = 'C183e <=>'

iAL1006.objective = 'bmOUT'

# save the modified models
cobra.io.write_sbml_model(iAL1006,'Documents/S_boulardii_modeling/models/ISME_revision/random_models_smetana/iAL1006_smetana.xml')
cobra.io.write_sbml_model(iMA871,'Documents/S_boulardii_modeling/models/ISME_revision/random_models_smetana/iMA871_smetana.xml')
cobra.io.write_sbml_model(iWV1314,'Documents/S_boulardii_modeling/models/ISME_revision/random_models_smetana/iWV1314_smetana.xml')
