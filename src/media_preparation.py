#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 14:16:14 2023

@author: mohammadmirhakkak

This file makes the media matrix for the pairwise analysis using CarveMe models
and makes a union of the media to make sure all the CarveMe models grow for
the sake of SMETANA analysis
"""


# function created for the yeast model to add lumen exchange reactions
# to simulate single growth (like in the main analysis code)
def join_lumen_rxns(model,mapping,Metabolite,Reaction):
    #Makes lumen exchange reactions and adds them to the model
    
    for i in range(len(model.reactions)):
        
        if len(model.reactions[i].metabolites)==1:

            model.reactions[i].lower_bound = 0

            met_id = list(model.reactions[i].metabolites)[0].id
            
            if met_id in mapping['yeast_REPLACEMENT_ID'].values:
                
                s = mapping.loc[mapping['yeast_REPLACEMENT_ID']==met_id,'CARVEME_abbr'].iloc[0]
                new_met = Metabolite(s+'_u',compartment = 'u')

                #add the lumen metabolite to the right side of the exchange reaction
                model.reactions[i].add_metabolites({new_met:1})

                model.reactions[i].bounds = (-1000,1000)

                #add an exchange reaction for the lumen metabolite (lumen <=> outside)
                r = Reaction("EX_u_"+s)
                model.add_reaction(r)
                r.add_metabolites({new_met:-1})
            
    return model


# change the yeast model for anaerobic growth
def anaerobicYeast(model):
    
    # change the coefficients of ATP, ADP, phosphate, H+ and water in pseudo
    # biomass reaction to anaerobic values (taken from the github of the 
    # original Yeast-GEM nature paper)
    model.reactions.r_4041.reaction = '30.49 s_0434 + 30.49 s_0803 + s_1096 + \
        s_3717 + s_3718 + s_3719 + s_3720 + s_4205 + s_4206 --> 30.49 s_0394 \
        + s_0450 + 30.49 s_0794 + 30.49 s_1322'
    
    # force NGAM (maintenance) reaction to not carry flux
    model.reactions.r_4046.bounds = (0,0)
    
    # heme a, NAD, NADH, NADP, NADPH, and CoA must be deleted from the
    # Cofactor pseudo reaction. They are not used under anaerobic condition
    # Only FAD, riboflavin, TDP, and THF are kept.
    model.reactions.r_4598.reaction = '9.99999974737875e-06 s_0687 + \
        0.000989999971352518 s_1405 + 1.20000004244503e-06 s_1475 + \
        6.34000025456771e-05 s_1487 --> s_4205'
    
    # Change media to anaerobic (no O2 uptake and allows sterol
    # and fatty acid exchanges)
    model.reactions.EX_u_o2.lower_bound = 0        #O2
    model.reactions.EX_u_ergst.lower_bound = -1000    #ergosterol
    model.reactions.EX_u_lanost.lower_bound = -1000    #lanosterol
    model.reactions.EX_u_hdcea.lower_bound = -1000    #palmitoleate
    model.reactions.EX_u_zymst.lower_bound = -1000    #zymosterol
    model.reactions.r_2134.lower_bound = -1000    #14-demethyllanosterol
    model.reactions.r_2137.lower_bound = -1000    #ergosta-5,7,22,24(28)-tetraen-3beta-ol
    model.reactions.EX_u_ocdcea.lower_bound = -1000    #oleate
    
    # Block pathways for proper glycerol production
    # Block oxaloacetate-malate shuttle (not present in anaerobic conditions)
    model.reactions.r_0713.lower_bound = 0 #Mithocondria
    model.reactions.r_0714.lower_bound = 0 #Cytoplasm
    #Block glycerol dehydroginase (only acts in microaerobic conditions)
    model.reactions.r_0487.upper_bound = 0
    #Block 2-oxoglutarate + L-glutamine -> 2 L-glutamate (alternative pathway)
    model.reactions.r_0472.upper_bound = 0
    
    return model

"""
Media Composition

SC media 
Component, amount (μg/L unless indicated)
Nitrogen Sources:
Ammonium sulfate, 5.0 g/L
 
Vitamins:
Biotin, 2.0
Calcium pantothenate, 400
Folic acid, 2.0
Inositol, 2.0 mg/L
Nicotinic acid, 400
p-Aminobenzoic acid, 200
Pyridoxine HCl, 400
Riboflavin, 200
Thiamine HCL, 400

SC media 
Component, amount (μg/L unless indicated)
Trace Elements:
Boric acid, 500
Copper sulfate, 40
Potassium iodide, 100
Ferric chloride, 200
Manganese sulfate, 400
Sodium molybdate, 200
Zinc sulfate, 400
 
Salts:
Potassium phosphate monobasic, 1.0 g/L
Magnesium sulfate, 0.5 g/L
Sodium chloride, 0.1 g/L
Calcium chloride, 0.1 g/L

SC media 
Component, amount (μg/L unless indicated)
Later added
Sodium acetate, 5 g/L 
Polysorbate 80, 1 mL/L 
Glucose, 20 g/L 
"""

import sys
sys.path
sys.path.append('/path/to/cplex/solver')

import pandas as pd
import glob
import cobra
from cobra import Metabolite,Reaction
from cobra.flux_analysis import flux_variability_analysis
import numpy as np

media_dic = dict()

# Nitrogen sources
# Ammonium sulfate -> SO4. 2NH4: 132.1395; separate ammonium and sulfate
media_dic['EX_so4_e'] = 5 * 1/132.1395 * (1000) * -1
media_dic['EX_nh4_e'] = 5 * 1/132.1395 * (1000) * -1

# Vitamins
media_dic['EX_btn_e'] = 2 * 1/244.3106 * (1/10**6) * (1000) * -1
# Calcium pantothenate -> (C9H16NO5)2. Ca: 476.5321; separate Calcium and pantothenate
media_dic['EX_ca2_e'] = 0.1 * 1/476.5321 * (1000) * -1
media_dic['EX_pnto__R_e'] = 400 * 1/476.5321 * (1/10**6) * (1000) * -1
media_dic['EX_fol_e'] = 2 * 1/441.3975 * (1/10**6) * (1000) * -1
media_dic['EX_inost_e'] = 2 * 1/180.1559 * (1/10**3) * (1000) * -1
media_dic['EX_nac_e'] = 400 * 1/123.1094 * (1/10**6) * (1000) * -1
media_dic['EX_anth_e'] = 200 * 1/137.136 * (1/10**6) * (1000) * -1
# Pyridoxine HCl -> (C9H16NO5)2. Ca: separate Pyridoxine and Cl
media_dic['EX_pydxn_e'] = 400 * 1/169.1778 * (1/10**6) * (1000) * -1
media_dic['EX_cl_e'] = 0.1 * 1/35.453 * (1000) * -1
media_dic['EX_ribflv_e'] = 200 * 1/376.3639 * (1/10**6) * (1000) * -1
media_dic['EX_thm_e'] = 400 * 1/265.3546 * (1/10**6) * (1000) * -1
media_dic['EX_thm_e'] = 400 * 1/265.3546 * (1/10**6) * (1000) * -1

# Trace elements
media_dic['EX_cu2_e'] = 40 * 1/63.546 * (1/10**6) * (1000) * -1
media_dic['EX_k_e'] = 1 * 1/39.0983 * (1000) * -1
media_dic['EX_i_e'] = 100 * 1/126.9045 * (1/10**6) * (1000) * -1
media_dic['EX_fe3_e'] = 200 * 1/55.845 * (1/10**6) * (1000) * -1
media_dic['EX_fe2_e'] = 200 * 1/55.845 * (1/10**6) * (1000) * -1
media_dic['EX_mn2_e'] = 400 * 1/54.938 * (1/10**6) * (1000) * -1
media_dic['EX_na1_e'] = 5 * 1/22.9898 * (1000) * -1
media_dic['EX_mobd_e'] = 200 * 1/161.9535 * (1/10**6) * (1000) * -1
media_dic['EX_zn2_e'] = 400 * 1/65.409 * (1/10**6) * (1000) * -1
media_dic['EX_pi_e'] =  1 * 1/97.9952 * (1000) * -1
media_dic['EX_mg2_e'] =  0.5 * 1/24.305 * (1000) * -1
media_dic['EX_ac_e'] =  5 * 1/60.052 * (1000) * -1
media_dic['EX_glc__D_e'] =  -10 * 20 / (5 + 20) # constrain glc uptake regarding
# Michaelis Menten equation

# Fatty acids are needed for Yeast's anaerobic growth
media_dic['EX_lanost_e'] = -1000
media_dic['EX_ergst_e'] = -1000
media_dic['EX_hdcea_e'] = -1000
media_dic['EX_zymst_e'] = -1000
media_dic['EX_ocdcea_e'] = -1000

# import the models and find the required additional nutrients to make them viable
# through FVA
model_dir = glob.glob("mohammadmirhakkak/S_boulardii_bacterial_communities/GEMs/*.xml")
model_dir = set(model_dir) - {'mohammadmirhakkak/S_boulardii_bacterial_communities/GEMs/yeast-GEM.xml',
                             'mohammadmirhakkak/S_boulardii_bacterial_communities/GEMs/yeast_smetana.xml'}

media = pd.DataFrame()
gr_rates = []
for i in model_dir:
    
    if i.endswith('.xml'):
        agora = cobra.io.read_sbml_model(i)
    elif i.endswith('.mat'):
        agora = cobra.io.load_matlab_model(i)
    gem_id = i.split('/')[-1][:-4]
    
    # constrain the GEMs using the defined diet and relax the other 
    # exchange bounds to identify which metabolites are required
    ex = []
    for i in agora.reactions:
        if i.id[:2]=='EX':
            if i.id in media_dic.keys():
                i.lower_bound = media_dic[i.id]
            else:
                ex.append(i)
                i.bounds = (-1000,1000)
    
    s = str(agora.objective.expression)
    loc1 = s.find('*')+1
    loc2 = s.find(' ')
    obj = s[loc1:loc2]
    gr = agora.slim_optimize()
    if gr >= 0.01:
        agora.reactions.get_by_id(obj).bounds = (0.01,0.01)
    else:
        agora.reactions.get_by_id(obj).bounds = (gr,gr)
                
    fva = flux_variability_analysis(agora,ex)
    fva = fva.query("minimum < -1e-6 and maximum < -1e-6")
    
    
    # check if still there are influxes required for the bacterial growth
    ex = [r for r in ex if r.id not in fva.index] # exclude identified reactions in FVA
    for r in ex:
        r.lower_bound = 0
    
    gr = agora.slim_optimize()
    if gr > 1e-6:
        table_common = pd.DataFrame({gem_id:media_dic.values()},index = media_dic.keys())
        table_essential1 = pd.DataFrame({gem_id:fva.maximum.values},index = fva.index)
        table_merged = pd.concat([table_common,table_essential1])
        media = pd.merge(media,table_merged, left_index=True, right_index=True,how = 'outer')
        continue
    #ex = random.sample(ex,len(ex))
    required = []
    counter = -1
    close = False
    limit = len(ex)
    while counter < (limit-1):
        
        counter+=1
        
        if not close:
            ex[counter].lower_bound = -0.2
            #print('open',counter,'until',limit)
        else:
            ex[counter].lower_bound = 0
            #print('close',counter,'until',limit)
            
        if agora.slim_optimize() > 1e-6 and not close:
            required.append(ex[counter])
            limit = counter
            counter = -1
            close = True
            #print('essential uptake is',counter)
            #print('growth rate is',agora.slim_optimize())
            continue
        elif (agora.slim_optimize() < 1e-6 or np.isnan(agora.slim_optimize())) and close:
            #print('closed', counter, 'uptake was essential. Open again')
            ex[counter].lower_bound = -0.2
            required.append(ex[counter])

    gr_rates.append(agora.slim_optimize())
    
    table_common = pd.DataFrame({gem_id:media_dic.values()},index = media_dic.keys())
    table_essential1 = pd.DataFrame({gem_id:fva.maximum.values},index = fva.index)
    table_essential2 = pd.DataFrame({gem_id:[-0.2]*len(required)},index = [r.id for r in required])
    table_essential = pd.concat([table_essential1,table_essential2])
    table_merged = pd.concat([table_common,table_essential])
    media = pd.merge(media,table_merged, left_index=True, right_index=True,how = 'outer')
    
    
    
            

media[media.isna()] = 0
# sort the dataframe; the common metabolites come first
# the essential metabolites come second
essen_comp = list(set(media.index) - set(media_dic.keys()))
common_comp = list(media_dic.keys())
media = media.loc[common_comp + essen_comp,]

# check if the yeast model is able to grow on the media
yeast = cobra.io.read_sbml_model("mohammadmirhakkak/S_boulardii_bacterial_communities/GEMs/yeast-GEM.xml")
id_mapping = pd.read_csv('mohammadmirhakkak/S_boulardii_bacterial_communities/dat/mets_map_v3_carveme.csv', index_col=0)
yeast = join_lumen_rxns(yeast,id_mapping,Metabolite,Reaction)
media_yeast_index = ['EX_u_' + i[3:-2] for i in media.index]
media_yeast = media.copy()
media_yeast.index = media_yeast_index
media_yeast = media_yeast.to_dict()

for i in media_yeast.keys():
    tested_media = media_yeast[i]
    
    for r in yeast.reactions:
        if len(r.metabolites)==1:
            if r.id in tested_media.keys():
                r.lower_bound = tested_media[r.id]
            else:
                r.lower_bound = 0
    
    yeast = anaerobicYeast(yeast)
    
    gr = yeast.slim_optimize()
    print(gr)
    
    '''
    if gr >= 0.01:
        yeast.reactions.get_by_id('r_2111').bounds = (0.01,0.01)
    else:
        yeast.reactions.get_by_id('r_2111').bounds = (gr,gr)
    
    fva = flux_variability_analysis(yeast,ex)
    fva = fva.query("minimum < -1e-7 and maximum < -1e-7")
    '''
                    
ind = [i[3:-2] for i in media.index]
media.index = ind
media.to_csv("mohammadmirhakkak/S_boulardii_bacterial_communities/res/media_carveme.csv")


# media for smetana

bigg = pd.read_table("mohammadmirhakkak/S_boulardii_bacterial_communities/dat/bigg_models_metabolites.txt")

compound_ex = [i+'_e' for i in media.index]
compound = [i for i in media.index]
c=-1
name = []
for i in compound_ex:
    c+=1
    try:
        name.append(bigg.loc[bigg.bigg_id == i,'name'].values[0])
    except:
        print(c)
        name.append('')
name[33] = 'Gly-Tyr'
name[97] = 'Gly-Cys'
name[107] = 'Gly-Asp'
    
smetana_media = pd.DataFrame({"medium":['SC']*len(name), "description":["SC+essentials"]*len(name), \
               "compound":compound, "name":name})
smetana_media.to_csv("mohammadmirhakkak/S_boulardii_bacterial_communities/res/smetana_media.tsv",sep = '\t', index = False)
