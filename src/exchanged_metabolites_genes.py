#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 13:27:52 2024

@author: mohammadmirhakkak
"""
import cobra
import pandas as pd
import glob
from cobra.flux_analysis import flux_variability_analysis

model_dir = glob.glob("mohammadmirhakkak/S_boulardii_bacterial_communities/models/*.xml")


gr_rules = []
model_ids = []
exchanged_mets = []
reaction_names = []
reaction_ids = []
compartments = []
for mdir in model_dir:
    
    #exclude original yeast model
    if 'yeast-GEM' in mdir:
        continue
    
    model = cobra.io.read_sbml_model(mdir)
    
    for r in model.exchanges:
        
        ex_met = list(r.metabolites)[0]
        ex_met_name = ex_met.name
        
        # find the metabolite in different compartments
        for met in model.metabolites:
            if met.name == ex_met_name:
                for ex_met_rxn in met.reactions:
                    gr_rules.append(ex_met_rxn.gene_reaction_rule)
                    model_ids.append(model.id)
                    compartments.append(model.compartments[met.compartment])
                    exchanged_mets.append(ex_met.name)
                    reaction_names.append(ex_met_rxn.name)
                    reaction_ids.append(ex_met_rxn.id)

ex_met_genes = pd.DataFrame({'GSMM ID':model_ids,
                             'Exchanged metabolite':exchanged_mets,
                             'Compartment':compartments,
                             'Metabolizing reaction ID':reaction_ids,
                             'Metabolizing reaction name':reaction_names,
                             'Gene-Reaction rule':gr_rules})


fva_all = pd.DataFrame()
for mdir in model_dir:
    
    #exclude original yeast model
    if 'yeast-GEM' in mdir:
        continue
    
    model = cobra.io.read_sbml_model(mdir)
    
    for ex in model.exchanges:
        ex.bounds = (-1000,1000)
    
    exchanges = [r for r in model.exchanges]
    ex2name = {r.id:(list(r.metabolites)[0].name+' exchange') for r in model.exchanges}
    
    fva = flux_variability_analysis(model,exchanges,fraction_of_optimum=0.95,pfba_factor = 1.1)
    
    fva['GSMM ID'] = [model.id]*fva.shape[0]
    
    ex_names = []
    for ex in fva.index:
        ex_names.append(ex2name[ex])
    
    fva['Reaction name'] = ex_names
    
    fva['Reaction ID'] = fva.index
    
    fva_all = pd.concat([fva_all,fva])
    
# sort the fva columns
fva_all = fva_all.loc[:,['GSMM ID','Reaction ID','Reaction name','minimum','maximum']]
fva_all.columns = ['GSMM ID','Reaction ID','Reaction name','Minimum flux','Maximum flux']
fva_all.index = range(fva_all.shape[0])

with pd.ExcelWriter("mohammadmirhakkak/S_boulardii_bacterial_communities/res/exchanged_metabolites.xlsx") as writer:
    ex_met_genes.to_excel(writer,sheet_name='Gene Association',index=False)
    fva_all.to_excel(writer,sheet_name = 'Flux Variability Analysis',index=False)
