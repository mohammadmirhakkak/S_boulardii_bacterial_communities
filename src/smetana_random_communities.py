#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 12:35:59 2023

@author: mohammadmirhakkak

Make dataframe of communities of three, four, and five members (random selection with and without fungal models) as input for
SMETANA
"""
import pandas as pd
import glob
import random
from itertools import combinations


model_dir = glob.glob("mohammadmirhakkak/S_boulardii_bacterial_communities/GEMs/random_models_smetana/*.xml")
model_file_names = [i.split('/')[-1][:-4] for i in model_dir]

fungal_model_file_names = [mdir for mdir in model_file_names if mdir.endswith('_smetana')]
bacterial_model_file_names = list(set(model_file_names) - set(fungal_model_file_names))



# repeat making the communities 20 times
random.seed(2024)
iteration_seeds = random.sample(range(1000),20)

for n in range(20):

    comm_names = []
    comm_names_with_fungi = []
    number = 0
    comm_members = []
    comm_members_with_fungi = []
    bacterial_combination = combinations(bacterial_model_file_names, 3)
    # make the community complete by adding every fungal model to every bacterial combination
    comm = []
    for bc in bacterial_combination:
        bc = list(bc)
        comm.append(bc + [fungal_model_file_names[0]])
        comm.append(bc + [fungal_model_file_names[1]])
        comm.append(bc + [fungal_model_file_names[2]])
    comm_number = len(comm)
    random.seed(iteration_seeds[n]) # set seed
    rand_index = random.sample(range(int(comm_number)),455) # 455 communities with 4 members (as many as yeast communities of 4 members in the paper)
    
    for j in rand_index:
        comm_iter = comm[j]
        comm_members = comm_members + comm_iter[:-1] # exclude the last one (fungal)
        comm_members_with_fungi = comm_members_with_fungi + comm_iter
        number += 1
        comm_names = comm_names + ['community' + str(number)] * len(comm_iter[:-1])
        comm_names_with_fungi = comm_names_with_fungi + ['community' + str(number) + '_with_fungi'] * len(comm_iter)
    
    df_comms_4 = pd.DataFrame({'1':comm_names,'2':comm_members})
    df_comms_4_with_fungi = pd.DataFrame({'1':comm_names_with_fungi,'2':comm_members_with_fungi})
    
    
    comm_names = []
    comm_names_with_fungi = []
    # number = 0; it should not reset to zero 5 member communities should have different IDs
    comm_members = []
    comm_members_with_fungi = []
    bacterial_combination = combinations(bacterial_model_file_names, 4)
    # make the community complete by adding every fungal model to every bacterial combination
    comm = []
    for bc in bacterial_combination:
        bc = list(bc)
        comm.append(bc + [fungal_model_file_names[0]])
        comm.append(bc + [fungal_model_file_names[1]])
        comm.append(bc + [fungal_model_file_names[2]])
    comm_number = len(comm)
    random.seed(iteration_seeds[n]) # set seed
    rand_index = random.sample(range(int(comm_number)),1365) # 1365 communities with 4 members (as many as yeast communities of 4 members in the paper)

    for j in rand_index:
        comm_iter = comm[j]
        comm_members = comm_members + comm_iter[:-1] # exclude the last one (fungal)
        comm_members_with_fungi = comm_members_with_fungi + comm_iter
        number += 1
        comm_names = comm_names + ['community' + str(number)] * len(comm_iter[:-1])
        comm_names_with_fungi = comm_names_with_fungi + ['community' + str(number) + '_with_fungi'] * len(comm_iter)
    
    df_comms_5 = pd.DataFrame({'1':comm_names,'2':comm_members})
    df_comms_5_with_fungi = pd.DataFrame({'1':comm_names_with_fungi,'2':comm_members_with_fungi})
    
    
    
    
    df_comms = pd.concat([df_comms_4,df_comms_5,df_comms_4_with_fungi,df_comms_5_with_fungi])
    
    df_comms.to_csv('mohammadmirhakkak/S_boulardii_bacterial_communities/res/random_communities_' + str(n+1) + '.tsv',
                       sep = '\t', header = False, index = False)
