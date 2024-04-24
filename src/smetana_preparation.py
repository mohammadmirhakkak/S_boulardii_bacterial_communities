#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 10:47:45 2023

@author: mohammadmirhakkak

This code makes a dataframe in which GEMs are grouped in communities. It will
be the input for SMETANA calculations
"""

import pandas as pd
import glob
import itertools

red = {"Documents/S_boulardii_modeling/models/20230206/smetana/L_gasseri.xml"}

yellow = {"Documents/S_boulardii_modeling/models/20230206/smetana/Lactococcus_lactis_subsp_lactis_Il1403_IL1403.xml"}

all_gems = set(glob.glob("Documents/S_boulardii_modeling/models/20230206/smetana/*.xml")) \
    - {'Documents/S_boulardii_modeling/models/20230206/smetana/yeast_smetana.xml'}
green = all_gems - red - yellow

bug2groups = dict()
bug2groups['red'] = ['L_gasseri']
bug2groups['yellow'] = ['Lactococcus_lactis_subsp_lactis_Il1403_IL1403']
bug2groups['green'] = [i.split('/')[-1][:-4] for i in all_gems]


smetana_gem = list()
smetana_comm = list()


all_tri = itertools.combinations(list(all_gems), 3)

c=1
for tri in all_tri:
    a=[]
    for one_dir in tri:
        bug = one_dir.split('/')[-1][:-4]
        smetana_gem.append(bug)
        if bug in bug2groups['red']:
            a.append('red')
        elif bug in bug2groups['yellow']:
            a.append('yellow')
        elif bug in bug2groups['green']:
            a.append('green')
    a.sort()
    a = set(a)
    sign = ''
    for i in a:
        sign = sign + i + '-'
    
    smetana_comm = smetana_comm + ['tri-'+sign+str(c)]*3
    
    c+=1
    


# three bacteria with boulardii
all_tri = itertools.combinations(list(all_gems), 3)

c=1
for tri in all_tri:
    a=[]
    for one_dir in tri:
        bug = one_dir.split('/')[-1][:-4]
        smetana_gem.append(bug)
        if bug in bug2groups['red']:
            a.append('red')
        elif bug in bug2groups['yellow']:
            a.append('yellow')
        elif bug in bug2groups['green']:
            a.append('green')
    a.sort()
    a = set(a)
    sign = ''
    for i in a:
        sign = sign + i + '-'
    
    smetana_gem.append('yeast_smetana')
    sign = sign + 'yeast-'
    
    smetana_comm = smetana_comm + ['tri-'+sign+str(c)]*4
    
    c+=1
    
    
    
    
all_four = itertools.combinations(list(all_gems), 4)

c=1
for four in all_four:
    a=[]
    for one_dir in four:
        bug = one_dir.split('/')[-1][:-4]
        smetana_gem.append(bug)
        if bug in bug2groups['red']:
            a.append('red')
        elif bug in bug2groups['yellow']:
            a.append('yellow')
        elif bug in bug2groups['green']:
            a.append('green')
    a.sort()
    a = set(a)
    sign = ''
    for i in a:
        sign = sign + i + '-'
    
    smetana_comm = smetana_comm + ['quad-'+sign+str(c)]*4
    
    c+=1
    


# four bacteria with boulardii
all_four = itertools.combinations(list(all_gems), 4)

c=1
for four in all_four:
    a=[]
    for one_dir in four:
        bug = one_dir.split('/')[-1][:-4]
        smetana_gem.append(bug)
        if bug in bug2groups['red']:
            a.append('red')
        elif bug in bug2groups['yellow']:
            a.append('yellow')
        elif bug in bug2groups['green']:
            a.append('green')
    a.sort()
    a = set(a)
    sign = ''
    for i in a:
        sign = sign + i + '-'
    
    smetana_gem.append('yeast_smetana')
    sign = sign + 'yeast-'
    
    smetana_comm = smetana_comm + ['quad-'+sign+str(c)]*5
    
    c+=1




smetana_table = pd.DataFrame({'community id':smetana_comm,'organism id':smetana_gem})

#smetana_table.to_csv('Documents/S_boulardii_modeling/results/smetana/communities.tsv', sep = '\t', header = False, index = False)
smetana_table.to_csv('Documents/S_boulardii_modeling/results/20230206/smetana/communities.tsv', sep = '\t', header = False, index = False)
