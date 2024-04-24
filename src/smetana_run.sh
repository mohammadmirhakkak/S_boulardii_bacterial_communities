#!/bin/bash


# overall results
smetana Documents/S_boulardii_modeling/models/20230206/smetana/*.xml -c \
Documents/S_boulardii_modeling/results/20230206/smetana/communities.tsv \
-m SC --mediadb Documents/S_boulardii_modeling/data/20230206/smetana_media.tsv \
--output Documents/S_boulardii_modeling/results/20230206/smetana/ --flavor ucsd --detailed


smetana Documents/S_boulardii_modeling/models/20230206/smetana/*.xml -c \
Documents/S_boulardii_modeling/results/20230206/smetana/communities_test.tsv \
-m SC --mediadb Documents/S_boulardii_modeling/data/20230206/smetana_media.tsv \
--output Documents/S_boulardii_modeling/results/20230206/smetana/ --flavor ucsd --detailed


smetana Documents/S_boulardii_modeling/models/20230206/smetana/*.xml -c \
Documents/S_boulardii_modeling/results/20230206/smetana/communities_test.tsv \
--output Documents/S_boulardii_modeling/results/20230206/smetana/ --flavor ucsd --detailed