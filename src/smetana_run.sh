#!/bin/bash


# overall results (default media)
smetana mohammadmirhakkak/S_boulardii_bacterial_communities/GEMs/*.xml -c \
mohammadmirhakkak/S_boulardii_bacterial_communities/res/communities.tsv \
--output mohammadmirhakkak/S_boulardii_bacterial_communities/res/ --flavor ucsd

# detailed results (default media)
smetana mohammadmirhakkak/S_boulardii_bacterial_communities/GEMs/*.xml -c \
mohammadmirhakkak/S_boulardii_bacterial_communities/res/communities.tsv \
--output mohammadmirhakkak/S_boulardii_bacterial_communities/res/ --flavor ucsd --detailed

# overall results (SC media)
smetana mohammadmirhakkak/S_boulardii_bacterial_communities/GEMs/*.xml -c \
mohammadmirhakkak/S_boulardii_bacterial_communities/res/communities.tsv \
-m SC --mediadb mohammadmirhakkak/S_boulardii_bacterial_communities/dat/smetana_media.tsv \
--output mohammadmirhakkak/S_boulardii_bacterial_communities/res/SC --flavor ucsd

# detailed results (SC media)
smetana mohammadmirhakkak/S_boulardii_bacterial_communities/GEMs/*.xml -c \
mohammadmirhakkak/S_boulardii_bacterial_communities/res/communities.tsv \
-m SC --mediadb mohammadmirhakkak/S_boulardii_bacterial_communities/dat/smetana_media.tsv \
--output mohammadmirhakkak/S_boulardii_bacterial_communities/res/SC --flavor ucsd --detailed
