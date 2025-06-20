# %% [markdown]
# # Extracting whole genome genotype data 
# 
# This is a script to extract genotypes from whole genomes of *An.funestus* collected during the LLINEUP-PBO trial that was conducted in Uganda from 2017-2019.  
# 
# 

# %%
import os
print(os.getcwd())

# %%
%pip install malariagen_data==9.0.1 numpy -qq

# %%
#Install and load packages
import malariagen_data
import os
import numpy as np
import pandas as pd
import allel
import xarray as xr
import glob

# %%
#define directory to save malariagen analysis results and load the prelease data
results_dir = "~/llineup-genomics-main"
os.makedirs(results_dir, exist_ok=True)
af1 = malariagen_data.Af1(pre=True,
                          gcs_cache='/home/namulil/lstm_projects/funestus_llineup/gcs_cache',
                          results_cache='home/namulil/lstm_projects/funestus_llineup/results_cache'
                          )

# %% [markdown]
# ##  Access Uganda PBO trial *An. funestus* genotypes

# %%

#function to extract biallelic genotypes

def extract_and_filter_snps(region, maf_threshold, output_path):
    
    array_snps = af1.snp_calls(region=region,
                               sample_sets=["1288-VO-UG-DONNELLY-VMF00219"],
                               sample_query= None,
                               site_mask='funestus' )
    
    gt = allel.GenotypeArray(array_snps['call_genotype'])
    
    no_missing = gt.count_missing(1) == 0
    gt_freq=gt.count_alleles().to_frequencies()
    which_pos = (np.max(gt_freq,1) < (1 - maf_threshold)) & no_missing
    gt_filtered = gt[which_pos,:]
   
    gt_biallelic = np.sum(gt_filtered>0,2)
    #convert to dataframe
    df_gt = pd.DataFrame(gt_biallelic)
    pos = array_snps['variant_position'][which_pos]
    
    chrom = np.array(array_snps.contigs)[array_snps['variant_contig']][which_pos]
    #snp_id = np.apply_along_axis(':'.join, 0, [chrom, pos.astype('str')])
    snp_id = np.apply_along_axis(lambda x: np.asarray(':'.join(x), dtype = 'object'), 0, [chrom, pos.values.astype('str')])
    df_gt.set_index(snp_id, inplace = True)
    df_gt.columns=array_snps.sample_id
    return(df_gt)
    

# %%
output_directory = '//home/namulil/lstm_projects/funestus_llineup/R analysis'

regions = ['2RL', '3RL', 'X']

# Dictionary comprehension to call the function for each region
gt = {region: extract_and_filter_snps(region, 0.02, os.path.join(output_directory, f'gt_{region}.csv')) for region in regions}


# %%
print (gt)



# %%
if not df_list:
    print("No CSV files found.")
else:
    df_gt = pd.concat(df_list, axis=0, ignore_index=True)
    df_gt.to_csv("//home/namulil/lstm_projects/funestus_llineup/R analysis")
    


