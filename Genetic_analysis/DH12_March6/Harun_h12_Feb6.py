from sys import argv
import numpy as np
import pandas as pd
import malariagen_data
import allel
from sys import stdout
from pathlib import Path
import os

if (len(argv) == 4):
    cohort = argv[1]
    contig = argv[2]
    num_randomisations = int(argv[3])
else:
    raise Exception("Fail. There should be three command line arguments (cohort, contig, num_randomisations).")


window_size = 1000

window_size = 1000

print("Step 2: Loading metadata...")
# load malariagen af1 data
## loading malariagen data
af1 = malariagen_data.Af1(pre=True,
                          gcs_cache='/home/namulil/lstm_projects/funestus_llineup/gcs_cache',
                          results_cache='home/namulil/lstm_projects/funestus_llineup/results_cache'
                          )

## load metadats
metadata = af1.sample_metadata(
    sample_sets=["1288-VO-UG-DONNELLY-VMF00219"], 
    sample_query = None
)
print("Step 3: Metadata loaded successfully.")

print("Now loading radomisation csv.")
randomisations = pd.read_csv("/home/namulil/lstm_projects/funestus_llineup/notebooks/DH12_permutations/net_round_randomisations_pop.csv", sep=None, engine="python")
print("Randomisations DataFrame (first 5 rows):")
print(randomisations.head())
print (randomisations.shape)
print(f"Is randomisations dataframe empty: {randomisations.empty}")

print ("Now filtering metadata on sample_id column in radomisation data")
#### make sure all metadata samples are in the randomisations file 
metadata = metadata.query(f"sample_id in {randomisations['sample_id'].to_list()}")# filter af1 metadata according to samples-ids in the created randomisation df
ids = metadata['sample_id'] # create the sample id list 
print("ids DataFrame (first 5 rows):")
print(ids.head())
print (ids.shape)
print (metadata.shape)
print(f"Is ids dataframe empty: {ids.empty}")


# make sure randomisations have only samples in the metadata
random = randomisations.query("sample_id in @ids")## Looks redandant to me but keeps rows only where sample id in ids dataframe also apperar in the randomisation column
cohort_ids = random.query(f"population == @cohort")['sample_id'].to_list()
print(f"Getting biallelic mask for {cohort} {contig}")


## calculate/ get haplotypes
ds_haps = af1.haplotypes(
    region=contig,
    analysis="funestus",
    sample_query=f"sample_id in {cohort_ids}",
    sample_sets=None,
)

## print the ds_hap before filter
print (f"ds_hap before filter:{ds_haps.dims}")# arrays do not have shape so using dimesions code
gt = allel.GenotypeArray(ds_haps["call_genotype"].data)## convert genotype data into a numpy array
print (gt.shape)

ac = gt.count_alleles()# count occurance of each allel across the dataset. ie total 1 and 0's
biallelic_mask = ac.is_biallelic() ## filter out only biallelic sites
print (biallelic_mask)
## in the first step we have converted ds genotype calls into a numpy array and filtred the biallic calls to a numpy array called bi-allelic


ds_haps = ds_haps.isel(variants=biallelic_mask)# Applying the fliter of bi allelic calls (according to array above using additional metadata) to our ds_hap that we calculated
ht = gt[biallelic_mask, :].to_haplotypes()
print(f"Checking shape of filtred ht:{ht.shape}") ## printing new filter checking new dataset ht
print(f"Any missing values in ht? {np.isnan(ht).any()}")# check to see any nNAN
missing_count = np.isnan(ht).sum()
print(f"Total missing values in ht: {missing_count}")# counting total missing values

metadata.set_index('sample_id', inplace = True) ## setting sample id as index
unique_sample_id = np.unique(ds_haps['sample_id'])

#dup hap ids
hap_ids = np.repeat(unique_sample_id,2)## since each sample has two haplotypes we are duplicating them out

pos = ds_haps["variant_position"].values # defing 'pos' argument to use later
midpoints = allel.moving_statistic(pos, statistic=np.mean, size=window_size)
startpoints = allel.moving_statistic(pos, statistic=np.min, size=window_size)
endpoints = allel.moving_statistic(pos, statistic=np.max, size=window_size)

#h12
print ("Now creating the pre_post dictionary")
pre_post_dict = {}
pre_post_dict['pre'] = random.query(f"Control_phase == 'Pre' & population == @cohort")['sample_id']
pre_post_dict['post']  = random.query(f"Control_phase == 'Post' & population == @cohort")['sample_id']
h12_df_dict = {}
for pheno in ['pre', 'post']: ## keys in the dictionary

    subsample_ht = ht[:, np.isin(hap_ids, pre_post_dict[pheno])]
    print(f"--------- Running H12 on {cohort} {pheno} | Chromosome {contig} ----------")
    ## checking ht_subset before calling the moving garuds
    print(f"subsample_ht shape: {subsample_ht.shape}")
    print(f"subsample_ht data:\n{subsample_ht}")

    h1, h12, h123, h2_h1 = allel.moving_garud_h(subsample_ht, size=window_size)
    h12_df_dict[pheno] = pd.DataFrame({'startpoint': startpoints, 
                                       'endpoint': endpoints, 
                                       'midpoint':midpoints, 
                                       'h12':h12})


#h12 randomisation

## start by checking columns available in the radom data
print(f"Available columns in random DataFrame: {random.columns}")

randomised_h12_dict = {'pre': dict(), 'post': dict()}
for i in np.arange(1, num_randomisations+1):
    i = '%0.4d' % i


    pre_post_dict['pre'] = random.query(f"`r{i}` == 'Pre' & population == @cohort")['sample_id']
    pre_post_dict['post']  = random.query(f"`r{i}` == 'Post' & population == @cohort")['sample_id']

    for pheno in ['pre', 'post']:

        subsample_ht = ht[:, np.isin(hap_ids, pre_post_dict[pheno])]
        print(f"--------- Running H12 on {cohort} {pheno} randomisation {i} | Chromosome {contig} ----------")
        h1, h12, h123, h2_h1 = allel.moving_garud_h(subsample_ht, size=window_size)
        randomised_h12_dict[pheno][f"r{i}"] = h1


#save
randomised_h12_df_pre = pd.DataFrame(randomised_h12_dict['pre']) ## save the pre results to df
randomised_h12_df_post = pd.DataFrame(randomised_h12_dict['post'])# save post h12 to df
h12_df_pre = pd.concat([h12_df_dict['pre'], randomised_h12_df_pre], axis = 1)
h12_df_post = pd.concat([h12_df_dict['post'], randomised_h12_df_post], axis = 1)

## checking if the pre df was created
print(f"Keys in h12_df_dict: {list(h12_df_dict.keys())}")


#h12_df_pre.to_csv("/home/namulil/lstm_projects/funestus_llineup/notebooks/DH12_permutations_{}.pre.{}.tsv".format(cohort, contig), sep="\t", index=False)
#h12_df_post.to_csv("/home/namulil/lstm_projects/funestus_llineup/notebooks/DH12_permutations_{}.post.{}.tsv".format(cohort, contig), sep="\t", index=False)

## improving code to see if we ca save straight in the DH12 folder, commented out code saves in notebooks
output_dir = "/home/namulil/lstm_projects/funestus_llineup/notebooks/DH12_permutations" # Define the directory

os.makedirs(output_dir, exist_ok=True) # Ensure the directory exists


# Save the file
file_path_pre = os.path.join(output_dir, "DH12_permutations_{}.pre.{}.tsv".format(cohort, contig))
file_path_post = os.path.join(output_dir, "DH12_permutations_{}.post.{}.tsv".format(cohort, contig))

h12_df_pre.to_csv(file_path_pre, sep="\t", index=False)
h12_df_post.to_csv(file_path_post, sep="\t", index=False)


