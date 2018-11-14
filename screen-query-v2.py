'''
Unable to get mBio calls. This data is in RESAMPLING datasets.
'''


import pandas as pd
import os
import sys

if len(sys.argv) == 0:
    genes_of_interest = "RVBD0667".split(",")
else:
    temp = sys.argv[1]
    # If file exists, then it's actually a path:
    if os.path.exists(temp):
        genes_of_interest = []
        for line in open(temp):
            # Get genes in tab-separated file in first column
            # (works with just one column and no tabs)
            gene = line.strip().split("\t")[0]
            genes_of_interest.append(gene)
    # If it doesn't exist, must be a comma-separated list of names:
    else:
        genes_of_interest = sys.argv[1].split(",")

csv_folder = "passaging and drug screen results"

def get_assays(csv_folder):
    assays = os.listdir(csv_folder)
    assays.sort()
    assay_list = []

    '''
    # Get the base of the assay
    for assay in assays:

            split_assay = assay.split(" mageck.csv")
            assay_base = split_assay[0]
            assay_list.append(assay_base)
    '''
    return assays

def trim_assay_names(assay_list):
    trimmed_assay_list = []

    for assay in assay_list:
        trimmed_assay_list.append(assay.split(".mageck")[0].split("result_")[1])

    return trimmed_assay_list


# NOTE: if you are getting an error message, it is probably because there are doubles of your gene
# Try adding a colon to the end of gene_of_interest
def gene_stats(gene_of_interest, csv_folder, assay_list, list_gene_df):

    gene_found = True

    trimmed_assay_list = trim_assay_names(assay_list)
    # Initializes the dataframe and all of the columns of interest
    gene_dataframe = pd.DataFrame(data=trimmed_assay_list,columns=["assay"])
    gene_dataframe['gene'] = gene_of_interest
#    gene_dataframe['mBio Call'] = "-"
    gene_dataframe['num sgRNAs'] = "-"
    gene_dataframe['mageck neg rank'] = "-"
    gene_dataframe['mageck neg p value'] = "-"
    gene_dataframe['mageck neg FDR'] = "-"
    gene_dataframe['mageck pos rank'] = "-"
    gene_dataframe['mageck pos p value'] = "-"
    gene_dataframe['mageck pos FDR'] = "-"
    gene_dataframe['mageck log fold change'] = "-"
    gene_dataframe['mageck fold change'] = "-"

    # Loops through each csv file in the folder
    for assay in assay_list:

        # Reads data into pandas
        mageck_assay_data = pd.read_csv(csv_folder + "/" + assay, delimiter="\t")

        # Gives you the row in the new dataframe where the current assay is
        assay_booleans = gene_dataframe['assay'].str.contains(assay.split(".mageck")[0].split("result_")[1])
        assay_info = gene_dataframe[assay_booleans]
        assay_index = assay_info.index

        # Gives you the row in the screen data where your gene of interest is
        gene_booleans = mageck_assay_data['id'].str.contains(gene_of_interest)
        gene_info = mageck_assay_data[gene_booleans]
        gene_index = gene_info.index

        try:
            lfc = float(mageck_assay_data['neg|lfc'].loc[gene_index])

        except:
            gene_found = False
            print (gene_of_interest + " not found.")
            break

        gene_dataframe.at[assay_index, 'num sgRNAs'] = int(mageck_assay_data['num'].loc[gene_index])

        gene_dataframe.at[assay_index,'mageck log fold change'] = lfc

        if lfc < 0 :
            gene_dataframe.at[assay_index, 'mageck fold change'] = -1/(2 ** lfc)

        else:
            gene_dataframe.at[assay_index,'mageck fold change'] = 2 ** lfc

        neg_p_value = float(mageck_assay_data['neg|p-value'].loc[gene_index])
        gene_dataframe.at[assay_index,'mageck neg p value'] = neg_p_value

        pos_p_value = float(mageck_assay_data['pos|p-value'].loc[gene_index])
        gene_dataframe.at[assay_index,'mageck pos p value'] = pos_p_value

        neg_rank = float(mageck_assay_data['neg|rank'].loc[gene_index])
        gene_dataframe.at[assay_index,'mageck neg rank'] = neg_rank

        pos_rank = float(mageck_assay_data['pos|rank'].loc[gene_index])
        gene_dataframe.at[assay_index,'mageck pos rank'] = pos_rank

        gene_dataframe.at[assay_index, 'mageck neg FDR'] = float(mageck_assay_data['neg|fdr'].loc[gene_index])
        gene_dataframe.at[assay_index, 'mageck pos FDR'] = float(mageck_assay_data['pos|fdr'].loc[gene_index])

        '''
        gene_resampling_booleans = resampling_assay_data['gene'].str.contains(gene_of_interest)
        gene_resampling_info = resampling_assay_data[gene_resampling_booleans]
        gene_resampling_index = gene_resampling_info.index
        '''
        gene_name = str(mageck_assay_data['id'].loc[gene_index])
        gene_name = gene_name.split("    ")
        gene_name = gene_name[1].split("\n")

        gene_dataframe['gene'].replace(to_replace=gene_of_interest, value=gene_name[0], inplace=True)
        '''
        call = str(resampling_assay_data['mBio Call'].loc[gene_resampling_index])
        call = call.split("    ")
        call = call[1].split("\n")

        gene_dataframe['mBio Call'].replace(to_replace="-", value=call[0], inplace=True)
        '''

    if not gene_found:
        return

    list_gene_df.append(gene_dataframe)
 #   gene_dataframe.to_csv(path_or_buf = gene_of_interest + " passaging and drug screen results.csv")

assay_list = get_assays(csv_folder)

list_gene_df = []

for gene_of_interest in genes_of_interest:
    print gene_of_interest
    gene_stats(gene_of_interest, csv_folder, assay_list, list_gene_df)

all_genes_df = pd.concat(list_gene_df)
all_genes_df.to_csv(path_or_buf = sys.argv[1] + " passaging and drug screen results.csv")