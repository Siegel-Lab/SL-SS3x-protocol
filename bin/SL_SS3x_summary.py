import pandas as pd
import os
import numpy as np
import csv
import scipy.io

def SS3x_summary_table_hf(SS3x_sheet_path,matrix_dir,neg_matrix_path,VSG_loc_path):
    # load the sheet with the information per well
    bc_wells_df = pd.read_csv(SS3x_sheet_path,sep="\t")
    
    # make a dictionary associating barcodes and wells 
    dic = dict(zip(bc_wells_df.BC_TAG, bc_wells_df.Identity_code))
    
    # Reindex bc_wells_df so that the "Identity_code" column is the index
    bc_wells_df.set_index('Identity_code',inplace=True)
    
    # creating a matrix with summary data
    mat = scipy.io.mmread(os.path.join(matrix_dir, "matrix.mtx"))
    # list of transcript ids, e.g. 'ENSG00000243485'
    features_path = os.path.join(matrix_dir, "features.tsv")
    feature_ids = [row[0] for row in csv.reader(open(features_path, mode="rt"), delimiter="\t")]

    # list of barcodes
    barcodes_path = os.path.join(matrix_dir, "barcodes.tsv")
    barcodes = [row[0] for row in csv.reader(open(barcodes_path, mode="rt"), delimiter="\t")]

    # transform table to pandas dataframe and label rows and columns
    matrix = pd.DataFrame.sparse.from_spmatrix(mat)
    matrix = matrix.sparse.to_dense()
    matrix.columns = barcodes

    # rename based on the dictionary so instead of having the barcodes as columns, we have the more informative "Identity code"

    matrix.rename(columns=dic, inplace=True)

    # sort columns based on the order in the barcode file list
    matrix = matrix.reindex(bc_wells_df.index, axis=1)

    # Remove the "Identity_code" as name of the columns
    matrix.columns.name = None

    # insert the column with the feature_ids

    matrix.insert(loc=0, column="feature_id", value=feature_ids)

    #subset matrix to only those ERCCs 

    matrix_ERCC = matrix.loc[matrix['feature_id'].str.contains('gSL_ERCC')]

    #remove ERCC from the original matrix

    matrix = matrix.loc[~matrix['feature_id'].str.contains('gSL_ERCC')]

    #set 'feature_id' as row index and remove the index name from the matrices

    matrix = matrix.set_index('feature_id')
    matrix_ERCC = matrix_ERCC.set_index('feature_id')
    matrix.index.name = None
    matrix_ERCC.index.name = None

    # Read in the VSG_location file into a dataframe

    VSG_location_df = pd.read_csv(VSG_loc_path, sep="\t")

    ##################### INDEX HOPPING REMOVAL SECTION ##############################

    neg_matrix = pd.read_csv(neg_matrix_path)

    ## format it in the same format that the standard matrix has

    t_neg_matrix = neg_matrix.transpose()


    #set column names
    t_neg_matrix.columns = t_neg_matrix.iloc[0]

    #drop the first row
    t_neg_matrix = t_neg_matrix.iloc[1:]

    # Remove the "Unnnamed: 0" as name of the columns
    t_neg_matrix.columns.name = None


    t_neg_matrix.rename(columns=dic, inplace=True)

    # sort columns based on the order in the barcode file list
    t_neg_matrix = t_neg_matrix.reindex(bc_wells_df.index, axis=1)

    # insert the column with the feature_ids (to do the filtering)
    t_neg_matrix['feature_id'] = t_neg_matrix.index


    #subset matrix to only those ERCCs

    t_neg_matrix_ERCC = t_neg_matrix.loc[t_neg_matrix['feature_id'].str.contains('gSL_ERCC')]

    #remove ERCC from the original matrix

    t_neg_matrix = t_neg_matrix.loc[~t_neg_matrix['feature_id'].str.contains('gSL_ERCC')]

    #set 'feature_id' as row index and remove the index name from the matrices

    t_neg_matrix = t_neg_matrix.set_index('feature_id')
    t_neg_matrix_ERCC = t_neg_matrix_ERCC.set_index('feature_id')

    t_neg_matrix.index.name = None
    t_neg_matrix_ERCC.index.name = None

    # substract the negative matrix from the original

    matrix_hf = matrix.sub(t_neg_matrix, fill_value=0)

    matrix_ERCC_hf = matrix_ERCC.sub(t_neg_matrix_ERCC, fill_value=0)


    # Reassign the values to the original matrices

    matrix = matrix_hf

    matrix_ERCC = matrix_ERCC_hf


    ##################### INDEX HOPPING REMOVAL SECTION ##############################

    
    ## Create a Summary table
    matrix_summary = pd.DataFrame(columns=matrix.columns)

    matrix_summary.loc['UMI_Total'] =  matrix.sum()

    matrix_summary.loc['Genes0_Total'] = matrix[matrix > 0].count()
    matrix_summary.loc['Genes1_Total'] = matrix[matrix > 1].count()
    matrix_summary.loc['Genes10_Total'] = matrix[matrix > 10].count()
    matrix_summary.loc['ERCC_Total'] = matrix_ERCC.sum()

	# Remove the ":pseudogene" from the Gene ID in the index of the matrix from the ones that have it, so i can find exact matches based on the Gene_ID in VSG_location_sorted_df , which IDs do not have the ":pseudogene" string

    matrix.index = matrix.index.str.replace(':pseudogene','')

    # Make a column with the total VSG UMI counts
    matrix_summary.loc['VSG_set_UMI_Total'] = matrix.loc[matrix.index.isin(VSG_location_df["Gene_ID"])].sum()
	#   matrix_summary.loc['VSG_set_UMI_Total'] = matrix.loc[matrix.index.str.contains('|'.join(VSG_location_df["Gene_ID"]))].sum()

    # Make a column with the ratio of VSG UMI counts over ERCC UMI counts
    matrix_summary.loc['ratio_VSG_over_ERCC_UMI_Total'] = matrix_summary.loc['VSG_set_UMI_Total'].astype(float).div(matrix_summary.loc['ERCC_Total'].astype(float)).replace(np.inf, np.nan)

    #  Make a column with the ratio of VSG UMI counts over genes UMI counts

    matrix_summary.loc['ratio_VSG_over_genes_UMI_Total'] = matrix_summary.loc['VSG_set_UMI_Total'].astype(float).div(matrix_summary.loc['UMI_Total'].astype(float)).replace(np.inf, np.nan)

    
    # Make a column with the total Non-VSG UMI counts (Total UMI counts - VSG UMI counts)

    matrix_summary.loc['non_VSG_UMI_Total'] = matrix_summary.loc['UMI_Total'] - matrix_summary.loc['VSG_set_UMI_Total']


    # Make a column with the ratio of Non-VSG UMI counts (Total UMI counts - VSG UMI counts) over ERCC UMI counts
    matrix_summary.loc['ratio_non_VSG_over_ERCC_UMI_Total'] = (matrix_summary.loc['UMI_Total'] - matrix_summary.loc['VSG_set_UMI_Total']).astype(float).div(matrix_summary.loc['ERCC_Total'].astype(float)).replace(np.inf, np.nan)


    matrix_summary_tp = matrix_summary.transpose()

    Gene_counts = pd.concat([matrix_summary_tp,bc_wells_df.iloc[:,[13,14,15,16,45]]],axis=1)

    return Gene_counts

def SS3x_summary_table(SS3x_sheet_path,matrix_dir,VSG_loc_path):
    # load the sheet with the information per well
    bc_wells_df = pd.read_csv(SS3x_sheet_path,sep="\t")
    
    # make a dictionary associating barcodes and wells 
    dic = dict(zip(bc_wells_df.BC_TAG, bc_wells_df.Identity_code))
    
    # Reindex bc_wells_df so that the "Identity_code" column is the index
    bc_wells_df.set_index('Identity_code',inplace=True)
    
    # creating a matrix with summary data
    mat = scipy.io.mmread(os.path.join(matrix_dir, "matrix.mtx"))
    # list of transcript ids, e.g. 'ENSG00000243485'
    features_path = os.path.join(matrix_dir, "features.tsv")
    feature_ids = [row[0] for row in csv.reader(open(features_path, mode="rt"), delimiter="\t")]

    # list of barcodes
    barcodes_path = os.path.join(matrix_dir, "barcodes.tsv")
    barcodes = [row[0] for row in csv.reader(open(barcodes_path, mode="rt"), delimiter="\t")]

    # transform table to pandas dataframe and label rows and columns
    matrix = pd.DataFrame.sparse.from_spmatrix(mat)
    matrix = matrix.sparse.to_dense()
    matrix.columns = barcodes

    # rename based on the dictionary so instead of having the barcodes as columns, we have the more informative "Identity code"

    matrix.rename(columns=dic, inplace=True)

    # sort columns based on the order in the barcode file list
    matrix = matrix.reindex(bc_wells_df.index, axis=1)

    # Remove the "Identity_code" as name of the columns
    matrix.columns.name = None

    # insert the column with the feature_ids

    matrix.insert(loc=0, column="feature_id", value=feature_ids)

    #subset matrix to only those ERCCs

    matrix_ERCC = matrix.loc[matrix['feature_id'].str.contains('gSL_ERCC')]

    #remove ERCC from the original matrix

    matrix = matrix.loc[~matrix['feature_id'].str.contains('gSL_ERCC')]

    #set 'feature_id' as row index and remove the index name from the matrices

    matrix = matrix.set_index('feature_id')
    matrix_ERCC = matrix_ERCC.set_index('feature_id')
    matrix.index.name = None
    matrix_ERCC.index.name = None

    # Read in the VSG_location file into a dataframe

    VSG_location_df = pd.read_csv(VSG_loc_path, sep="\t")

    
    ## Create a Summary table
    matrix_summary = pd.DataFrame(columns=matrix.columns)

    matrix_summary.loc['UMI_Total'] =  matrix.sum()

    matrix_summary.loc['Genes0_Total'] = matrix[matrix > 0].count()
    matrix_summary.loc['Genes1_Total'] = matrix[matrix > 1].count()
    matrix_summary.loc['Genes10_Total'] = matrix[matrix > 10].count()
    matrix_summary.loc['ERCC_Total'] = matrix_ERCC.sum()


	# Remove the ":pseudogene" from the Gene ID in the index of the matrix from the ones that have it, so i can find exact matches based on the Gene_ID in VSG_location_sorted_df , which IDs do not have the ":pseudogene" string

    matrix.index = matrix.index.str.replace(':pseudogene','')

	# Make a column with the total VSG UMI counts
    matrix_summary.loc['VSG_set_UMI_Total'] = matrix.loc[matrix.index.isin(VSG_location_df["Gene_ID"])].sum()

    # Make a column with the ratio of VSG UMI counts over ERCC UMI counts
    matrix_summary.loc['ratio_VSG_over_ERCC_UMI_Total'] = matrix_summary.loc['VSG_set_UMI_Total'].astype(float).div(matrix_summary.loc['ERCC_Total'].astype(float)).replace(np.inf, np.nan)

    #  Make a column with the ratio of VSG UMI counts over genes UMI counts

    matrix_summary.loc['ratio_VSG_over_genes_UMI_Total'] = matrix_summary.loc['VSG_set_UMI_Total'].astype(float).div(matrix_summary.loc['UMI_Total'].astype(float)).replace(np.inf, np.nan)
    
    # Make a column with the total Non-VSG UMI counts (Total UMI counts - VSG UMI counts)

    matrix_summary.loc['non_VSG_UMI_Total'] = matrix_summary.loc['UMI_Total'] - matrix_summary.loc['VSG_set_UMI_Total']


    # Make a column with the ratio of Non-VSG UMI counts (Total UMI counts - VSG UMI counts) over ERCC UMI counts
    matrix_summary.loc['ratio_non_VSG_over_ERCC_UMI_Total'] = (matrix_summary.loc['UMI_Total'] - matrix_summary.loc['VSG_set_UMI_Total']).astype(float).div(matrix_summary.loc['ERCC_Total'].astype(float)).replace(np.inf, np.nan)

    matrix_summary_tp = matrix_summary.transpose()

    Gene_counts = pd.concat([matrix_summary_tp,bc_wells_df.iloc[:,[13,14,15,16,45]]],axis=1)

    return Gene_counts

def SS3x_matrix(SS3x_sheet_path,matrix_dir):
    # load the sheet with the information per well
    bc_wells_df = pd.read_csv(SS3x_sheet_path,sep="\t")
    
    # make a dictionary associating barcodes and wells 
    dic = dict(zip(bc_wells_df.BC_TAG, bc_wells_df.Identity_code))
    
    # Reindex bc_wells_df so that the "Identity_code" column is the index
    bc_wells_df.set_index('Identity_code',inplace=True)
    
    # creating a matrix with summary data
    mat = scipy.io.mmread(os.path.join(matrix_dir, "matrix.mtx"))
    # list of transcript ids, e.g. 'ENSG00000243485'
    features_path = os.path.join(matrix_dir, "features.tsv")
    feature_ids = [row[0] for row in csv.reader(open(features_path, mode="rt"), delimiter="\t")]

    # list of barcodes
    barcodes_path = os.path.join(matrix_dir, "barcodes.tsv")
    barcodes = [row[0] for row in csv.reader(open(barcodes_path, mode="rt"), delimiter="\t")]

    # transform table to pandas dataframe and label rows and columns
    matrix = pd.DataFrame.sparse.from_spmatrix(mat)
    matrix = matrix.sparse.to_dense()
    matrix.columns = barcodes

    # rename based on the dictionary so instead of having the barcodes as columns, we have the more informative "Identity code"

    matrix.rename(columns=dic, inplace=True)

    # sort columns based on the order in the barcode file list
    matrix = matrix.reindex(bc_wells_df.index, axis=1)

    # Remove the "Identity_code" as name of the columns
    matrix.columns.name = None

    # insert the column with the feature_ids

    matrix.insert(loc=0, column="feature_id", value=feature_ids)

    #subset matrix to only those ERCCs

    matrix_ERCC = matrix.loc[matrix['feature_id'].str.contains('gSL_ERCC')]

    #remove ERCC and external from the original matrix

    matrix = matrix.loc[~matrix['feature_id'].str.contains('gSL_ERCC')]

    #set 'feature_id' as row index and remove the index name from the matrices

    matrix = matrix.set_index('feature_id')
    matrix_ERCC = matrix_ERCC.set_index('feature_id')
    matrix.index.name = None
    matrix_ERCC.index.name = None

    return matrix, matrix_ERCC


def SS3x_matrix_hf(SS3x_sheet_path,matrix_dir,neg_matrix_path):
    # load the sheet with the information per well
    bc_wells_df = pd.read_csv(SS3x_sheet_path,sep="\t")
    
    # make a dictionary associating barcodes and wells 
    dic = dict(zip(bc_wells_df.BC_TAG, bc_wells_df.Identity_code))
    
    # Reindex bc_wells_df so that the "Identity_code" column is the index
    bc_wells_df.set_index('Identity_code',inplace=True)
    
    # creating a matrix with summary data
    mat = scipy.io.mmread(os.path.join(matrix_dir, "matrix.mtx"))
    # list of transcript ids, e.g. 'ENSG00000243485'
    features_path = os.path.join(matrix_dir, "features.tsv")
    feature_ids = [row[0] for row in csv.reader(open(features_path, mode="rt"), delimiter="\t")]

    # list of barcodes
    barcodes_path = os.path.join(matrix_dir, "barcodes.tsv")
    barcodes = [row[0] for row in csv.reader(open(barcodes_path, mode="rt"), delimiter="\t")]

    # transform table to pandas dataframe and label rows and columns
    matrix = pd.DataFrame.sparse.from_spmatrix(mat)
    matrix = matrix.sparse.to_dense()
    matrix.columns = barcodes

    # rename based on the dictionary so instead of having the barcodes as columns, we have the more informative "Identity code"

    matrix.rename(columns=dic, inplace=True)

    # sort columns based on the order in the barcode file list
    matrix = matrix.reindex(bc_wells_df.index, axis=1)

    # Remove the "Identity_code" as name of the columns
    matrix.columns.name = None

    # insert the column with the feature_ids

    matrix.insert(loc=0, column="feature_id", value=feature_ids)

    #subset matrix to only those ERCCs

    matrix_ERCC = matrix.loc[matrix['feature_id'].str.contains('gSL_ERCC')]


    #remove ERCC from the original matrix

    matrix = matrix.loc[~matrix['feature_id'].str.contains('gSL_ERCC')]

    #set 'feature_id' as row index and remove the index name from the matrices

    matrix = matrix.set_index('feature_id')
    matrix_ERCC = matrix_ERCC.set_index('feature_id')
    matrix.index.name = None
    matrix_ERCC.index.name = None


    ##################### INDEX HOPPING REMOVAL SECTION ##############################

    neg_matrix = pd.read_csv(neg_matrix_path)

    ## format it in the same format that the standard matrix has

    t_neg_matrix = neg_matrix.transpose()


    #set column names
    t_neg_matrix.columns = t_neg_matrix.iloc[0]

    #drop the first row
    t_neg_matrix = t_neg_matrix.iloc[1:]

    # Remove the "Unnnamed: 0" as name of the columns
    t_neg_matrix.columns.name = None


    t_neg_matrix.rename(columns=dic, inplace=True)

    # sort columns based on the order in the barcode file list
    t_neg_matrix = t_neg_matrix.reindex(bc_wells_df.index, axis=1)

    # insert the column with the feature_ids (to do the filtering)
    t_neg_matrix['feature_id'] = t_neg_matrix.index


    #subset matrix to only those ERCCs

    t_neg_matrix_ERCC = t_neg_matrix.loc[t_neg_matrix['feature_id'].str.contains('gSL_ERCC')]

    #remove ERCC from the original matrix

    t_neg_matrix = t_neg_matrix.loc[~t_neg_matrix['feature_id'].str.contains('gSL_ERCC')]

    #set 'feature_id' as row index and remove the index name from the matrices

    t_neg_matrix = t_neg_matrix.set_index('feature_id')
    t_neg_matrix_ERCC = t_neg_matrix_ERCC.set_index('feature_id')

    t_neg_matrix.index.name = None
    t_neg_matrix_ERCC.index.name = None

    # substract the negative matrix from the original

    matrix_hf = matrix.sub(t_neg_matrix, fill_value=0)

    matrix_ERCC_hf = matrix_ERCC.sub(t_neg_matrix_ERCC, fill_value=0)

    # Reassign the values to the original matrices

    matrix = matrix_hf

    matrix_ERCC = matrix_ERCC_hf



    ##################### INDEX HOPPING REMOVAL SECTION ##############################

    
    return matrix, matrix_ERCC
