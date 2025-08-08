import pandas as pd
import os
import numpy as np
import csv
import scipy.io

def SL_SS3x_summary_table_hf(SS3x_sheet_path,matrix_dir,neg_matrix_path):
    # load the sheet with the information per well
    bc_wells_df = pd.read_csv(SS3x_sheet_path,sep="\t")

    # make a dictionary associating barcodes and wells
    dic = dict(zip(bc_wells_df.BC_TAG, bc_wells_df.Well_ID))
    
    # Reindex bc_wells_df so that the "Well_ID" column is the index
    bc_wells_df.set_index('Well_ID',inplace=True)
    
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
    
    # rename based on the dictionary so instead of having the barcodes as columns, we have the more informative "Well_ID"
    
    matrix.rename(columns=dic, inplace=True)
    
    # sort columns based on the order in the barcode file list
    matrix = matrix.reindex(bc_wells_df.index, axis=1)
    
    # Remove the "Well_ID" as name of the columns
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
    
    
    ## Create a Summary table
    matrix_summary = pd.DataFrame(columns=matrix.columns)
    
    matrix_summary.loc['UMI_Total'] =  matrix.sum()
    
    matrix_summary.loc['Genes0_Total'] = matrix[matrix > 0].count()
    matrix_summary.loc['Genes1_Total'] = matrix[matrix > 1].count()
    matrix_summary.loc['Genes10_Total'] = matrix[matrix > 10].count()
    matrix_summary.loc['ERCC_Total'] = matrix_ERCC.sum()
    
    matrix_summary_tp = matrix_summary.transpose()
    
    Gene_counts = pd.concat([matrix_summary_tp,bc_wells_df.loc[:,["Plate","Type","Sorting","Library_ID"]]],axis=1)
 
    return Gene_counts

def SL_SS3x_summary_table(SS3x_sheet_path,matrix_dir):
    # load the sheet with the information per well
    bc_wells_df = pd.read_csv(SS3x_sheet_path,sep="\t")

    # make a dictionary associating barcodes and wells
    dic = dict(zip(bc_wells_df.BC_TAG, bc_wells_df.Well_ID))
    
    # Reindex bc_wells_df so that the "Well_ID" column is the index
    bc_wells_df.set_index('Well_ID',inplace=True)
    
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
    
    # rename based on the dictionary so instead of having the barcodes as columns, we have the more informative "Well_ID"
    
    matrix.rename(columns=dic, inplace=True)
    
    # sort columns based on the order in the barcode file list
    matrix = matrix.reindex(bc_wells_df.index, axis=1)
    
    # Remove the "Well_ID" as name of the columns
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

    ## Create a Summary table
    matrix_summary = pd.DataFrame(columns=matrix.columns)
    
    matrix_summary.loc['UMI_Total'] =  matrix.sum()
    
    matrix_summary.loc['Genes0_Total'] = matrix[matrix > 0].count()
    matrix_summary.loc['Genes1_Total'] = matrix[matrix > 1].count()
    matrix_summary.loc['Genes10_Total'] = matrix[matrix > 10].count()
    matrix_summary.loc['ERCC_Total'] = matrix_ERCC.sum()
    
    matrix_summary_tp = matrix_summary.transpose()
    
    Gene_counts = pd.concat([matrix_summary_tp,bc_wells_df.loc[:,["Plate","Type","Sorting","Library_ID"]]],axis=1)
 
    return Gene_counts

def SL_SS3x_matrix(SS3x_sheet_path,matrix_dir):
    # load the sheet with the information per well
    bc_wells_df = pd.read_csv(SS3x_sheet_path,sep="\t")

    # make a dictionary associating barcodes and wells
    dic = dict(zip(bc_wells_df.BC_TAG, bc_wells_df.Well_ID))
    
    # Reindex bc_wells_df so that the "Well_ID" column is the index
    bc_wells_df.set_index('Well_ID',inplace=True)
    
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
    
    # rename based on the dictionary so instead of having the barcodes as columns, we have the more informative "Well_ID"
    
    matrix.rename(columns=dic, inplace=True)
    
    # sort columns based on the order in the barcode file list
    matrix = matrix.reindex(bc_wells_df.index, axis=1)
    
    # Remove the "Well_ID" as name of the columns
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

    return matrix, matrix_ERCC


def SL_SS3x_matrix_hf(SS3x_sheet_path,matrix_dir,neg_matrix_path):
    # load the sheet with the information per well
    bc_wells_df = pd.read_csv(SS3x_sheet_path,sep="\t")

    # make a dictionary associating barcodes and wells
    dic = dict(zip(bc_wells_df.BC_TAG, bc_wells_df.Well_ID))
    
    # Reindex bc_wells_df so that the "Well_ID" column is the index
    bc_wells_df.set_index('Well_ID',inplace=True)
    
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
    
    # rename based on the dictionary so instead of having the barcodes as columns, we have the more informative "Well_ID"
    
    matrix.rename(columns=dic, inplace=True)
    
    # sort columns based on the order in the barcode file list
    matrix = matrix.reindex(bc_wells_df.index, axis=1)
    
    # Remove the "Well_ID" as name of the columns
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
