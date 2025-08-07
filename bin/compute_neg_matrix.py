import argparse
import subprocess
import numpy as np
import pandas as pd
import time


def bam2sam(bamfile, outname, samtools_path=None):
    if samtools_path is None:
        subprocess.call(f"samtools view -h -o {outname}.sam {bamfile}", shell=True)
    else:
        subprocess.call(f"{samtools_path} view -h -o {outname}.sam {bamfile}", shell=True)


def sam2tsv(input_file, output_file):
    cmd = f'''
    #!/usr/bin/bash
    HEADER="gene_id\\tindex1\\tindex2\\ttag\\tumi"
    INPUT_FILE={input_file}
    OUTPUT_FILE={output_file}

    echo -e $HEADER > $OUTPUT_FILE
    # extract rows with reads (skip metadata) from sam file 
    # extract columns GN, CB, UB 
    # extract rows without "-" in columns 
    # remove STAR attribute tags| parse CB to index1, index2, TAG
    grep -E 'GN:Z' $INPUT_FILE|cut -f17,19,20|grep -v ':-'|sed -r 's/(..:Z:)//g'|awk '{{s1=substr($2,1,8)}}{{s2=substr($2,9,8)}}{{s3=substr($2,17,11)}}{{print $1"\\t"s1"\\t"s2"\\t"s3"\\t"$3}}'>>$OUTPUT_FILE
    rm {input_file}
    '''
    subprocess.call(cmd, shell=True)


def generate_negative_matrix(filename, threshold=0.8):
    df = pd.read_csv(filename, sep='\t')
    barcodes = sorted(set(df.index1 + df.index2 + df.tag))
    genes = sorted(set(df.gene_id))
    tags = sorted(set(df.tag))
    
    res = df.groupby(['tag', 'umi', 'gene_id', 'index1', 'index2']).size().to_frame()
    res.rename(columns={0: 'counts'}, inplace=True)
    res['tot1'] = res.assign(tot1=lambda x: x.groupby(['tag', 'umi', 'gene_id', 'index1'])['counts'].sum())['tot1']
    res = res.swaplevel('index1', 'index2').assign(tot2=lambda x: x.groupby(['tag', 'umi', 'gene_id', 'index2'])['counts'].sum())
    res['filter'] = res['counts'] < threshold * (res['tot1'] + res['tot2'] - res['counts'])
    res = res[res['filter']].reset_index().drop(columns=['counts', 'tot1', 'tot2', 'filter'])
    res['CB'] = res['index1'] + res['index2'] + res['tag']
    res.drop(columns=['index1', 'index2', 'tag'], inplace=True)
    res = res.groupby(['CB', 'gene_id']).nunique()
    res = res.unstack(fill_value=0)
    res.columns = res.columns.droplevel()
    res.index.name = None
    res.columns.name = None
    count_mtx = res.reindex(index=barcodes, columns=genes, copy=True)
    count_mtx.fillna(0, inplace=True)
    subprocess.call(f"rm {filename}", shell=True)
    return count_mtx


if __name__ == '__main__':
    start_time = time.time()
    parser = argparse.ArgumentParser(
                    prog='ComputeNegativeMatrix',
                    description='Generate the negative count matrix for hopped reads in sequencing experiments',
                    epilog='It is recommended to run with screen or srun')
    parser.add_argument('filename', help='bam file with raw sequencing counts')
    parser.add_argument('outname', help='output filename')
    parser.add_argument('--samtools_path', help='path to samtools tool', type=str)
    parser.add_argument('--threshold', help='filtering threshold', type=float)
    args = parser.parse_args()
    bam2sam(args.filename, args.outname, args.samtools_path)
    sam2tsv(f"{args.outname}.sam", f"{args.outname}.tsv")
    threshold = 0.8 if args.threshold is None else args.threshold
    neg_matrix = generate_negative_matrix(f"{args.outname}.tsv", threshold=threshold)
    neg_matrix.to_csv(f"{args.outname}_negative_matrix.csv")
    print(f'Elapsed time: {time.time() - start_time} seconds.')