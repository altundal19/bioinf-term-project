import io 
import pandas as pd 
import gzip 
from datetime import datetime


ref_genome = "data/Homo_sapiens_assembly38.fasta"

default_vcf_header = f"""##fileformat=VCFv4.1
##fileDate={datetime.now().strftime("%Y%m%d")}
##source=myImputationProgramV3.1
##reference={ref_genome}
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"""

def read_vcf(path):
    if path[-3:] == ".gz": 
        with io.TextIOWrapper(gzip.open(path,'r')) as f:
            lines = [l for l in f if not l.startswith('##')]
    else:
        with open(path,"r") as f:
            lines = [l for l in f if not l.startswith('##')]

    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

def df_to_vcf(df, outdir, header=default_vcf_header):
    with open(outdir, 'w') as vcf:
        vcf.write(header)

    df.to_csv(outdir, sep="\t", mode='a', index=False, header=False)

def read_vcf_header(path):
    if path[-3:] == ".gz": 
        with io.TextIOWrapper(gzip.open(path,'r')) as f:
            lines = [l for l in f if l.startswith('#')]
    else:
        with open(path,"r") as f:
            lines = [l for l in f if l.startswith('#')]

    return ''.join(lines)

def filter(path):
    
    vcf_data = read_vcf(path)
    header = read_vcf_header(path)
    filtered_vcf = vcf_data[(vcf_data['FILTER'] == 'PASS') | (vcf_data['FILTER'] == '.')]
    df_to_vcf(filtered_vcf, "out.vcf", header)

def calc_scores(path_pred, path_data):
    pred = read_vcf(path_pred)
    data = read_vcf(path_data)

    pred["key"] = pred['CHROM'] +"-"+ pred['POS'].astype(str) +"-"+ pred['REF']  +"-"+ pred["ALT"]
    data["key"] = data['CHROM']  +"-"+  data['POS'].astype(str) +"-"+  data['REF']  +"-"+ data["ALT"]

    predkey = set(pred['key'])
    datakey = set(data['key'])

    # intersection tp
    tp = predkey.intersection(datakey)

    # # false pos pred-data
    fp = predkey.difference(datakey)

    # # false neg data-pred
    fn = datakey.difference(predkey)

    n_tp = len(tp)
    n_fp = len(fp)
    n_fn = len(fn)
    n = len(predkey)

    # accuracy tp/all
    acc = n_tp/n
    recall = n_tp/(n_tp+n_fn) 
    precision = n_tp/(n_tp+n_fp)
    f1_score = 2*precision*recall/(precision+recall)

    print("Accuracy: ", acc)
    print("Recall: ", recall)
    print("Precision: ", precision)
    print("F1 score: ", f1_score)

calc_scores("data/hc_bed_filtered.recode.vcf", "out.vcf")
