import os
from scipy.stats import norm
import numpy
import math
import argparse
import time
import sys
import pandas as pd
import sqlite3


#parse arguments
parser = argparse.ArgumentParser(description = "CHEERS computes the disease enrichment within functional annotations by taking into account quantitative changes in read counts within peaks", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--input", help = "Text file containing peak coordinates and specificity scores for each of the analyzed samples")
parser.add_argument("--ld", help = "Directory with LD information for each SNP if CHEERS is used on all SNPs in LD")
parser.add_argument("--snp_list", help = "list of SNPs if CHEERS is used on finemapped set")
parser.add_argument("--trait", help = "Name of the analyzed trait")
parser.add_argument("--outdir", help = "Directory where to output results")
args = parser.parse_args()

#set the timer
start_time1 = time.time()

#stop the execution if all the parameters are not there
if args.input is None:
    print ("Please provide the file that is output of CHEERS_normalize.py.")
    sys.exit()

if (args.ld is None and args.snp_list is None):
    print ("Please provide path to LD or SNP-list.")
    sys.exit(1)

if (args.ld is not None and args.snp_list is not None):
    print ("Please provide only path to LD or SNP-list.")
    sys.exit(1)

if args.trait is None:
    print("Please provide the trait name for the output files.")
    sys.exit(1)

if args.outdir is None:
    print("Please provide the path to directory to output files.")
    sys.exit(1)

#set variables
samplesList = []

'''
open file containing normalized counts per region. 
output of the code: CHEERS_normalize.py
'''

#load data

conn = sqlite3.connect('cheers_tmp.db')

#@profile
def load_sample_data(readin, samplesList):
    df = pd.read_csv(readin, sep='\t')
    rowCount = len(df.index)
    for i in range(3, len(df.columns)):
        sample = df.columns[i]
        samplesList.append(sample)
        # set ranking for each value for each sample
        df['rank_' + str(sample)] = df[sample].rank(method='min').astype(int)
        df.rename(columns={sample: 'value_' + str(sample)}, inplace=True)
    melt = pd.wide_to_long(df, ['value', 'rank'], i=['chr', 'start', 'end'], j='sample', sep='_', suffix='\\w+').reset_index() 
    melt.to_sql('peaks', conn, index=False, if_exists='replace')
    return rowCount

#@profile
def load_snps(snps):
    df = pd.read_csv(snps, sep='\t', header=None, names=['name', 'snp_chr', 'pos'])
    df.to_sql('snps', conn, index=False, if_exists='replace')

sql_query = '''
       select peaks.*, snps.name, snps.snp_chr, snps.pos
       from peaks join snps on pos between start and end 
       and peaks.chr = snps.snp_chr
       '''
    
#@profile
def merge_peaks_with_snps(sql_query, conn):
    merged = pd.read_sql_query(sql_query, conn)
    merged.drop(columns=['snp_chr'], inplace=True)
    return merged

print("loading peak data...")
N = load_sample_data(args.input, samplesList)
print("loading snp data...")
load_snps(args.snp_list)
print("merging data...")
overlap = merge_peaks_with_snps(sql_query, conn)

print("writing snp output")
snp_out_name = str(args.outdir) + str(args.trait) + '_SNPsOverlappingPeaks.txt'
SNPsOverlappingPeaks = overlap.groupby(['name', 'chr', 'pos', 'start', 'end','sample'])['rank'].min().unstack(fill_value='Null').reset_index()
SNPsOverlappingPeaks.to_csv(snp_out_name, index=False)

print("writing unique peaks output")
uniquePeaks = overlap.groupby(['chr', 'start', 'end', 'sample'])['rank'].min().unstack(fill_value='Null').reset_index()
unique_out_file = str(args.outdir) + str(args.trait) + '_uniquePeaks.txt'
uniquePeaks.to_csv(unique_out_file, index=False)

print("calculating means and p-values")
samplesList = [{sample: {'observedMean': uniquePeaks[sample].mean()}} for sample in samplesList]

n = len(uniquePeaks.index)

mean_sd = math.sqrt((N**2-1)/(12*n))
mean_mean = (1+N)/2

for sample in samplesList:
    for key, value in sample.iteritems():
        value['pValue'] = 1-norm.cdf(value['observedMean'], loc=mean_mean, scale=mean_sd)

'''
create the txt file with the p-values
'''
pValueName = str(args.outdir) + str(args.trait) + '_disease_enrichment_pValues.txt'
with open(pValueName, 'w') as f:
    for sample in samplesList:
        for key, value in sample.iteritems():
            f.write(key + '\t' + str(value['pValue']) + '\n')

'''
create the txt file with the mean ranks
'''
meanName = str(args.outdir) + str(args.trait) + '_disease_enrichment_observedMeanRank.txt'
with open(meanName, 'w') as f:
    for sample in samplesList:
        for key, value in sample.iteritems():
            f.write(key + '\t' + str(value['observedMean']) + '\n')

'''
output in the log file
'''

end_time1 = time.time()
running_time = (end_time1 - start_time1)

logfileName = str(args.outdir) + str(args.trait) + ".log"
with open(logfileName, "w") as logfile:
    print >> logfile, 'Total number of peaks\t%s' % (str(N))
    print >> logfile, 'Number of overlapping peaks\t%s' % (str(n))
    print >> logfile, 'Number of SNPs overlapping peaks\t%s' % (str(len(SNPsOverlappingPeaks.index)))
    print >> logfile, 'Distribution mean\t%s' %  (str(mean_mean))
    print >> logfile, 'Distribution sd\t%s' %  (str(mean_sd))
    print >> logfile, 'Running time in seconds\t%s' % (running_time)

