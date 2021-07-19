import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import regex
import scipy.stats as stats
import scipy.optimize as optimize
import scipy.signal as signal
from scipy import cluster
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq 
from Bio.SeqUtils import GC  
import random

sns.set_style('white')
plt.rcParams['xtick.labelsize']=15
plt.rcParams['ytick.labelsize']=15

def loading_fasta_gbk(file_name,typeoffile):
    """reads either fasta or gbk files, file type needs to be given as 'fasta' or 'genbank' """
    loaded=SeqIO.read(file_name, typeoffile)
    return loaded

class Genome:
    def __init__(self, genome_list,genome_annotation, start, end, strand, length):
        self.name=genome_list #list with every gene name such as CCNA_00001
        self.annotation=genome_annotation # gene annotation if there is one, if none stores NA
        self.start=start # stores translational start position for each gene
        self.end=end #stores end position of each gene
        self.strand=strand # + or - strand (+1 or -1)
        self.length=length # length of gene
        
def reading_gbk_new(gbk_file, features_to_extract):
    """function that will load from the gbk file: the start, end, strand and length of gene as well as the name and annotated name/function. 
    Returns one array and 2 lists """
    
    genome_gene=[]
    genome_gene_name=[]
    
    genome_start=[]
    genome_end=[]
    genome_strand=[]
    genome_length=[]
    for i in range(0,len(gbk_file.features)):
        isfeature=False
        for j in range(len(features_to_extract)):
            if gbk_file.features[i].type == features_to_extract[j]:
                isfeature=True
            
        if isfeature==True:
            
            genome_gene.append(gbk_file.features[i].qualifiers['locus_tag'][0])

            if 'product' in gbk_file.features[i].qualifiers:
                genome_gene_name.append(gbk_file.features[i].qualifiers['product'][0]) 
            else:
                genome_gene_name.append('NA')

            if gbk_file.features[i].location.strand < 0 :
                genome_start.append(gbk_file.features[i].location.end)
                genome_end.append(gbk_file.features[i].location.start)
                genome_strand.append(-1)
                genome_length.append(abs(gbk_file.features[i].location.end-gbk_file.features[i].location.start)+1)
            else:
                genome_start.append(gbk_file.features[i].location.start)
                genome_end.append(gbk_file.features[i].location.end)
                genome_strand.append(1)
                genome_length.append(abs(gbk_file.features[i].location.end-gbk_file.features[i].location.start)+1)

    genome = Genome(genome_gene,genome_gene_name,genome_start,genome_end,genome_strand,genome_length) 
    
    return genome    

genome_folder = './data/scer/'
genomefasta = {}
for i in range(1,10):
    x = loading_fasta_gbk(genome_folder + 'chr0{}.fsa'.format(i),'fasta')
    genomefasta[x.name] = x
for i in range(10,17):
    x = loading_fasta_gbk(genome_folder + 'chr{}.fsa'.format(i),'fasta')
    genomefasta[x.name] = x

l = {}
for c in ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI']:
    chrom = 'chr'+c
    l[chrom]=len(genomefasta[chrom].seq)

la = {'chrI': 0,
 'chrII': 230218,
 'chrIII': 1043402,
 'chrIV': 1360022,
 'chrV': 2891955,
 'chrVI': 3468829,
 'chrVII': 3738990,
 'chrVIII': 4829930,
 'chrIX': 5392573,
 'chrX': 5832461,
 'chrXI': 6578212,
 'chrXII': 7245028,
 'chrXIII': 8323205,
 'chrXIV': 9247636,
 'chrXV': 10031969,
 'chrXVI': 11123260}

##calculating ATcontent
def ATcontent(genome, start, end):
    content=100-(GC(genome.seq[start:end]))
    return content

### code a sliding window to record AT content
def sliding_window(genome,window_length):
    sliding_array=np.zeros([1,len(genome.seq)])
    for i in range(int(window_length/2), len(genome.seq)-int(window_length/2)):
        start=i-int(window_length/2)
        end=i+int(window_length/2)
        sliding_array[0][i]=ATcontent(genome,start, end)
    return sliding_array

def formatGenomeDict(genomedict,genomefasta):
    l = {}
    for c in ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI']:
        chrom = 'chr'+c
        l[chrom]=len(genomefasta[chrom].seq)
    genomedict['chrI']['fullstart']=genomedict['chrI']['start']
    genomedict['chrI']['fullend']=genomedict['chrI']['end']

    chrom='chrII'
    genomedict[chrom]['fullstart']=genomedict[chrom]['start']+l['chrI']
    genomedict[chrom]['fullend']=genomedict[chrom]['end']+l['chrI']

    chrom='chrIII'
    genomedict[chrom]['fullstart']=genomedict[chrom]['start']+l['chrI']+l['chrII']
    genomedict[chrom]['fullend']=genomedict[chrom]['end']+l['chrI']+l['chrII']

    chrom='chrIV'
    genomedict[chrom]['fullstart']=genomedict[chrom]['start']+l['chrI']+l['chrII']+l['chrIII']
    genomedict[chrom]['fullend']=genomedict[chrom]['end']+l['chrI']+l['chrII']+l['chrIII']

    chrom='chrV'
    genomedict[chrom]['fullstart']=genomedict[chrom]['start']+l['chrI']+l['chrII']+l['chrIII']+l['chrIV']
    genomedict[chrom]['fullend']=genomedict[chrom]['end']+l['chrI']+l['chrII']+l['chrIII']+l['chrIV']

    chrom='chrVI'
    genomedict[chrom]['fullstart']=genomedict[chrom]['start']+l['chrI']+l['chrII']+l['chrIII']+l['chrIV']+l['chrV']
    genomedict[chrom]['fullend']=genomedict[chrom]['end']+l['chrI']+l['chrII']+l['chrIII']+l['chrIV']+l['chrV']

    chrom='chrVII'
    genomedict[chrom]['fullstart']=genomedict[chrom]['start']+l['chrI']+l['chrII']+l['chrIII']+l['chrIV']+l['chrV']+l['chrVI']
    genomedict[chrom]['fullend']=genomedict[chrom]['end']+l['chrI']+l['chrII']+l['chrIII']+l['chrIV']+l['chrV']+l['chrVI']

    chrom='chrVIII'
    genomedict[chrom]['fullstart']=genomedict[chrom]['start']+l['chrI']+l['chrII']+l['chrIII']+l['chrIV']+l['chrV']+l['chrVI']+l['chrVII']
    genomedict[chrom]['fullend']=genomedict[chrom]['end']+l['chrI']+l['chrII']+l['chrIII']+l['chrIV']+l['chrV']+l['chrVI']+l['chrVII']

    chrom='chrIX'
    genomedict[chrom]['fullstart']=genomedict[chrom]['start']+l['chrI']+l['chrII']+l['chrIII']+l['chrIV']+l['chrV']+l['chrVI']+l['chrVII']+l['chrVIII']
    genomedict[chrom]['fullend']=genomedict[chrom]['end']+l['chrI']+l['chrII']+l['chrIII']+l['chrIV']+l['chrV']+l['chrVI']+l['chrVII']+l['chrVIII']

    chrom='chrX'
    genomedict[chrom]['fullstart']=genomedict[chrom]['start']+l['chrI']+l['chrII']+l['chrIII']+l['chrIV']+l['chrV']+l['chrVI']+l['chrVII']+l['chrVIII']+l['chrIX']
    genomedict[chrom]['fullend']=genomedict[chrom]['end']+l['chrI']+l['chrII']+l['chrIII']+l['chrIV']+l['chrV']+l['chrVI']+l['chrVII']+l['chrVIII']+l['chrIX']

    chrom='chrXI'
    genomedict[chrom]['fullstart']=genomedict[chrom]['start']+l['chrI']+l['chrII']+l['chrIII']+l['chrIV']+l['chrV']+l['chrVI']+l['chrVII']+l['chrVIII']+l['chrIX']+l['chrX']
    genomedict[chrom]['fullend']=genomedict[chrom]['end']+l['chrI']+l['chrII']+l['chrIII']+l['chrIV']+l['chrV']+l['chrVI']+l['chrVII']+l['chrVIII']+l['chrIX']+l['chrX']

    chrom='chrXII'
    genomedict[chrom]['fullstart']=genomedict[chrom]['start']+l['chrI']+l['chrII']+l['chrIII']+l['chrIV']+l['chrV']+l['chrVI']+l['chrVII']+l['chrVIII']+l['chrIX']+l['chrX']+l['chrXI']
    genomedict[chrom]['fullend']=genomedict[chrom]['end']+l['chrI']+l['chrII']+l['chrIII']+l['chrIV']+l['chrV']+l['chrVI']+l['chrVII']+l['chrVIII']+l['chrIX']+l['chrX']+l['chrXI']

    chrom='chrXIII'
    genomedict[chrom]['fullstart']=genomedict[chrom]['start']+l['chrI']+l['chrII']+l['chrIII']+l['chrIV']+l['chrV']+l['chrVI']+l['chrVII']+l['chrVIII']+l['chrIX']+l['chrX']+l['chrXI']+l['chrXII']
    genomedict[chrom]['fullend']=genomedict[chrom]['end']+l['chrI']+l['chrII']+l['chrIII']+l['chrIV']+l['chrV']+l['chrVI']+l['chrVII']+l['chrVIII']+l['chrIX']+l['chrX']+l['chrXI']+l['chrXII']

    chrom='chrXIV'
    genomedict[chrom]['fullstart']=genomedict[chrom]['start']+l['chrI']+l['chrII']+l['chrIII']+l['chrIV']+l['chrV']+l['chrVI']+l['chrVII']+l['chrVIII']+l['chrIX']+l['chrX']+l['chrXI']+l['chrXII']+l['chrXIII']
    genomedict[chrom]['fullend']=genomedict[chrom]['end']+l['chrI']+l['chrII']+l['chrIII']+l['chrIV']+l['chrV']+l['chrVI']+l['chrVII']+l['chrVIII']+l['chrIX']+l['chrX']+l['chrXI']+l['chrXII']+l['chrXIII']

    chrom='chrXV'
    genomedict[chrom]['fullstart']=genomedict[chrom]['start']+l['chrI']+l['chrII']+l['chrIII']+l['chrIV']+l['chrV']+l['chrVI']+l['chrVII']+l['chrVIII']+l['chrIX']+l['chrX']+l['chrXI']+l['chrXII']+l['chrXIII']+l['chrXIV']
    genomedict[chrom]['fullend']=genomedict[chrom]['end']+l['chrI']+l['chrII']+l['chrIII']+l['chrIV']+l['chrV']+l['chrVI']+l['chrVII']+l['chrVIII']+l['chrIX']+l['chrX']+l['chrXI']+l['chrXII']+l['chrXIII']+l['chrXIV']

    chrom='chrXVI'
    genomedict[chrom]['fullstart']=genomedict[chrom]['start']+l['chrI']+l['chrII']+l['chrIII']+l['chrIV']+l['chrV']+l['chrVI']+l['chrVII']+l['chrVIII']+l['chrIX']+l['chrX']+l['chrXI']+l['chrXII']+l['chrXIII']+l['chrXIV']+l['chrXV']
    genomedict[chrom]['fullend']=genomedict[chrom]['end']+l['chrI']+l['chrII']+l['chrIII']+l['chrIV']+l['chrV']+l['chrVI']+l['chrVII']+l['chrVIII']+l['chrIX']+l['chrX']+l['chrXI']+l['chrXII']+l['chrXIII']+l['chrXIV']+l['chrXV']
        
    return genomedict

def downloadGEO(filename):
    if ~os.path.exists('./data'):
        os.system('mkdir ./data')
    if ~os.path.exists('./data/{}'.format(filename)):
        print('downloading {}.gz from GEO'.format(filename))
        f = filename.split('_')
        os.system("wget -P ./data/ https://ftp.ncbi.nlm.nih.gov/geo/samples/{}nnn/{}/suppl/{}.gz".format(f[0][:-3],f[0],filename))
        print('unzipping {}.gz'.format(filename))
        os.system('gzip -d ./data/{}.gz'.format(filename))
        
def downloadMNase(filename = 'GSM3069971_MNase_YPD30_WT-B_166_40U.sgr'):
    downloadGEO(filename)
    mnase = pd.read_csv('./data/GSM3069971_MNase_YPD30_WT-B_166_40U.sgr',sep = '\t',header = None)
    mnase.columns = ['chr','pos','val'] #the mnase data is 1 nt short

    return mnase

def downloadScc1(filename = 'GSM3668254_Scc1-WT_Tension.calibrated.txt'):
    downloadGEO(filename)
    Scc1 = pd.read_csv('./data/GSM3668254_Scc1-WT_Tension.calibrated.txt',sep=',',index_col=0)
    
    return Scc1

def loadyeastRNAseqData():
    #import RNA-seq wigs
    downloadGEO('GSM5001907_D20-252008_nodup_plus_all.txt')
    RNAseq_pl = pd.read_csv('./data/GSM5001907_D20-252008_nodup_plus_all.txt',sep = ',',index_col=0)
    downloadGEO('GSM5001907_D20-252008_nodup_minus_all.txt')
    RNAseq_min = pd.read_csv('./data/GSM5001907_D20-252008_nodup_minus_all.txt',sep = ',',index_col=0)
    RNAseq_me = RNAseq_pl.copy()
    RNAseq_me['rev'] = RNAseq_min.val_norm
    RNAseq_me = RNAseq_me.drop(columns = 'value')
    RNAseq_me.columns = ['chr','pos','fwd','rev']
    RNAseq_me['merged'] = RNAseq_pl['val_norm'].values.copy()+RNAseq_min['val_norm'].values.copy()
    
    return RNAseq_me

def loadyeastRNAseqaFData():
    #import RNA-seq wigs
    downloadGEO('GSM5001909_D20-252007_nodup_plus_all.txt')
    RNAseq_pl = pd.read_csv('./data/GSM5001909_D20-252007_nodup_plus_all.txt',sep = ',',index_col=0)
    downloadGEO('GSM5001909_D20-252007_nodup_minus_all.txt')
    RNAseq_min = pd.read_csv('./data/GSM5001909_D20-252007_nodup_minus_all.txt',sep = ',',index_col=0)
    RNAseq_me = RNAseq_pl.copy()
    RNAseq_me['rev'] = RNAseq_min.val_norm
    RNAseq_me = RNAseq_me.drop(columns = 'value')
    RNAseq_me.columns = ['chr','pos','fwd','rev']
    RNAseq_me['merged'] = RNAseq_pl['val_norm'].values.copy()+RNAseq_min['val_norm'].values.copy()
    
    return RNAseq_me

def loadyeastGlyRNAseqData():
    #import RNA-seq wigs
    downloadGEO('GSM5001911_D20-3448_plus_all.txt')
    RNAseq_pl = pd.read_csv('./data/GSM5001911_D20-3448_plus_all.txt',sep = ',',index_col=0)
    downloadGEO('GSM5001911_D20-3448_minus_all.txt')
    RNAseq_min = pd.read_csv('./data/GSM5001911_D20-3448_minus_all.txt',sep = ',',index_col=0)
    RNAseq_me = RNAseq_pl.copy()
    RNAseq_me['rev'] = RNAseq_min.val_norm
    RNAseq_me = RNAseq_me.drop(columns = 'value')
    RNAseq_me.columns = ['chr','pos','fwd','rev']
    RNAseq_me['merged'] = RNAseq_pl['val_norm'].values.copy()+RNAseq_min['val_norm'].values.copy()
    
    return RNAseq_me

def loadraffChIP():
    chip = {}
    downloadGEO('GSM5001899_D20-5952_all_nodup.txt')
    yGapR_raff = pd.read_csv('./data/GSM5001899_D20-5952_all_nodup.txt',sep = ',', index_col=0) #gapR chIP, in rep1
    normalization_factor = sum(yGapR_raff[yGapR_raff['chr']!='chrXII'].value.values)/1000000
    yGapR_raff['val_norm_no12'] = yGapR_raff.value/normalization_factor
    yGapR_raff['smooth'] = yGapR_raff.val_norm_no12.rolling(250,center=True).mean()
    
    return yGapR_raff

def loadraff2ChIP():
    chip = {}
    downloadGEO('GSM5001900_D20-5953_all_nodup.txt')
    yGapR_raff = pd.read_csv('./data/GSM5001900_D20-5953_all_nodup.txt',sep = ',', index_col=0) #gapR chIP, in rep1
    normalization_factor = sum(yGapR_raff[yGapR_raff['chr']!='chrXII'].value.values)/1000000
    yGapR_raff['val_norm_no12'] = yGapR_raff.value/normalization_factor
    yGapR_raff['smooth'] = yGapR_raff.val_norm_no12.rolling(250,center=True).mean()
    
    return yGapR_raff

def loadraffnegChIP():
    chip = {}
    downloadGEO('GSM5001905_D20-261001_all_nodup.txt')
    negIP = pd.read_csv('./data/GSM5001905_D20-261001_all_nodup.txt',sep = ',', index_col=0) #neg ctrl
    normalization_factor = sum(negIP[negIP['chr']!='chrXII'].value.values)/1000000
    negIP['val_norm_no12'] = negIP.value/normalization_factor
    negIP['smooth'] = negIP.val_norm_no12.rolling(250,center=True).mean()
    
    return negIP

def loadraffINP():
    downloadGEO('GSM5001903_D20-6528_all_nodup.txt')
    negINP = pd.read_csv('./data/GSM5001903_D20-6528_all_nodup.txt',sep = ',', index_col=0)
    normalization_factor = sum(negINP[negINP['chr']!='chrXII'].value.values)/1000000
    negINP['val_norm_no12'] = negINP.value/normalization_factor
    negINP['smooth'] = negINP.val_norm_no12.rolling(250,center=True).mean()

    return negINP
    
def loadraffaFINP():
    downloadGEO('GSM5001904_D20-6530_all_nodup.txt')
    negINPaF = pd.read_csv('./data/GSM5001904_D20-6530_all_nodup.txt',sep = ',', index_col=0)
    normalization_factor = sum(negINPaF[negINPaF['chr']!='chrXII'].value.values)/1000000
    negINPaF['val_norm_no12'] = negINPaF.value/normalization_factor
    negINPaF['smooth'] = negINPaF.val_norm_no12.rolling(250,center=True).mean()
    
    return negINPaF

def loadChIPFold(IP,negIP):
    chip_fold = IP.copy()
    chip_fold['fold_nolog'] = (IP.smooth+0.01)/(negIP.smooth+0.01)
    chip_fold.drop(columns=['value','val_norm','val_norm_no12','smooth'],inplace=True)
    
    return chip_fold

def loadraffaFChIP():
    chip = {}
    downloadGEO('GSM5001901_D20-5954_all_nodup.txt')
    yGapR_raffaF = pd.read_csv('./data/GSM5001901_D20-5954_all_nodup.txt',sep = ',', index_col=0) #gapR chIP, in rep1
    normalization_factor = sum(yGapR_raffaF[yGapR_raffaF['chr']!='chrXII'].value.values)/1000000
    yGapR_raffaF['val_norm_no12'] = yGapR_raffaF.value/normalization_factor
    yGapR_raffaF['smooth'] = yGapR_raffaF.val_norm_no12.rolling(250,center=True).mean()
    
    return yGapR_raffaF

def loadaFnegChIP():
    chip = {}
    downloadGEO('GSM5001906_D20-261002_all_nodup.txt')
    yGapR_aFneg = pd.read_csv('./data/GSM5001906_D20-261002_all_nodup.txt',sep = ',', index_col=0) #gapR chIP, in rep1
    normalization_factor = sum(yGapR_aFneg[yGapR_aFneg['chr']!='chrXII'].value.values)/1000000
    yGapR_aFneg['val_norm_no12'] = yGapR_aFneg.value/normalization_factor
    yGapR_aFneg['smooth'] = yGapR_aFneg.val_norm_no12.rolling(250,center=True).mean()
    
    return yGapR_aFneg

def loadglyChIP():
    chip = {}
    downloadGEO('GSM4628318_D19_5482_all.txt')
    gly = pd.read_csv('./data/GSM4628318_D19_5482_all.txt',sep = ',', index_col=0) #gapR chIP, in rep1
    normalization_factor = sum(gly[gly['chr']!='chrXII'].value.values)/1000000
    gly['val_norm_no12'] = gly.value/normalization_factor
    gly['smooth'] = gly.val_norm_no12.rolling(250,center=True).mean()
    
    return gly

def loadglynegChIP():
    chip = {}
    downloadGEO('GSM4628316_D19_5480_all.txt')
    neg = pd.read_csv('./data/GSM4628316_D19_5480_all.txt',sep = ',', index_col=0) #gapR chIP, in rep1
    normalization_factor = sum(neg[neg['chr']!='chrXII'].value.values)/1000000
    neg['val_norm_no12'] = neg.value/normalization_factor
    neg['smooth'] = neg.val_norm_no12.rolling(250,center=True).mean()
    
    return neg

half_window = 100

def enrichedRegions(chip_data, cutoff):
    x = chip_data[chip_data > cutoff]
    out = []
    z = x.index[0]
    start = x.index[0]
    for i in range(1,len(x.index)):
        if x.index[i] == z+1:
            z = x.index[i]
        else:
            end = z
            out.append([start,end])
            z = x.index[i]
            start = x.index[i]
    out.append([start,z])
    return out 

def getSequence(chipRegions, chip_data,genomeSequence, saveFile = False):
    chipMax = []
    sequence = []
    for i in chipRegions:
        cmax = chip_data.loc[i[0]:i[1]].val_norm.idxmax()
        chrom = chip_data.loc[cmax].chr
        chrom_pos = chip_data.loc[cmax].pos
#        print(cmax, i[0],i[1])
        chipMax.append(str(cmax)+'_'+chrom)
        sequence.append(str(genomeSequence[chrom].seq[chrom_pos-half_window:chrom_pos+half_window]))
    if saveFile != False:
        a = open(saveFile, 'w')
        for i in range(len(sequence)):
            a.write('>loc{}\n'.format(chipMax[i]))
            a.write(sequence[i] + '\n')
        a.close()
    
    return sequence

def unenrichedRegions(chip_data, cutoff):
    x = chip_data[chip_data < cutoff]
    out = []
    z = x.index[0]
    start = x.index[0]
    for i in range(1,len(x.index)):
        if x.index[i] == z+1:
            z = x.index[i]
        else:
            end = z
            out.append([start,end])
            z = x.index[i]
            start = x.index[i]
    out.append([start,z])
    return out 

def getuSequence(chipRegions, chip_data,genomeSequence, saveFile = False):
    chipMin = []
    sequence = []
    for i in chipRegions:
        if (i[1] - i[0]) > 100:
            cmin = chip_data.loc[i[0]:i[1]].val_norm.idxmin()
            chrom = chip_data.loc[cmin].chr
            chrom_pos = chip_data.loc[cmin].pos
            if len(str(genomeSequence[chrom].seq[chrom_pos-half_window:chrom_pos+half_window])) > 100:
                chipMin.append(str(cmin)+'_'+chrom)
                sequence.append(str(genomeSequence[chrom].seq[chrom_pos-half_window:chrom_pos+half_window]))
    if saveFile != False:
        a = open(saveFile, 'w')
        for i in range(len(sequence)):
            a.write('>loc{}\n'.format(chipMin[i]))
            a.write(sequence[i] + '\n')
        a.close()
    
    return sequence

def extract_CEN(gbk_file, features_to_extract):
    """function that will load from the gbk file: the start, end, strand and length of gene as well as the name and annotated name/function. 
    Returns one array and 2 lists """
    
    genome_gene=[]
    genome_gene_name=[]
    
    genome_start=[]
    genome_end=[]
    genome_strand=[]
    genome_length=[]
    for i in range(0,len(gbk_file.features)):
        isfeature=False
        for j in range(len(features_to_extract)):
            if gbk_file.features[i].type == features_to_extract[j]:
                isfeature=True
            
        if isfeature==True:
            
            genome_gene.append(gbk_file.features[i].qualifiers['note'][0].split("; ")[0])

            if 'note' in gbk_file.features[i].qualifiers:
                genome_gene_name.append(gbk_file.features[i].qualifiers['note'][0].split("; ")[-1]) 
            else:
                genome_gene_name.append('NA')

            if gbk_file.features[i].location.strand < 0 :
                genome_start.append(gbk_file.features[i].location.end)
                genome_end.append(gbk_file.features[i].location.start)
                genome_strand.append(-1)
                genome_length.append(abs(gbk_file.features[i].location.end-gbk_file.features[i].location.start)+1)
            else:
                genome_start.append(gbk_file.features[i].location.start)
                genome_end.append(gbk_file.features[i].location.end)
                genome_strand.append(1)
                genome_length.append(abs(gbk_file.features[i].location.end-gbk_file.features[i].location.start)+1)

    genome = Genome(genome_gene,genome_gene_name,genome_start,genome_end,genome_strand,genome_length) 
    
    return genome   

def extract_ars(gbk_file, features_to_extract):
    """function that will load from the gbk file: the start, end, strand and length of gene as well as the name and annotated name/function. 
    Returns one array and 2 lists """
    
    genome_gene=[]
    genome_gene_name=[]
    
    genome_start=[]
    genome_end=[]
    genome_strand=[]
    genome_length=[]
    for i in range(0,len(gbk_file.features)):
        isfeature=False
        for j in range(len(features_to_extract)):
            if gbk_file.features[i].type == features_to_extract[j]:
                isfeature=True
            
        if isfeature==True:
            
            genome_gene.append(gbk_file.features[i].qualifiers['note'][0].split("; ")[0])

            if 'note' in gbk_file.features[i].qualifiers:
                genome_gene_name.append(gbk_file.features[i].qualifiers['note'][0].split("; ")[-1]) 
            else:
                genome_gene_name.append('NA')

            if gbk_file.features[i].location.strand < 0 :
                genome_start.append(gbk_file.features[i].location.end)
                genome_end.append(gbk_file.features[i].location.start)
                genome_strand.append(-1)
                genome_length.append(abs(gbk_file.features[i].location.end-gbk_file.features[i].location.start)+1)
            else:
                genome_start.append(gbk_file.features[i].location.start)
                genome_end.append(gbk_file.features[i].location.end)
                genome_strand.append(1)
                genome_length.append(abs(gbk_file.features[i].location.end-gbk_file.features[i].location.start)+1)

    genome = Genome(genome_gene,genome_gene_name,genome_start,genome_end,genome_strand,genome_length) 
    
    return genome   

chr_lengths = pd.read_csv(genome_folder + 'scer.genome',sep = '\t',header=None)
chr_lengths.columns = ['chromosome','length']

def calccumsum(chip):
    gapR_fold_nolog_cumsum = {}
    for chrom in chr_lengths.chromosome:
        gapR_fold_nolog_cumsum[chrom] = np.cumsum(chip[chip.chr==chrom].fold_nolog.values)
        
    return gapR_fold_nolog_cumsum

def calccumsumtelos(chip):
    gapR_fold_nolog_cumsum = {}
    for chrom in chr_lengths.chromosome:
        gapR_fold_nolog_cumsum[chrom] = np.cumsum(chip[chip.chr==chrom].fold_nolog_nosmooth.values)
        
    return gapR_fold_nolog_cumsum