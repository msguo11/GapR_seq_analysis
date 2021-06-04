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

sns.set_style('white')
plt.rcParams['xtick.labelsize']=15
plt.rcParams['ytick.labelsize']=15

#functions for reading a fasta and calculating the AT content
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

def readGenome(a, skip = 0):
    genomeFile = open(a, 'r')
    out = ''
    if skip != 0:
        for i in range(0,skip,1):
                genomeFile.readline()
    line = genomeFile.readline()
    while line != '':
        out = out + line[:-1]
        line = genomeFile.readline()
    return out

def readCDSMG1655(annoteFile, skip = 0):
    a = open(annoteFile, 'r')
    gtype, start, end, strand, funct, bNum, gene = [], [], [], [], [], [], []
    for i in range(0,skip):
        a.readline()
    line = a.readline()
    while line != '':
        if regex.findall('CDS', line):
            z = line.split('\t')
            b = z[8].split('ID=')
            c = b[1].split(':')[0]
            gtype.append(z[2])
            start.append(z[3])
            end.append(z[4])
            strand.append(z[6])
            if regex.findall('product', line):
                zz = line.split('product=')[1]
                funct.append(zz.split(';')[0])
            else:
                funct.append('n/a')
            y = line.split('locus_tag=')[1]
            bNum.append(y.split(';')[0])
            gene.append(c.split('\"')[0])
        line = a.readline()
    out = np.array([gtype, start, end, strand, funct, bNum, gene])
    out = pd.DataFrame(out).transpose()
    out.columns = ['gtype', 'start', 'end', 'strand', 'function', 'bNum', 'geneName']
    return out

genome_folder = './data/'
ecolifasta=loading_fasta_gbk(genome_folder + 'NC000913_2.fasta','fasta')
ecoligbk=loading_fasta_gbk(genome_folder + 'NC_000913_2.gbk','genbank')
genome=reading_gbk_new(ecoligbk,['CDS','tRNA','rRNA','ncRNA'])

def downloadGEO(filename):
    if ~os.path.exists('./data/{}'.format(filename)):
        print('downloading {}.gz from GEO'.format(filename))
        f = filename.split('_')
        os.system("wget -P ./data/ https://ftp.ncbi.nlm.nih.gov/geo/samples/{}nnn/{}/suppl/{}.gz".format(f[0][:-3],f[0],filename))
        print('unzipping {}.gz'.format(filename))
        os.system('gzip -d ./data/{}.gz'.format(filename))

def loadChipData():
    chip = pd.DataFrame()
    
    downloadGEO('GSM4628313_D18-11475-3531L_norm.wig')
    chip['high']=pd.read_csv('./data/GSM4628313_D18-11475-3531L_norm.wig', sep = '\t', header = None, skiprows=2)[1]

    downloadGEO('GSM4628314_D19-11570-4278G_MG1655_norm.wig')
    chip['high_rep2']=pd.read_csv('./data/GSM4628314_D19-11570-4278G_MG1655_norm.wig', sep = '\t', header = None, skiprows=2)[1]
    
    downloadGEO('GSM4628315_D18-11479-3531L_norm.wig')
    chip['Rif_high']=pd.read_csv('./data/GSM4628315_D18-11479-3531L_norm.wig', sep = '\t', header = None, skiprows=2)[1]
 
    downloadGEO('GSM4628311_D19-5504-3883G_norm.wig')
    chip['wtH']=pd.read_csv('./data/GSM4628311_D19-5504-3883G_norm.wig', sep = '\t', header = None, skiprows=2)[1]

    downloadGEO('GSM4628312_D19-11573-4278G_MG1655_norm.wig')
    chip['wtH_rep2']=pd.read_csv('./data/GSM4628312_D19-11573-4278G_MG1655_norm.wig', sep = '\t', header = None, skiprows=2)[1]

    downloadGEO('GSM4989042_D20-5423-4700M_norm.wig')
    chip['76_rep1']=pd.read_csv('./data/GSM4989042_D20-5423-4700M_norm.wig', sep = '\t', header = None, skiprows=2)[1]
    
    downloadGEO('GSM4989043_D20-5424-4700M_norm.wig')
    chip['76_rep2']=pd.read_csv('./data/GSM4989043_D20-5424-4700M_norm.wig', sep = '\t', header = None, skiprows=2)[1]
    
    return chip

def loadRNAseqData():
    #import RNA-seq wigs
    downloadGEO('GSM4628309_D19-11574-4278G_MG1655_norm_fw.wig')
    RNAseqf_me = pd.read_csv('./data/GSM4628309_D19-11574-4278G_MG1655_norm_fw.wig',sep = '\t',header = None,skiprows=2, index_col=0)
    downloadGEO('GGSM4628309_D19-11574-4278G_MG1655_norm_rv.wig')
    RNAseqr_me = pd.read_csv('./data/GSM4628309_D19-11574-4278G_MG1655_norm_rv.wig',sep = '\t',header = None,skiprows=2, index_col=0)
    RNAseq_me = RNAseqf_me.reindex(RNAseqf_me.index,fill_value=0)
    RNAseq_me['rev'] = RNAseqr_me.reindex(RNAseqf_me.index,fill_value=0)
    RNAseq_me.columns = ['fwd','rev']
    
    return RNAseq_me

##calculating ATcontent
def ATcontent(start, end):
    from Bio.Seq import Seq, MutableSeq 
    from Bio.SeqUtils import GC  
    content=100-(GC(ecolifasta.seq[start:end]))
    return content

### code a sliding window to record AT content
def sliding_window(window_length):
    sliding_array=np.zeros([1,len(ecolifasta.seq)])
    for i in range(int(window_length/2), len(ecolifasta.seq)-int(window_length/2)):
        start=i-int(window_length/2)
        end=i+int(window_length/2)
        sliding_array[0][i]=ATcontent(start, end)
    return sliding_array

#mask rRNA loci
def maskrRNA(chipFile):
    chipFile[4034004:4038929] = np.nan
    chipFile[4166428:4170080] = np.nan
    chipFile[3939350:3945307] = np.nan
    chipFile[3421216:3427258] = np.nan
    chipFile[4203834:4221970] = np.nan
    chipFile[2723768:2729041] = np.nan
    chipFile[223408:229167] = np.nan
    return chipFile

#Guo and Haakonsen et al caulo data
def loading_fasta_gbk(file_name,typeoffile):
    """reads either fasta or gbk files, file type needs to be given as 'fasta' or 'genbank' """
    loaded=SeqIO.read(file_name, typeoffile)
    return loaded

caulobactergbk=loading_fasta_gbk('./data/' + 'NC_011916.gbk','genbank')
caulobacterfasta=loading_fasta_gbk('./data/' + 'NC_011916.fna','fasta')

def gaussian_smooth(data_array, index_step, sigma):
    mu = 0
    bins = np.arange(-4*sigma, 4*sigma, index_step, dtype=np.float32)
    gaussian = index_step*1/(sigma * np.sqrt(2 * np.pi))*np.exp( - (bins - mu)**2 / (2 * sigma**2) )
    return signal.convolve(data_array,gaussian,mode='same')

def initialize_rpm(data,index_step,length):
    
    smoothed_data = np.zeros([1, length], dtype=np.float32)
    smoothed_data[0][:] = gaussian_smooth(data[0],index_step,5*index_step)
    rpm = np.zeros([1, length])
    normalization_factor=float(sum(smoothed_data[0][:]))/float(1000000)
    for i in range(0,length):
        rpm[0][i]=float(smoothed_data[0][i])/normalization_factor
        
    return rpm

### function for loading data
def loading_chip_data(name):
    coordinate=np.loadtxt(name, dtype=np.int32, delimiter='\t', skiprows=2)
    chip_data=np.zeros([1,len(caulobacterfasta.seq)])
    ind=coordinate[:,0]<len(caulobacterfasta.seq)
    chip_data[0][coordinate[ind,0]-1]=coordinate[ind,1]
    chip_rpm=initialize_rpm(chip_data, 10, len(caulobacterfasta.seq))
    return chip_rpm

def loadCauloChipData():
    #negative control
    downloadGEO('GSM2690549_neg_chip.wig')
    chip_neg=loading_chip_data('./data/GSM2690549_neg_chip.wig')
    #GapR-3xFLAG
    downloadGEO('GSM2690550_GapR_chip.wig')
    chip_356=loading_chip_data('./data/GSM2690550_GapR_chip.wig')

    caulo = pd.DataFrame(chip_neg[0])
    caulo['356'] = chip_356[0]
    caulo.columns = ['wt','gapR']
    
    return caulo

##calculating ATcontent, caulos
def ATcontent_caulo(start, end):
    from Bio.Seq import Seq, MutableSeq 
    from Bio.SeqUtils import GC  
    content=100-(GC(caulobacterfasta.seq[start:end]))
    return content

### code a sliding window to record AT content
def sliding_window_caulo(window_length):
    sliding_array=np.zeros([1,len(caulobacterfasta.seq)])
    for i in range(int(window_length/2), len(caulobacterfasta.seq)-int(window_length/2)):
        start=i-int(window_length/2)
        end=i+int(window_length/2)
        sliding_array[0][i]=ATcontent_caulo(start, end)
    return sliding_array

#recover 200 bp around most enriched GapR-3xFLAG regions
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
        cmax = chip_data.loc[i[0]:i[1]].idxmax()
#        print(cmax, i[0],i[1])
        chipMax.append(cmax)
        sequence.append(str(genomeSequence.seq[cmax-half_window:cmax+half_window]))
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
            cmin = chip_data.loc[i[0]:i[1]].idxmin()
            chipMin.append(cmin)
            sequence.append(str(genomeSequence.seq[cmin-half_window:cmin+half_window]))
    if saveFile != False:
        a = open(saveFile, 'w')
        for i in range(len(sequence)):
            a.write('>loc{}\n'.format(chipMin[i]))
            a.write(sequence[i] + '\n')
        a.close()
    
    return sequence

#compress intergenic regions
def readCDSMG1655(annoteFile, skip = 0):
    a = open(annoteFile, 'r')
    gtype, start, end, strand, funct, bNum, gene = [], [], [], [], [], [], []
    for i in range(0,skip):
        a.readline()
    line = a.readline()
    while line != '':
        if regex.findall('CDS', line):
            z = line.split('\t')
            b = z[8].split('ID=')
            c = b[1].split(':')[0]
            gtype.append(z[2])
            start.append(z[3])
            end.append(z[4])
            strand.append(z[6])
            if regex.findall('product', line):
                zz = line.split('product=')[1]
                funct.append(zz.split(';')[0])
            else:
                funct.append('n/a')
            y = line.split('locus_tag=')[1]
            bNum.append(y.split(';')[0])
            gene.append(c.split('\"')[0])
        line = a.readline()
    out = np.array([gtype, start, end, strand, funct, bNum, gene])
    out = pd.DataFrame(out).transpose()
    out.columns = ['gtype', 'start', 'end', 'strand', 'function', 'bNum', 'geneName']
    return out

MG1655annoteFile = './data/gffEditedNoDup.txt'
annote = readCDSMG1655(MG1655annoteFile, 3)

def returnAnnote(RNAseq):
    meanS = []
    sumS = []
    for ind, ann in annote.iterrows():
        if ann.strand == '+':
            meanS.append(np.mean(RNAseq['fwd'][int(ann.start)-1:int(ann.end)]))
            sumS.append(np.sum(RNAseq['fwd'][int(ann.start)-1:int(ann.end)]))
        else:
            meanS.append(np.mean(RNAseq['rev'][int(ann.start)-1:int(ann.end)]))
            sumS.append(np.sum(RNAseq['rev'][int(ann.start)-1:int(ann.end)]))

    annote['newRPK'] = meanS / sum(sumS) * 1000000 * 1000 #convert to rpkm
    annote['start']=annote['start'].astype(int)
    annote['end']=annote['end'].astype(int)
    
    return annote

def next_N_sum(data,N):
    """Function to populate an array that contains the sliding window sum of the N previous bps"""
    data_len = len(data)
    cumsum = data.astype(np.float64).cumsum().values
    next_N_sum = np.zeros([1, data_len], dtype=np.float64)
    next_N_sum[0,:N] = cumsum[:N]
    next_N_sum[0,N:] = cumsum[N:] - cumsum[:data_len-N]
    return next_N_sum