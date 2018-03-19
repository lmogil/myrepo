##############
######Lauren S. Mogil
#####python3.6 
'''For this challenge, you are asked to prototype a variant annotation tool. We will provide you with
a VCF file, and you will create a small software program to output a table annotating each
variant in the file. Each variant must be annotated with the following pieces of information:
1. Type of variation (Substitution, Insertion, Silent, Intergenic, etc.) If there are multiple
possibilities, annotate with the most deleterious possibility.
2. Depth of sequence coverage at the site of variation.
3. Number of reads supporting the variant.
4. Percentage of reads supporting the variant versus those supporting reference reads.
5. Allele frequency of variant from Broad Institute ExAC Project API
(API documentation is available here: http://exac.hms.harvard.edu/)
6. Additional optional information from ExAC that you feel might be relevant.
For this project please upload all relevant code (written in whatever language you like) along
with the annotated VCF file to a Github account and provide the link to the below email address.
Please note that work will be assessed based on quality of code and documentation more-so
than the annotation.'''

#used ExAC release 0.3.1 and downloaded ExAC.r0.3.1.sites.vep.vcf and legacy-exacv1_downloads-release0.3.1-manuscript_data-ExAC.r0.3.1.sites.vep.gene.table

import gzip
import re
import sys
import argparse
import os
from collections import defaultdict

##user input 
#input directory and appropriate files to be used to annotate variants
def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Script to filter imputed VCF')
    parser.add_argument('-i', '--inputdir',
                        help='directory containing VCF file',
                        required='True'
                        )
    parser.add_argument('-c', '--chl',
                        help='challenge_vcf_file',
			type=str,
                        required='True'
                        )
    parser.add_argument('-e', '--exac',
			help='exac vcf file from broad',
			type=str,
			required='True'
			)
    parser.add_argument('-g', '--gene',
                        help='gene file from exac project',
                        type=str,
                        required='True')
    return parser.parse_args(args)
    
#retrieve command line arguments
args = check_arg(sys.argv[1:])
chrpath = args.inputdir
chlfile = args.chl
exacfile = args.exac
genefile = args.gene


#set files 
chrfile = chrpath +  chlfile  
exac_vcf = chrpath +  exacfile
gene_info= chrpath +  genefile

#find genes associated with variants using file from ExAC 
#make a dictionary of variant ids and genes to be used later when formating annotation file
#gene_info='/Users/laurenmogil/Downloads/legacy-exacv1_downloads-release0.3.1-manuscript_data-ExAC.r0.3.1.sites.vep.gene.table'
#gene_info is the file that contains the variants and genes associated along with other information
genes={}
for line in open(gene_info):
    if(line.startswith('CHROM')):	
    	find=line.strip().split('\t').index('Gene')
    	continue
    arr = line.strip().split('\t')
    (chr, pos) = arr[0:2]
    gene=arr[find]
    var=chr+':'+pos
    genes[var] = gene


#exac_vcf='/Users/laurenmogil/Downloads/ExAC.r0.3.1.sites.vep.vcf'
#initialize dictionary key=chr:pos value= the set of relevant info from the ExAC VCF
exac_info = defaultdict(set)
for line in gzip.open(exac_vcf):
    if(line.startswith('##')): #skip over header
            continue
    arr = line.strip().split()
    if(line.startswith('#CHROM')): #only one line should match #CHROM
        continue
    (chr, pos, id, ref, alt, qual, filter, info) = arr[0:8] 
    infos=info.split(';') #split info by semi-colon 
    AFs=[item for item in infos if re.search('AF=',item)] #find the allele freq. in info list
    af=AFs[0].strip().split("=")[1] #pull allele freq. value
    #find clinvar id if it check if it exists and add value if not add place holder
    cvid=[item for item in infos if re.search('clinvar_measureset_id',item)]
    if len(cvid)==0:
    	cv="---"
    else:
    	cv=cvid[0].strip().split("=")[1]
    exac2add=alt+' '+af+' '+cv
    vars=chr+':'+pos #make variant id chr:pos
    exac_info[vars].add(exac2add) #add information to dictionary 
    

#chrfile= '/Users/laurenmogil/Downloads/Challenge_data.vcf'
#open a file to output annotation 
outfile = open(chrpath+"vcf_challenge_annotation_use.txt","w")
outfile.write('Variant_ID'+'\t'+'Type_vars'+'\t'+'Seq_depth'+'\t'+'reads_vars'+'\t'+'percent_reads_support_vars'+'\t'+'AF_challenge'+'\t'+'AF_ExAC'+'\t'+'Clinvar_id'+'\t'+'Gene'+'\n')


for line in open(chrfile):
    if(line.startswith('##')): #skip over header
            continue
    arr = line.strip().split()
    if(line.startswith('#CHROM')): #only one line should match #CHROM
        continue
    (chr, pos, id, ref, alt, qual, filter, info, format) = arr[0:9] #define columns
    infos=info.split(";") #split info column by semi-colon 
    #find relevant column information from INFO column of VCF and define
    TYPE=[item for item in infos if re.search('TYPE=',item)]
    RO=[item for item in infos if re.search('RO=',item)]
    AO=[item for item in infos if re.search('AO=',item)]
    DP=[item for item in infos if re.search('DP=',item)]
    AF=[item for item in infos if re.search('AF=',item)]
    #pull value of variable from above
    typ=TYPE[0].split("=")[1]
    altcount=AO[0].split("=")[1]
    refcount=RO[0].split("=")[1]
    depth=DP[0].split("=")[1]
    af=AF[0].split("=")[1]
    vrs=chr+':'+pos #create variant id
     #if there's more than one type of mutation set as the most deleterious
    #next set of if statements are to rank by deleterious mutation and corresponding values 
    #pull the altcount (number of reads supporting vars) for that variant
    #pull allele frequency in challenge VCF file
    #if theres only one value for type of mutation keep that and proceed
    #if information missing or unavailable added '---' or 'NA'
    if len(typ.strip().split(',')) and len(altcount.strip().split(',')) <2:
    	typ=typ
    	altcount1=altcount
    	af1=af
    #if theres more than one type of mutation use most deleterious
    else:
    	if(bool(re.search('del',typ)) == True):
    		ind=typ.strip().split(',').index("del")
    		typ='del'
    		altcount1=altcount.strip().split(',')[ind]
    		af1=af.strip().split(',')[ind]
    	elif(bool(re.search('ins',typ)) == True):
    		ind=typ.strip().split(',').index("ins")
    		typ='ins'
    		altcount1=altcount.strip().split(',')[ind]
    		af1=af.strip().split(',')[ind]
    	elif(bool(re.search('mnp',typ)) == True):
    		ind=typ.strip().split(',').index("mnp")
    		typ='mnp'
    		altcount1=altcount.strip().split(',')[ind]
    		af1=af.strip().split(',')[ind]
    	elif(bool(re.search('complex',typ)) == True):
    		ind=typ.strip().split(',').index("complex")
    		typ='complex'
    		altcount1=altcount.strip().split(',')[ind]
    		af1=af.strip().split(',')[ind]
    	elif(bool(re.search('snp',typ)) == True):
    		ind=typ.strip().split(',').index("snp")
    		typ='snp'
    		altcount1=altcount.strip().split(',')[ind]
    		af1=af.strip().split(',')[ind]
    if vrs in genes: #if the variant id is the key in genes dictionary then pull value
    	addgene = genes[vrs]
    else:
    	addgene='---'#if not add empty place holder
    if vrs in exac_info:#find variant id and pull correct info from dictionary
    	addexac=exac_info[vrs]
    	addexac=str(addexac)
    	addexac=addexac.replace("{", "").replace("}", "").replace("'", "").replace(' ', '\t') 
    	altex=addexac.split()[0]      
    	afs=addexac.split()[1]
    	ccid=addexac.split()[2]
    	if alt in altex.strip().split(','): #match variant and pull AF data
    		a=altex.strip().split(',').index(alt)
    		af2=afs.strip().split(',')[a]
    		ccid=ccid
    	else:
    		af2='NA'
    		ccid=ccid
    else:
    	af2='---'
    	ccid='---'
    #calculate the percent supporting the variant
    percent= round(100*(float(altcount1)/(float(altcount1)+float(refcount))),2)
    #format output as tab delimited file
    out= vrs+'\t'+typ+'\t'+depth+'\t'+altcount1+'\t'+str(percent)+'\t'+af1+'\t'+af2+'\t'+ccid+'\t'+addgene+'\n'
    outfile.write(out)
   



outfile.close()
