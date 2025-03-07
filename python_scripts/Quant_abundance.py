#!/usr/bin/env python3
import argparse
import sys, re, os
import subprocess
import logging
import pysam
import numpy as np
import pandas as pd
import random
from datetime import datetime
from Bio import SeqIO
from sklearn.ensemble import IsolationForest
from concurrent.futures import ProcessPoolExecutor, as_completed

####### SeqIO will ignore the blank in the fasta id !!!!!!! #############
# obtain current time
def time_current():
    current_time = f"[{str(datetime.now().replace(microsecond=0))}]"
    return current_time
def check_path(file_path_list, content_list=False):
    if content_list != False:
        # find a file with mutiple styles (e.g., in .gz or not) from mutiple dirs, and return the dir, else return false
        for each_path in list(file_path_list):
            for content in list(content_list):
                check_content = os.path.join(each_path, content)
                if os.path.exists(check_content):
                    if os.path.isdir(check_content):
                        return check_content
                    elif os.path.isfile(check_content):
                        if os.stat(check_content).st_size > 0:
                            return check_content
                        else:
                            return False
        return False
    else:
        # if file/dir not exists, return false; if file is exists but is blank, return false; if file size > 0 or dir is exists, return true
        if os.path.exists(file_path_list): # file_path_list here is the path of a folder or a file, not a list
            if os.path.isdir(file_path_list):
                return True
            elif os.path.isfile(file_path_list):
                if os.stat(file_path_list).st_size > 0:
                    return True
                else:
                    return False
        else:
            return False
def only_title(input_file):
    # check if file only have a title
    try:
        df = pd.read_table(input_file, header=0)
    except:
        return False
    if df.empty:
        return True
    else:
        return False
# a read may have multiple regions that mapped to the same or different reference contigs.
# samtools or pysam only report the alignment number rather the reads number
# if a read mapped more than one time to a reference contig, the mapped count will only be recorded as 1
# if a read mapped to multiple contigs, the mapped counts of each contig will be storaged
def count_paired_reads_per_contig(align_file_path, file_type):
    paired_reads_count = {}
    if file_type == 'bam':
        with pysam.AlignmentFile(align_file_path, "rb") as align_file:
            for read in align_file:
                ref_contig = align_file.get_reference_name(read.reference_id) # ref_contig is contig id
                if ref_contig not in paired_reads_count:
                    paired_reads_count[ref_contig] = 0
                # remove unpapped, secondary alignments, and supplementary alignments, which are useful for SNP calling but will 
                # get some false-positive counts
                if (read.is_paired) and (not read.is_unmapped) and (read.is_read1) and (not read.mate_is_unmapped) and (not read.is_secondary) and (not read.is_supplementary):
                    paired_reads_count[ref_contig] += 1
    elif file_type == 'sam':
        with pysam.AlignmentFile(align_file_path, "r") as align_file:
            for read in align_file:
                ref_contig = align_file.get_reference_name(read.reference_id) # ref_contig is contig id
                if ref_contig not in paired_reads_count:
                    paired_reads_count[ref_contig] = 0
                if (read.is_paired) and (not read.is_unmapped) and (read.is_read1) and (not read.mate_is_unmapped) and (not read.is_secondary) and (not read.is_supplementary):
                    paired_reads_count[ref_contig] += 1
    return paired_reads_count
def remove_outlier(data, calculate_method='IQR', calculate_threshold=3):
    # data is a dic, key is contig, value is length
    length_values = list(data.values())
    contigs = list(data.keys())
    # Z-score
    output = []
    if calculate_method == 'Z-score':
        std_v = np.std(length_values)
        mean_v = np.mean(length_values)
        for key in contigs:
            z_score = np.abs((float(data[key]) - mean_v) / std_v)
            if z_score > (calculate_threshold * std_v):
                contigs.remove(key)
                length_values.remove(data[key])
                del data[key]
        output = length_values
    # IQR (1.5) or Tukey's fence (1.5 or 3)
    elif calculate_method == 'IQR':
        q1 = np.percentile(length_values, 25)
        q3 = np.percentile(length_values, 75)
        iqr = q3 - q1
        lower_bound = q1 - (calculate_threshold * iqr)
        upper_bound = q3 + (calculate_threshold * iqr)
        for key in contigs:
            if (data[key] < lower_bound) or (data[key] > upper_bound):
                contigs.remove(key)
                length_values.remove(data[key])
                del data[key]
        output = length_values
    # Isolation Forest
    elif calculate_method == 'IF':
        # initialize an Isolation Forest module
        n = int(len(length_values) * 0.2)
        model = IsolationForest()
        temp_test = random.sample(length_values, n)
        df_test = np.array(temp_test).reshape(-1,1)
        model.fit(df_test)
        # calculation
        df = np.array(length_values).reshape(-1,1)
        scores = model.decision_function(df)
        output = list(df[np.where(scores > calculate_threshold)].flatten())
        for key in contigs:
            if data[key] not in output:
                contigs.remove(key)
                del data[key]
    final_length = np.mean(output)
    return [contigs, final_length]
def check_file_ending(filename):
    with open(filename, 'rb') as file:
        file.seek(-1, os.SEEK_END)
        last_byte = file.read(1)
        file.seek(-2, os.SEEK_END)
        last_two_byte = file.read(2)
        if (last_byte == b'\n') and (last_two_byte != b'\r\n'):
            return True
        elif (last_byte == b'\n') and (last_two_byte == b'\r\n'):
            return 'w'
        else:
            return False
def dict_to_df(input:dict):
    # input is a dict like: {column_1:{index_1:a, index_2:b}, column_2:{index_1:c, index_2:d}}, column is sample id and index is MAGs/contigs/gene_set/virus id
    # return df is:
    #           column_1  column2
    # index_1   a         c
    # index_2   b         d
    out_dict = {}; index = False; old = False
    columns = input.keys()
    for column in columns:
        out_dict[column] = list(input[column].values())
        index = input[column].keys()
        if old == False:
            old = index
            continue
        if index != old:
            print(f'The number of genomes is different across samples.')
            exit(1)
        else:
            old = index
    temp_df = pd.DataFrame(data=out_dict, index=index)
    return temp_df
# calculate total bases and reads using readfq
def run_readfq(input_fastq):
    try:
        cmd = f"readfq {input_fastq}"
        process = subprocess.run(cmd, stdout=subprocess.PIPE, check=True, text=True, shell=True)
        stdout = process.stdout.strip().split('\t')
        reads = re.sub(' ', '', stdout[0].split(':')[1])
        bases = re.sub(' ', '', stdout[1].split(':')[1])
        return [int(reads), int(bases)]
    except subprocess.CalledProcessError as e:
        print(f"ERROR when processing {input_fastq}: {e.stderr}")
        exit(1)
def compare_two_fastq(sample, forward, reverse):
    f_reads, f_bases = run_readfq(forward)
    r_reads, r_bases = run_readfq(reverse)
    if f_reads != r_reads:
        logger.info(f'{time_current()}     The reads number of forward ({f_reads}) and reverse ({r_reads}) is unequal, use average values.')
        reads = np.mean([f_reads, r_reads])
    else:
        reads = f_reads
    bases = f_bases + r_bases
    return [sample, reads, bases]
def process_samples(samples:dict, max_works):
    results = []
    with ProcessPoolExecutor(max_workers=max_works) as executor:
        futures = []
        for sample_name, (forward_file, reverse_file) in samples.items():
            future = executor.submit(compare_two_fastq, sample_name, forward_file, reverse_file)
            futures.append(future)
        for future in as_completed(futures):
            try:
                result = future.result()
                if result:
                    results.append(result)
            except Exception as e:
                print(f"EOOR when processing samples: {e}")
                exit(1)
    return results

##obtain options
parser = argparse.ArgumentParser(description="This scripts calculate the abundance of MAGs or multiple gene sets in samples. Pandas, Bio, minimap2, samtools, readfq were needed. Please make sure no repetitive ids were present in different MAGs/gene set.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-d', help='''Directories contain parie-end reads in fastq format. Sample_1.fastq/.fastq.gz/.fq/.fq.gz are acceptable.
                    The "Sample" reperesents the sample ID, and '_1' means the forward reads. Reverse reads should be represented by '_2'. 
                    You can specify one or more directories, script will search raw reads in all of them.''', type=str, required=True, default=False, nargs='+')
parser.add_argument('-l', '--list', help='''Sample list file. One sample per line. The title and order are based on this file.''', type=str, required=True, default=False)
parser.add_argument('-a', '--assembly', help='''Assembly results in fasta format. If you use this option, please make sure that 
                    the seq id of MAGs should all exist in the assembly, which means that you must use the unressembly MAGs. 
                    If you use the assembly module, do not give assembly files us this option, and you should put files in 
                    the folder given by "-g/--genome". PLEASE make sure that no blank present in the sequence id.''', type=str, required=False, default=False)
parser.add_argument('-g', '--genome', help='''Dir contains your MAGs/genes in fasta format. All suyffix are acceptable.
                    Make sure that the seq ids of MAGs/genes/assembly are not duplicated. If input is genes, sequences of the same kind
                    of genes should in the same fasta file, e.g., amt1.fasta, PII.fasta. In the gene module, pipeline will normalized
                    using the average length of the gene length. If input is gene clusters, you must concentrate genes in the cluster
                    to a super gene in one line. Each sequences in the fasta is a gene cluster. In the "assembly" module, 
                    each file as the assembly file of each sample. PLEASE make sure that no blank present in the sequence id.''', type=str, required=True)
parser.add_argument('--analyzed_files', help='''MAG/gene_set list file in the gene or MAG module. Pipeline will only copy files contained in 
                    the list file rather than copy all files in the dir specified by option "-g/--genome".''', type=str, required=False, default=False)
parser.add_argument('-o', help="Output dir.", type=str, required=True, default=False)
parser.add_argument('--module', help='''Choose 'gene' (each gene per file containing all sequences of this gene) module which normalize abundance
                    using the average value of sequence length in the file; 'MAG' module (each MAG per file) which normalize abundance
                    using the MAG size in bp; 'assembly' module (each sequence is an individual in all input files) which normalize abundance 
                    using the length of the sequence. Each input file is a sample in the 'assembly' module, and this module will NOT perform 
                    'between-sample' mapping (only sequences from their original sample will used as the input reference file for minimap2 mapping
                    for the sample); or choose 'virus' module (each sequence is an individual in all input files) which is very similar
                    with the 'assembly' module but it will perform 'between-sample' mapping (all input sequences will concentrate together as the 
                    reference of minimap2 mapping). ''', type=str, choices=['gene', 'MAG', 'assembly', 'virus'], required=False, default='MAG')
parser.add_argument('--method', help="Choose outlier remove method in the gene module.", type=str, choices=['IQR', 'Z-score', 'IF'], required=False, default='IQR')
parser.add_argument('--threshold', help="Threshold for outlier removing in the gene module.", type=float, required=False, default=3)
parser.add_argument('-t', '--threads', help="Threads.", type=int, required=False, default=8)
parser.add_argument('--skip_input_check', help='''Do not check input MAGs/gene_set/assembly/virus are fasta, and do not check seq id. Please make sure 
                    that input are fasta and sequence ids are renamed by the file name, like 'SY15_1_2', amoung which SY15 is the file name, 
                    '_1' is the contig id, and '_2' is the protein id.''', required=False, default=False, action='store_true')
parser.add_argument('--retain_fastqc', help="Retain temp fastqc files rather than delete them.", required=False, default=False, action='store_true')
parser.add_argument('--retain_bam', help="Retain temp bam files rather than delete them.", required=False, default=False, action='store_true')
parser.add_argument('--retain_counts', help="Retain temp paired-end mapped counts results rather than delete them.", required=False, default=False, action='store_true')
parser.add_argument('--retain_depth', help="Retain temp depth files rather than delete them.", required=False, default=False, action='store_true')
parser.add_argument('--retain_assembly', help="Retain temp concentrated assembly file rather than delete it.", required=False, default=False, action='store_true')
parser.add_argument('--update_contig_abundance', help='''Re-calculate the contig abundance using sequencing depth and mapped reads,
                    rather than using exists abundance in the contig abundance file. This option is only used within our lab to resolve 
                    problems produced by previous version. Please ignore it if you are using version v2.0.0 or higher.''', required=False, default=False, action='store_true')
parser.add_argument('--update_sample_inf', help='''Re-calculate the total base and total reads using readfq. This option is only used within our lab to resolve 
                    problems produced by previous version. Please ignore it if you are using version v2.0.0 or higher.''', required=False, default=False, action='store_true')
parser.add_argument('-v', '--version', help="Show version and exist.", action='version', version='v2.0.0')
args = parser.parse_args()
if args.update_sample_inf == True:
    args.update_contig_abundance == True
module = args.module
analyzed_files = args.analyzed_files
genomes_dir = args.genome
outdir = args.o
threads = int(args.threads)
assembly = args.assembly
raw_reads_dirs = args.d
method = args.method
threshold = args.threshold
skip_input_check = args.skip_input_check
## Set up the logger
if not check_path(outdir):
    os.mkdir(outdir)
list_file_name = os.path.basename(args.list)
log_file = os.path.join(outdir,f'Calculate_{module}_abundance_for_{list_file_name}.log')
if check_path(log_file):
    os.remove(log_file)
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)
#####################################################
########          pipeline begin          ###########
#####################################################
## get sample list
logger.info(f'{time_current()} First, carry out preparatory works.')
logger.info(f'{time_current()}   Get sample list.')
with open(args.list, 'r') as input:
    lines = input.read().splitlines()
sample_list = []
for line in lines:
    if str(line) not in sample_list:
        sample_list.append(str(line))
sample_number = len(sample_list)
logger.info(f'{time_current()}   Finished getting sample list.')

## check input raw reads files
logger.info(f'{time_current()}   Checking input raw reads...')
fastq_dict = {}
for sample in sample_list:
    forward_list = [f'{sample}_1.fastq', f'{sample}_1.fq', f'{sample}_1.fastq.gz', f'{sample}_1.fq.gz']
    reverse_list = [f'{sample}_2.fastq', f'{sample}_2.fq', f'{sample}_2.fastq.gz', f'{sample}_2.fq.gz']
    raw_reads_path_1 = check_path(raw_reads_dirs, forward_list)
    raw_reads_path_2 = check_path(raw_reads_dirs, reverse_list)
    if raw_reads_path_1 == False:
        logger.info(f'{time_current()}   Forward reads in fastq/fastq.gz/fq/fq.gz format of sample {sample} was not found in {raw_reads_dirs}. Exit...')
        exit(1)
    if raw_reads_path_2 == False:
        logger.info(f'{time_current()}   Reverse reads in fastq/fastq.gz/fq/fq.gz format of sample {sample} was not found in {raw_reads_dirs}. Exit...')
        exit(1)
    logger.info(f'{time_current()}     Found {os.path.basename(raw_reads_path_1)} and {os.path.basename(raw_reads_path_2)}.')
    # obtain total reads and sequencing bases of clean reads using readfq
    fastq_dict[sample] = [raw_reads_path_1, raw_reads_path_2]
logger.info(f'{time_current()}   Finished checking raw reads.')

# check input MAGs, make sure all files are fasta, copy MAGs which contained in the 
dir_of_analyzed_genomes = os.path.join(outdir, '00_input_files')
if not check_path(genomes_dir):
    logger.info(f'{time_current()}   Could not found {genomes_dir}, exit.')
    exit(1)
if skip_input_check == False:
    logger.info(f'{time_current()}   Get input MAGs/genes, and make sure all files are fasta.')
    if not check_path(dir_of_analyzed_genomes):
        os.mkdir(dir_of_analyzed_genomes)
    for root, dir, files in os.walk(genomes_dir):
        if root != genomes_dir:
            break
        if analyzed_files != False:
            with open(analyzed_files, 'r') as temp_in:
                analyzed_list = list(temp_in.read().splitlines()) # bin name with no suffix
        else:
            analyzed_list = list(files) # file name
        for file in files:
            bin_path = os.path.join(genomes_dir, file)
            bin_name = os.path.splitext(file)[0]
            file_suffix = os.path.splitext(file)[1]
            file_type = False
            if file in analyzed_list:
                file_type = 'f'
            elif bin_name in analyzed_list:
                file_type = 'b'
            if file_type != False:
                suffix_list = ['.fasta', '.fa', '.fna']
                if file_suffix not in suffix_list:
                    logger.info(f'{time_current()}   Found file which was not in fasta format or wa not a DNA sequence file. Please check input genome dir. Exit...')
                    exit(1)
                # open the bin file and check the final string
                end_with_n = check_file_ending(bin_path)
                if not check_path(os.path.join(dir_of_analyzed_genomes, file)):
                    if end_with_n == True:
                        os.system(f"cp {bin_path} {dir_of_analyzed_genomes}")
                    elif end_with_n == False:
                        outbin_path = open(os.path.join(dir_of_analyzed_genomes, file), 'w')
                        with open(bin_path, 'r') as temp_in:
                            content = temp_in.read()
                            outbin_path.writelines(f"{content}\n")
                        outbin_path.close()
                    else:
                        logger.info(f'{time_current()}   Found {bin_path} in windows type, exit.')
                        exit(1)
                if file_type == 'f':
                    analyzed_list.remove(file)
                elif file_type == 'b':
                    analyzed_list.remove(bin_name)
        if len(analyzed_list) > 0:
            logger.info(f'{time_current()}   Query files "{analyzed_list}" were not found in the input dir, please check input list file or input dir. Exit...')
            exit(1)
elif skip_input_check == True:
    logger.info(f'{time_current()}   Skip input check, please make sure no other files exists in the input dir.')
    dir_of_analyzed_genomes = genomes_dir
## Get concentrated assembly from MAGs or just using the specified assembly file, for gene module and MAG module only.
if (assembly != False) and (module == 'gene' or module == 'assembly'):
    logger.info(f'{time_current()}   ERROR! You could not give assembly file using "-a" option in the gene or assembly module, exit.')
    exit(1)
elif (assembly != False) and (module == 'MAG' or module == 'virus'):
    assembly_file = assembly
    logger.info(f'{time_current()}   Concentrated assembly file was given, pipeline will use it.')
elif (assembly == False) and (module == 'MAG' or module == 'gene' or module == 'virus'):
    logger.info(f'{time_current()}   Get concentrated assembly from MAGs or just using the specified assembly file.')
    assembly_file = os.path.join(outdir, f'concentrated_file_of_{module}.fa')
    if not check_path(assembly_file):
        exit_value = os.system(f'cat {dir_of_analyzed_genomes}/* > {assembly_file}')
        if exit_value == 0:
            logger.info(f'{time_current()}   Finished concentrating assembly from MAGs.')
        else:
            logger.info(f'{time_current()}   Could not concentrate assembly using your MAGs.')
            exit(1)
    else:
        logger.info(f'{time_current()}   Found previous concentrated assembly file.')

## Read MAGs to obtain all sequence ids, length, and host MAG. 
contigid_2_genomelen = {}
contigid_2_genomeid = {}
contig_2_lenth = {}
genome_list = []
genome_2_contigs = {}
contig_set = set()
logger.info(f'{time_current()}   Get contig id list and genome size/gene_length/contig_length for query MAGs/gene_set/assembly.')
for root,dirs,files in os.walk(dir_of_analyzed_genomes):
    if root != dir_of_analyzed_genomes:
        break
    for file in files: # files represent query MAGs or genes, each MAG/gene in each file
        input = os.path.join(dir_of_analyzed_genomes, file)
        genome = os.path.splitext(file)[0] # MAG_id/gene_set_id/assembly_id
        genome_list.append(genome) # get the MAG/gene_set/assembly list, will also get file list in virus module but useless
        if (module == 'MAG') and (genome not in contigid_2_genomelen):
            contigid_2_genomelen[genome] = 0
        if genome not in genome_2_contigs:
            genome_2_contigs[genome] = []
        if (module == 'gene' or module == 'assembly') and (genome not in contigid_2_genomelen):
            contigid_2_genomelen[genome] = {} # in gene/assembly module, key is gene_name/sample_id, value is a dic whose key is contig id and value is contig length
        for records in SeqIO.parse(input, 'fasta'):
            contigid_2_genomeid[records.id] = genome
            if module == 'MAG':
                contig_set.add(records.id) # in MAG module, just storage all contig ids
                genome_2_contigs[genome].append(records.id) # obtain contigs list for each MAG
                contig_2_lenth[records.id] = int(len(records.seq))
                contigid_2_genomelen[genome] += contig_2_lenth[records.id] 
            elif module == 'gene' or module == 'assembly':
                # for gene module, 'genome' is gene set id (the input file name without suffix), for assembly module, 'genome' is sample/assembly id
                contig_2_lenth[records.id] = contigid_2_genomelen[genome][records.id] = int(len(records.seq))
                genome_2_contigs[genome].append(records.id)
            elif module == 'virus':
                # in the virus module, just storage the length of each sequences
                contig_set.add(records.id)
                contig_2_lenth[records.id] = int(len(records.seq))
logger.info(f'{time_current()}   Finished getting contig id list and genome size for each MAG.')
# remove outliers in the gene module
if module == 'gene':
    logger.info(f'{time_current()}   Remove outliers in the gene set.')
    for genome in genome_list:
        if genome not in genome_2_contigs:
            genome_2_contigs[genome] = []
        temp_list = []
        # remove_outlier return a list, list[0] is a list containing contig ids with out outliers, and list[1] is mean length of the gene
        temp_list = remove_outlier(contigid_2_genomelen[genome], calculate_method=method, calculate_threshold=threshold)
        contigid_2_genomelen[genome] = temp_list[1] # value of old dic is a dic, now it is the average length of the gene
        contig_set.update(temp_list[0]) # obatin all contig ids without outliers in the 'gene' module 
        genome_2_contigs[genome] = temp_list[0]
    # remove outlier sequences from concentrated assembly
    os.remove(assembly_file)
    with open(assembly_file, 'w') as tempout:
        for records in SeqIO.parse(assembly_file, 'fasta'):
            if records.id in contig_set:
                tempout.writelines(f'>{records.id}\n{records.seq}\n')
    logger.info(f'{time_current()}   Finished removing outliers.')

# if assembly was given, we need to know the sequences which were not belong to any MAGs.
# Pipeline will calculate the TPM/RPKM/DPGB/PPB for unbinned sequences. 
# all these unbinned sequences will be considered as a 'Unbinned_organisms'
if (assembly != False) and (module == 'MAG' or module == 'virus'):
    logger.info(f'{time_current()}   Getting extra sequences ids from given assembly.')
    temp_contig_dic = {}
    for records in SeqIO.parse(assembly_file, 'fasta'):
        temp_contig_dic[records.id] = int(len(records.seq))
    temp_keys_set = set(temp_contig_dic.keys())
    if temp_keys_set.issuperset(contig_set):
        different_ids_set = temp_keys_set.difference(contig_set)
    else:
        logger.info(f'{time_current()}   Some contigs "{contig_set.difference(temp_keys_set)}" in MAGs was not present in assembly, exit.')
        exit(1)
    if module == 'MAG':
        if 'Unbinned_organisms' not in genome_list:
            genome_list.append('Unbinned_organisms')
        contigid_2_genomelen['Unbinned_organisms'] = 0
        genome_2_contigs['Unbinned_organisms'] = list(different_ids_set)
        for each_diff in different_ids_set:
            contigid_2_genomelen['Unbinned_organisms'] += temp_contig_dic[each_diff]
            contigid_2_genomeid[each_diff] = 'Unbinned_organisms'
    contig_2_lenth = temp_contig_dic
    contig_set = temp_keys_set # transfer binned contig list to total contig list
    del temp_contig_dic
    del temp_keys_set
    del different_ids_set
    logger.info(f'{time_current()}   Finished getting extra sequences ids from given assembly.')

# make need dirs 
align_dir = os.path.join(outdir, '01_align_files')
if not check_path(align_dir): os.mkdir(align_dir)
depth_dir = os.path.join(outdir, '02_depth_files')
if not check_path(depth_dir): os.mkdir(depth_dir)
mapped_count_dir = os.path.join(outdir, '03_mapped_count_results')
if not check_path(mapped_count_dir): os.mkdir(mapped_count_dir)
sample_inf_dir = os.path.join(outdir, '04_sample_information')
if not check_path(sample_inf_dir): os.mkdir(sample_inf_dir)
abundance_dir = os.path.join(outdir, '05_contig_abundance')
if not check_path(abundance_dir): os.mkdir(abundance_dir)
mapped_reads_of_sample = {} # storage the total mapped reads counts (value) for each sample (key)
total_depth_of_sample = {} # key is sample id, value is total depth of assembly of this sample
final_results = {}
final_results['contig_mapped_counts'] = {} # key is sample id, value is a dic whose key is contig id, value is mapped count
final_results['contig_depth'] = {} # key is sample id, value is a dic whose key is contig id, value is total depth of contig
final_results['contig_length'] = {}
final_results['contig_DPGB_in_samples'] = {} # key is sample id, value is a dic whose key is contig id, value is DPGB
final_results['contig_PPB_in_samples'] = {} # key is sample id, value is a dic whose key is contig id, value is PPB
final_results['contig_RPKM_in_samples'] = {} # key is sample id, value is a dic whose key is contig id, value is RPKM
final_results['contig_TPM_in_samples'] = {} # key is sample id, value is a dic whose key is contig id, value is TPM

## run align step to obtain mapped reads counts
logger.info(f'{time_current()} Second, performing abundance calculation steps.')
logger.info(f'{time_current()}   Mapping parir-end reads to concentrated assembly, and calculating abundance of contigs...')
if (module == 'MAG') or (module == 'gene') or (module == 'virus'):
    if not check_path(f"{assembly_file}.fai"):
        logger.info(f'{time_current()}     This is {module} module, now perform faidx for the concentrated assembly.')
        # first, faidx the concentrated assembly file for gene and MAG module
        exit_value = os.system(f'samtools faidx {assembly_file}')
        if exit_value == 0:
            logger.info(f'{time_current()}     Finished faidx for the concentrated assembly.')
        else:
            logger.info(f'{time_current()}     Something wrong for faidx the concentrated assembly, exit.'); exit(1)
    else:
        logger.info(f'{time_current()}     The faidx for the concentrated assembly was already done.')

total_reads_of_sample = {}
total_bases_of_sample = {}
sample_inf_file = os.path.join(sample_inf_dir, 'sample_inf.txt')
temp_sample_inf_file = os.path.join(sample_inf_dir, 'sample_inf-temp.txt')
logger.info(f'{time_current()}   calculating contig abundance...')
if (not check_path(sample_inf_file)) or (only_title(temp_sample_inf_file)) or (args.update_sample_inf == True):
    logger.info(f'{time_current()}   Updating total reads and total bases...')
    sample_inf = open(sample_inf_file, 'w')
    sample_inf.writelines(f"Sample_id\tTotal_mapped_reads (pair-end)\tTotal_mapped_depth\tTotal_reads\tTotal_bases\n")
    sample_inf.close()
    if check_path(temp_sample_inf_file) and (not only_title(temp_sample_inf_file)): # empty file or a file only containing a title
        logger.info(f'{time_current()}     Reading total reads and total bases from previous result...')
        temp_df = pd.read_table(temp_sample_inf_file, header=0)
        for index, row in temp_df.iterrows():
            sample_id = row['Sample_id']
            reads = row['Total_reads']
            bases = row['Total_bases']
            total_reads_of_sample[sample_id] = reads
            total_bases_of_sample[sample_id] = bases
    else:
        logger.info(f'{time_current()}     Calculating total reads and total bases...')
        information = process_samples(fastq_dict, threads)
        temp_sample_inf = open(temp_sample_inf_file, 'w')
        temp_sample_inf.writelines(f"Sample_id\tTotal_reads\tTotal_bases\n")
        for each_list in information:
            sample_id, reads, bases = each_list
            total_reads_of_sample[sample_id] = reads
            total_bases_of_sample[sample_id] = bases
            temp_sample_inf.writelines(f"{sample_id}\t{reads}\t{bases}\n")
        temp_sample_inf.close()
        logger.info(f'{time_current()}   Total reads and total bases were writen to {temp_sample_inf_file}.')
else:
    logger.info(f'{time_current()}   Reading total reads and total bases from previous result...')
    if check_path(sample_inf_file) and (not only_title(sample_inf_file)):
        temp_df = pd.read_table(sample_inf_file, header=0)
    elif check_path(temp_sample_inf_file) and (not only_title(temp_sample_inf_file)):
        temp_df = pd.read_table(temp_sample_inf_file, header=0)
    else:
        logger.info(f'{time_current()}   Could not found samples inf file, please add "--update_sample_inf" option.')
        exit(1)
    if 'Total_bases' not in temp_df:
        logger.info(f'{time_current()}   Did not find total bases column in sample inf file, please add "--update_sample_inf" option.')
        exit(1)
    for index, row in temp_df.iterrows():
        sample_id = row['Sample_id']
        reads = row['Total_reads']
        bases = row['Total_bases']
        total_reads_of_sample[sample_id] = reads
        total_bases_of_sample[sample_id] = bases

for sample in sample_list:
    reads_f, reads_r = fastq_dict.get(sample)
    align_file = os.path.join(align_dir, sample)
    depth_file = os.path.join(depth_dir, f'{sample}.depth.txt')
    abundance_file = os.path.join(abundance_dir, f'Contig_abundance_inf_of_{sample}.txt')
    mapped_count_file = os.path.join(mapped_count_dir, f'{sample}.mapped_counts.txt')
    # total mapped count of the sample
    sample_reads = float(total_reads_of_sample[sample])
    sample_bases = float(total_bases_of_sample[sample])
    if module == 'assembly': contig_set = genome_2_contigs[sample]
    # mapping
    if (check_path(abundance_file)) and (not only_title(abundance_file)):
        logger.info(f'{time_current()}     Seems that the abundance calculation step for sample {sample} was already done.')
        if args.update_contig_abundance == True:
            logger.info(f'{time_current()}     Updating exists abundance results of sample {sample}.')
            inf_df = pd.read_table(abundance_file, header=0)
            if 'Total_depth' in inf_df.columns:
                inf_df = inf_df.rename(columns={'Total_depth':'Total_mapped_depth'})
            if 'DPGM' in inf_df.columns:
                inf_df = inf_df.rename(columns={'DPGM':'DPGB'})
            if 'PPM' in inf_df.columns:
                inf_df = inf_df.rename(columns={'PPM':'PPB'})
            if 'dTPM' in inf_df.columns:
                inf_df = inf_df.rename(columns={'dTPM':'PPB'})
            #(f"Name\tsample_id\tMAG_id\tlength(bp)\tMapped_reads\tTotal_mapped_depth\tTPM\tRPKM\tPPB\tDPGB\n")
            temp = os.path.join(abundance_dir, f'Contig_abundance_inf_of_{sample}_old.txt')
            os.system(f"mv {abundance_file} {temp}")
            inf_df['DPGB'] = inf_df['Total_mapped_depth'] / (inf_df['length(bp)'] * sample_bases) * 1e+12
            inf_df['RPKM'] = inf_df['Mapped_reads'] / (inf_df['length(bp)'] * sample_reads) * 1e+9
            total_DPGB = inf_df['DPGB'].sum()
            total_RPKM = inf_df['RPKM'].sum()
            sample_total_mapped_reads = inf_df['Mapped_reads'].sum()
            sample_total_depth = inf_df['Total_mapped_depth'].sum()
            if total_DPGB > 0:
                inf_df['PPB'] = (inf_df['DPGB'] / total_DPGB) * 1e+9
            else:
                inf_df['PPB'] = 0
            if total_RPKM > 0:
                inf_df['TPM'] = (inf_df['RPKM'] / total_RPKM) * 1e+6
            else:
                inf_df['TPM'] = 0
            inf_df.to_csv(abundance_file, header=True, index=False, mode='w', sep='\t', na_rep=0)
        if args.update_sample_inf == True:
            sample_inf = open(sample_inf_file, 'a')
            sample_inf.writelines(f"{sample}\t{sample_total_mapped_reads}\t{sample_total_depth}\t{sample_reads}\t{sample_bases}\n")
            sample_inf.close()
        continue
    if check_path(f'{align_file}.sort.bam') or check_path(depth_file) or check_path(mapped_count_file):
        logger.info(f'{time_current()}     Seems that the mapping step for sample {sample} was already done.')
    else:
        # align reads to concentrated assembly
        if (module == 'MAG') or (module == 'gene') or (module == 'virus'):
            logger.info(f'{time_current()}     Mapping reads to concentrated assembly of sample {sample} using minimap2.')
        elif module == 'assembly':
            # in the assembly module, only perform the within sample mapping
            logger.info(f'{time_current()}     Mapping raw reads to the assembly or assembly genes file of sample {sample} using minimap2.')
            assembly_file = os.path.join(dir_of_analyzed_genomes, f'{sample}{file_suffix}')
            # fadxi for each assembly in the assembly module
            logger.info(f'{time_current()}     Perform faidx for the assembly (or DNA sequences of assembly gene repdiction file) of sample {sample}.')
            exit_value = os.system(f'samtools faidx {assembly_file}')
            if exit_value == 0:
                logger.info(f'{time_current()}     Finished faidx for the concentrate assembly.')
            else:
                logger.info(f'{time_current()}     Something wrong for faidx step of sample {sample}.'); exit(1)
        # mapping using minimap2
        exit_value = os.system(f'minimap2 -ax sr -t {threads} -o {align_file}.sam {assembly_file} {reads_f} {reads_r}')
        if exit_value == 0:
            logger.info(f'{time_current()}     Finished mapping raw reads for sample {sample}.')
        else:
            logger.info(f'{time_current()}     Something wrong for mapping raw reads for sample {sample}.'); exit(1)
        # transfering sam to bam file
        logger.info(f'{time_current()}     Tranfering sam to bam for sample {sample}.')
        exit_value = os.system(f'samtools view -@ {threads} -bt {assembly_file}.fai -o {align_file}.bam {align_file}.sam && rm -f {align_file}.sam')
        if exit_value == 0:
            # remove the faidx file in the assembly module
            if module == 'assembly': os.system(f"rm -f {assembly_file}.fai")
            logger.info(f'{time_current()}     Finished tranfering sam to unsorted bam for sample {sample}.')
        else:
            logger.info(f'{time_current()}     Something wrong for tranfering sam to unsorted bam for sample {sample}.'); exit(1)
        # sortting bam file
        logger.info(f'{time_current()}     Tranfering sam to sorted bam for sample {sample}.')
        exit_value = os.system(f'samtools sort -@ {threads} -O BAM -o {align_file}.sort.bam {align_file}.bam && rm -f {align_file}.bam')
        if exit_value == 0:
            logger.info(f'{time_current()}     Finished sortting bam file for sample {sample}.')
        else:
            logger.info(f'{time_current()}     Something wrong for sortting bam file for sample {sample}.'); exit(1)
        logger.info(f'{time_current()}     Finished mapping steps for sample {sample}.')
    
    # calculate the mapped counts and average depth for each contig
    if check_path(depth_file) or check_path(mapped_count_file):
        logger.info(f'{time_current()}     Seems that the "samtools depth" calculation for sample {sample} were already done.')
    else:
        logger.info(f'{time_current()}     Calculating depth of query MAGs/gene_set in sample {sample}.')
        exit_value = os.system(f'samtools depth -aa -@ {threads} -o {depth_file} {align_file}.sort.bam')
        if exit_value == 0:
            logger.info(f'{time_current()}     Finished "samtools depth" calculation for sample {sample}.')
        else:
            logger.info(f'{time_current()}     Something wrong for "samtools depth" calculation for sample {sample}.'); exit(1)

    # calculate mapped reads counts for each contig
    logger.info(f'{time_current()}     Calculating mapped reads counts for each query contig in sample {sample}.')
    # write mapped reads count to file rather than dictionary.
    if check_path(mapped_count_file):
        logger.info(f'{time_current()}     Seems that the mapped counts for sample {sample} were already wrote to file {mapped_count_file}.')
    else:
        dic_of_mapped_counts = count_paired_reads_per_contig(f'{align_file}.sort.bam', 'bam')
        with open(f'{mapped_count_file}', 'w') as temp_mapped_counts:
            for contig in contig_set:
                # for thses module, we need abundan information of all contigs in all samples
                temp_mapped_counts.writelines(f'{contig}\t')
                if contig in dic_of_mapped_counts:
                    temp_mapped_counts.writelines(f'{dic_of_mapped_counts[contig]}\n')
                else:
                    temp_mapped_counts.writelines('0\n')
        if args.retain_bam == False: os.remove(f'{align_file}.sort.bam')
        logger.info(f'{time_current()}     Finished calculating mapped reads counts for each query contig in sample {sample}.')

    ### begin calculation ###
    logger.info(f'{time_current()}     Calculating contig abundance (TPM/RPKM/DPGB/PPB) for sample {sample}...')
    output = open(abundance_file, 'w')
    sample_inf = open(sample_inf_file, 'a')
    if not check_path(sample_inf_file):
        sample_inf.writelines(f"Sample_id\tTotal_mapped_reads (pair-end)\tTotal_mapped_depth\tTotal_reads\tTotal_bases\n")
    output.writelines(f"Name\tsample_id\tMAG_id\tlength(bp)\tMapped_reads\tTotal_mapped_depth\tTPM\tRPKM\tPPB\tDPGB\n")
    final_results['contig_mapped_counts'][sample] = {} # key is sample id, value is a dic whose key is contig id, value is mapped count
    final_results['contig_depth'][sample] = {} # key is contig id, value is total depth of contig
    final_results['contig_length'][sample] = {}
    final_results['contig_DPGB_in_samples'][sample] = {} # key is contig id, value is DPGB of contig
    final_results['contig_PPB_in_samples'][sample] = {} # key is contig id, value is PPB of contig
    final_results['contig_RPKM_in_samples'][sample] = {} # key is contig id, value is RPKM of contig
    final_results['contig_TPM_in_samples'][sample] = {} # key is contig id, value is TPM of contig
    '''
    In the new version, calculate (1) the contig/MAG length in samtools results, (2) depth, both using pandas directly rather
    than use for loop. In the old version, I use for loop to traverse the depth file. In the new version, I traverse genome list,
    obtain the contig list of genomes, then just calculate these using pandas.
    There are three colomns in the depth file of samtools, [0] is contig id, [1] is position, [2] is depth.
    '''
    total_RPKM = 0 # Initialize the total RPKM for each sample 
    total_DPGB = 0 # Initialize the total DPGB for each sample
    mapped_reads_of_sample[sample] = 0 # key is sample id, value is total mapped counts of this sample
    total_depth_of_sample[sample] = 0 # key is sample id, value is total depth of assembly of this sample

    ############### get mapped counts information #########
    mapped_count_file = os.path.join(mapped_count_dir, f'{sample}.mapped_counts.txt')
    logger.info(f'{time_current()}     Reading mapped count result {mapped_count_file}...')
    # read the mapped count file, get mapped count for each contig, make a dictionary
    df_counts = pd.read_table(mapped_count_file, header=None)
    final_results['contig_mapped_counts'][sample] = df_counts.set_index(0)[1].to_dict()
    # get mapped count for sample
    mapped_reads_of_sample[sample] = np.sum(df_counts.iloc[:,1])
    ################# get depth information ##################
    logger.info(f'{time_current()}     Reading depth result {depth_file}...')
    depth_file = os.path.join(depth_dir, f'{sample}.depth.txt')
    chunksize = 1e+6
    sample_total_depth = 0
    depth_dict = {}; length_dict = {}
    for chunk in pd.read_table(depth_file, header=None, chunksize=chunksize):
        sample_total_depth += np.sum(chunk.iloc[:,2])
        grouped_df = chunk.groupby(by=0, as_index=True, sort=False)
        dict_b = grouped_df[2].sum().to_dict()
        for key_b, value_b in dict_b.items():
            if key_b in final_results['contig_depth'][sample]:
                final_results['contig_depth'][sample][key_b] += value_b
            else:
                final_results['contig_depth'][sample][key_b] = value_b
        dict_a = grouped_df[1].max().to_dict()
        for key_a, value_a in dict_a.items():
            if key_a in final_results['contig_length'][sample]:
                if value_a > final_results['contig_length'][sample][key_a]:
                    final_results['contig_length'][sample][key_a] = value_a
            else:
                final_results['contig_length'][sample][key_a] = value_a
    total_depth_of_sample[sample] = sample_total_depth

    # write sample information, including sample ID, total mapped reads, total depth, and total reads
    sample_inf.writelines(f"{sample}\t{mapped_reads_of_sample[sample]}\t{total_depth_of_sample[sample]}\t{int(sample_reads)}\t{int(sample_bases)}\n")

    # calculate abundance for contigs, in the MAG or gene module, all contigs were used to perform compare between samples
    logger.info(f'{time_current()}     Calculating the abundance of each sequence in the sample {sample}...')
    for contig in contig_set:
        # in gene/MAG module, contig_set is a set contain all contig ids; in 'assembly' module, it only contains contig ids in this sample
        length_of_contig = final_results['contig_length'][sample][contig]
        if contig_2_lenth[contig] != length_of_contig:
            logger.info(f'{time_current()}     PLEASE NOTE, length of contig {contig} is {length_of_contig} in samtools results but {contig_2_lenth[contig]} in fasta file. Pipeline will use the value in fasta file.')
            length_of_contig = float(contig_2_lenth[contig])
        final_results['contig_DPGB_in_samples'][sample][contig] = (float(final_results['contig_depth'][sample][contig]) / (length_of_contig * float(sample_bases))) * 1e+12
        final_results['contig_RPKM_in_samples'][sample][contig] = (float(final_results['contig_mapped_counts'][sample][contig]) / (length_of_contig * float(sample_reads))) * 1e+9
    total_DPGB = np.sum(list(final_results['contig_DPGB_in_samples'][sample].values()))
    total_RPKM = np.sum(list(final_results['contig_RPKM_in_samples'][sample].values()))
    if total_DPGB == 0:
        logger.info(f'{time_current()}     PLEASE NOTE! The total DPGB was zero in sample {sample}, which means no reads of this sample were mapped to any MAG/gene_set.')
    if total_RPKM == 0:
        logger.info(f'{time_current()}     PLEASE NOTE! The total RPKM was zero in sample {sample}, which means no reads of this sample were mapped to any MAG/gene_set.')
    for contig in contig_set:
        if total_DPGB > 0:
            final_results['contig_PPB_in_samples'][sample][contig] = (final_results['contig_DPGB_in_samples'][sample][contig] / total_DPGB) * 1e+9
        else:
            final_results['contig_PPB_in_samples'][sample][contig] = 0
        if total_RPKM > 0:
            final_results['contig_TPM_in_samples'][sample][contig] = (final_results['contig_RPKM_in_samples'][sample][contig] / total_RPKM) * 1e+6
        else:
            final_results['contig_TPM_in_samples'][sample][contig] = 0
    logger.info(f'{time_current()}     Writting abundance results of sample {sample}...')
    for contig in final_results['contig_mapped_counts'][sample]:
        contig_len = contig_2_lenth[contig]
        mag = contigid_2_genomeid[contig]
        depth = final_results['contig_depth'][sample][contig]
        reads = final_results['contig_mapped_counts'][sample][contig]
        tpm = final_results['contig_TPM_in_samples'][sample][contig]
        rpkm = final_results['contig_RPKM_in_samples'][sample][contig]
        PPB = final_results['contig_PPB_in_samples'][sample][contig]
        dpgm = final_results['contig_DPGB_in_samples'][sample][contig]
        output.writelines(f"{contig}\t{sample}\t{mag}\t{contig_len}\t{reads}\t{depth}\t{tpm}\t{rpkm}\t{PPB}\t{dpgm}\n")
    output.close()
    sample_inf.close()
    del df_counts
    if args.retain_depth == False: os.remove(depth_file)
    if args.retain_counts == False: os.remove(mapped_count_file)
    logger.info(f'{time_current()}     Finished calcaulting contig abundance for sample {sample}.')
# remove the faidx file for the concentrated assembly in the MAG and gene module
logger.info(f'{time_current()} Finished mapping and contig abundance calculation steps.')
del final_results
'''
In the assembly modle, we just need contig abundance in their own sample, in the MAG or gene module, we need calculate the abundance 
of a MAG or a gene_set. In the virus module, we just need combine contig abundance files of all samples.
'''
if module == 'virus':
    file_dict = {'DPGB':f'{outdir}/Abundance-depth_DPGB.txt', 'PPB':f'{outdir}/Percentage-depth_PPB.txt', 
                'RPKM':f'{outdir}/Abundance-counts_RPKM.txt', 'TPM':f'{outdir}/Percentage-counts_TPM.txt', 
                'Mapped_reads':f'{outdir}/Information-count.txt', 'Total_mapped_depth':f'{outdir}/Information-depth.txt'}
    # read all abundance df
    temp_dict = {}
    for sample in sample_list:
        abundance_file = os.path.join(abundance_dir, f'Contig_abundance_inf_of_{sample}.txt')
        logger.info(f'{time_current()}     Reading contig abundance from {abundance_file}...')
        inf_df = pd.read_table(abundance_file, header=0)
        if 'Total_depth' in inf_df.columns:
            inf_df = inf_df.rename(columns={'Total_depth':'Total_mapped_depth'})
        if 'DPGM' in inf_df.columns:
            inf_df = inf_df.rename(columns={'DPGM':'DPGB'})
        if 'PPM' in inf_df.columns:
            inf_df = inf_df.rename(columns={'PPM':'PPB'})
        if 'dTPM' in inf_df.columns:
            inf_df = inf_df.rename(columns={'dTPM':'PPB'})
        logger.info(f'{time_current()}     Merging abundance values for each methods...')
        for each_method in file_dict.keys():
            if each_method not in temp_dict:
                temp_dict[each_method] = pd.DataFrame()
            temp_df = inf_df[['Name', each_method]]
            temp_df.rename(columns={each_method:sample}, inplace=True)
            if len(temp_dict[each_method]) == 0:
                temp_dict[each_method] = temp_df
            else:
                temp_dict[each_method] = pd.merge(temp_dict[each_method], temp_df, on='Name', how='left')
    for each_method in file_dict.keys():
        temp_dict[each_method].to_csv(file_dict[each_method], sep='\t', header=True, index=False, na_rep=0, mode='w')
        logger.info(f'{time_current()}   Write {each_method} results to file: {file_dict[each_method]}.')
## calculate abundance for MAG/gene_set
elif (module == 'MAG') or (module == 'gene'):
    DPGB_file = f'{outdir}/Abundance-depth_DPGB.txt'
    PPB_file = f'{outdir}/Percentage-depth_PPB.txt'
    RPKM_file = f'{outdir}/Abundance-counts_RPKM.txt'
    TPM_file = f'{outdir}/Percentage-counts_TPM.txt'
    depth_inf_file = f'{outdir}/Information-depth.txt'
    count_inf_file = f'{outdir}/Information-count.txt'
    final_results = {}
    logger.info(f'{time_current()} Third, calculating abundance for MAG/gene_set/virus.')
    logger.info(f'{time_current()}   Reading sample information from {sample_inf_file}...')
    # clean old dics
    mapped_reads_of_sample = {}
    mapped_bases_of_sample = {}
    total_reads_of_sample = {}
    total_bases_of_sample = {}
    df = pd.read_table(sample_inf_file, header=0)
    if 'Total_depth' in df.columns:
        df = df.rename(columns={'Total_depth':'Total_mapped_depth'})
    for index, rows in df.iterrows():
        # header: Sample_id\tTotal_mapped_reads (pair-end)\tTotal_mapped_depth\tTotal_reads\tTotal_bases\n
        mapped_reads_of_sample[str(rows['Sample_id'])] = rows['Total_mapped_reads (pair-end)']
        mapped_bases_of_sample[str(rows['Sample_id'])] = rows['Total_mapped_depth']
        total_bases_of_sample[str(rows['Sample_id'])] = rows['Total_bases']
        total_reads_of_sample[str(rows['Sample_id'])] = rows['Total_reads']
    final_results['genome_size'] = {} # key is genome id, value is the genome size
    final_results['MAG_total_counts'] = {} # key is genome id, value is the total mapped reads number of the MAG
    final_results['MAG_total_depth'] = {} # key is genome id, value is the total depth of all posiitions of the MAG
    final_results['MAG_RPKM_in_samples'] = {} # key is genome id, value is a dictionary whose key is sample id and value is the RPKM value of the genome
    final_results['MAG_TPM_in_samples'] = {} # key is genome id, value is a dictionary whose key is sample id and value is the TPM value of the genome
    final_results['MAG_DPGB_in_samples'] = {} # key is genome id, value is a dictionary whose key is sample id and value is the DPGB value of the genome
    final_results['MAG_PPB_in_samples'] = {} # key is genome id, value is a dictionary whose key is sample id and value is the PPB value of the genome
    logger.info(f'{time_current()}   Calculating the MAG abundance (RPKM/TPM/DPGB/PPB) between samples...')
    for sample in sample_list:
        abundance_file = os.path.join(abundance_dir, f'Contig_abundance_inf_of_{sample}.txt')
        sample_reads = float(total_reads_of_sample[sample])
        sample_bases = float(total_bases_of_sample[sample])
        total_RPKM = 0 # Initialize the total RPKM for each sample 
        total_DPGB = 0 # Initialize the total DPGB for each sample
        final_results['genome_size'][sample] = {} # key is genome id, value is the genome size
        final_results['MAG_total_counts'][sample] = {} # key is genome id, value is the total mapped reads number of the MAG
        final_results['MAG_total_depth'][sample] = {} # key is genome id, value is the total depth of all posiitions of the MAG
        final_results['MAG_RPKM_in_samples'][sample] = {} # key is genome id, value is a dictionary whose key is sample id and value is the RPKM value of the genome
        final_results['MAG_TPM_in_samples'][sample] = {} # key is genome id, value is a dictionary whose key is sample id and value is the TPM value of the genome
        final_results['MAG_DPGB_in_samples'][sample] = {} # key is genome id, value is a dictionary whose key is sample id and value is the DPGB value of the genome
        final_results['MAG_PPB_in_samples'][sample] = {} # key is genome id, value is a dictionary whose key is sample id and value is the PPB value of the genome
        for genome_id in genome_list:
            final_results['MAG_total_depth'][sample][genome_id] = 0
            final_results['MAG_total_counts'][sample][genome_id] = 0 # initialization
            if module == 'MAG':
                final_results['genome_size'][sample][genome_id] = 0 # initialization
            elif module == 'gene':
                final_results['genome_size'][sample][genome_id] = [] # initialization
            final_results['genome_size'][sample][genome_id] = 0 # genome size of the MAG in the result of samtools, we will compare this value with the value obtained from the fasta file of MAG.
            final_results['MAG_total_depth'][sample][genome_id] = 0 # key is sample ID, value is the total depth of all positions of the MAG in the sample
            final_results['MAG_total_counts'][sample][genome_id] = 0 # key is sample ID, value is the total mapped reads number of the MAG in the sample
            final_results['MAG_DPGB_in_samples'][sample][genome_id] = 0 # key is sample ID, value is the DPGB of each MAG in the sample
            final_results['MAG_PPB_in_samples'][sample][genome_id] = 0 # key is sample ID, value is the PPB of each MAG in the sample
            final_results['MAG_RPKM_in_samples'][sample][genome_id] = 0 # key is sample ID, value is the RPKM of each MAG in the sample
            final_results['MAG_TPM_in_samples'][sample][genome_id] = 0 # key is sample ID, value is the TPM of each MAG in the sample
        # get information for each genome
        '''
        in the genome/gene_set size calculation, 'MAG' module need total length of contigs of the MAG,
        while 'gene' module need the average length of the gene. The last position is the length of the contig in
        the samtools depth result.
        '''
        # read abundance information result
        logger.info(f'{time_current()}     Reading contig abundance from {abundance_file}...')
        inf_df = pd.read_table(abundance_file, header=0)
        if 'Total_depth' in inf_df.columns:
            inf_df = df.rename(columns={'Total_depth':'Total_mapped_depth'})
        if 'DPGM' in inf_df.columns:
            inf_df = inf_df.rename(columns={'DPGM':'DPGB'})
        if 'PPM' in inf_df.columns:
            inf_df = inf_df.rename(columns={'PPM':'PPB'})
        if 'dTPM' in inf_df.columns:
            inf_df = inf_df.rename(columns={'dTPM':'PPB'})
        logger.info(f'{time_current()}     Calculatting the MAGs abundance in sample {sample}...')
        # group by genome/gene_set id
        grouped_df = inf_df.groupby(by='MAG_id', as_index=True, sort=True)
        final_results['MAG_total_depth'][sample] = grouped_df['Total_mapped_depth'].sum().to_dict()
        final_results['MAG_total_counts'][sample] = grouped_df['Mapped_reads'].sum().to_dict()
        if module == 'MAG':
            final_results['genome_size'][sample] = grouped_df['length(bp)'].sum().to_dict()
        elif module == 'gene':
            final_results['genome_size'][sample] = grouped_df['length(bp)'].mean().to_dict()
        
        ## begin abundance calculation ##
        for genome_id in genome_list:
            # RPKM, 'contigid_2_genomelen[genome_id] contains the average length of gene or the genome size of MAG'
            genome_size_of_MAG = float(final_results['genome_size'][sample][genome_id])
            if genome_size_of_MAG != contigid_2_genomelen[genome_id]:
                logger.info(f'{time_current()}     PLEASE NOTE, the size of genome {genome_id} is {genome_size_of_MAG} in contig abundance results but {contigid_2_genomelen[genome_id]} in fasta file. Pipeline will use the value in fasta file.')
                genome_size_of_MAG = float(contigid_2_genomelen[genome_id])
            final_results['MAG_RPKM_in_samples'][sample][genome_id] = (float(final_results['MAG_total_counts'][sample][genome_id]) / (genome_size_of_MAG * sample_reads)) * 1e+9
            # DPGB
            final_results['MAG_DPGB_in_samples'][sample][genome_id] = (float(final_results['MAG_total_depth'][sample][genome_id]) / (genome_size_of_MAG * sample_bases)) * 1e+12
        # calculate the PPB value
        total_RPKM = float(np.sum(list(final_results['MAG_RPKM_in_samples'][sample].values())))
        total_DPGB = float(np.sum(list(final_results['MAG_DPGB_in_samples'][sample].values())))
        for genome_id in genome_list:
            # TPM
            if total_RPKM > 0:
                final_results['MAG_TPM_in_samples'][sample][genome_id] = (float(final_results['MAG_RPKM_in_samples'][sample][genome_id]) / total_RPKM) * 1e+6
            else:
                final_results['MAG_TPM_in_samples'][sample][genome_id] = 0
            # PPB
            if total_DPGB > 0:
                final_results['MAG_PPB_in_samples'][sample][genome_id] = (float(final_results['MAG_DPGB_in_samples'][sample][genome_id]) / total_DPGB) * 1e+9
            else:
                final_results['MAG_PPB_in_samples'][sample][genome_id] = 0
    logger.info(f'{time_current()}   Finished MAG abundance (RPKM/TPM/DPGB/PPB) calculation between samples.')
    logger.info(f'{time_current()} Finished abundance calculation for MAG/gene_set.')
    ####################################### write final results ####################################################
    logger.info(f'{time_current()} Finally, writting all results to files...')
    result_df = dict_to_df(final_results['MAG_DPGB_in_samples'])
    result_df.to_csv(DPGB_file, sep='\t', index=True, header=True, na_rep=0, mode='w')
    logger.info(f'{time_current()}   Write DPGB results to file: {DPGB_file}.')

    result_df = dict_to_df(final_results['MAG_RPKM_in_samples'])
    result_df.to_csv(RPKM_file, sep='\t', index=True, header=True, na_rep=0, mode='w')
    logger.info(f'{time_current()}   Write RPKM results to file: {RPKM_file}.')

    result_df = dict_to_df(final_results['MAG_PPB_in_samples'])
    result_df.to_csv(PPB_file, sep='\t', index=True, header=True, na_rep=0, mode='w')
    logger.info(f'{time_current()}   Write PPB results to file: {PPB_file}.')

    result_df = dict_to_df(final_results['MAG_TPM_in_samples'])
    result_df.to_csv(TPM_file, sep='\t', index=True, header=True, na_rep=0, mode='w')
    logger.info(f'{time_current()}   Write TPM results to file: {TPM_file}.')

    result_df = dict_to_df(final_results['MAG_total_depth'])
    result_df.to_csv(depth_inf_file, sep='\t', index=True, header=True, na_rep=0, mode='w')
    logger.info(f'{time_current()}   Write depth informations of depth to file: {depth_inf_file}.')

    result_df = dict_to_df(final_results['MAG_total_counts'])
    result_df.to_csv(count_inf_file, sep='\t', index=True, header=True, na_rep=0, mode='w')
    logger.info(f'{time_current()}   Write mapped counts informations of depth to file: {count_inf_file}.')

if (check_path(f'{outdir}/concentrated_MAG_assembly.fa')) or (check_path(f'{outdir}/concentrated_MAG_assembly.fa.fai')) and (args.retain_assembly == False):
    os.system(f'rm {outdir}/concentrated_MAG_assembly.fa*')
if check_path(dir_of_analyzed_genomes) and (skip_input_check == False): os.system(f"rm -r {dir_of_analyzed_genomes}")
os.system(f"cp {sample_inf_file} {outdir}")
logger.info(f'{time_current()}   Mapping information of all samples has been written to file "{outdir}/sample_inf.txt".')
logger.info(f'{time_current()} Finished all steps, thanks for using this software. --Liping Qu')
