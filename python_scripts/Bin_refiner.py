#!/usr/bin/env python3
import argparse
import sys
import re
import os
import logging
import itertools
import pandas as pd
from datetime import datetime
from Bio import SeqIO
'''Note: Refinement only worked for bin sets from the same assembly, and binner should not change the sequence id. First, we assign the quality 
score to genomes. Second, the one in a MAG pair (share identify higher than a threshould) with lower score will be removed. Third, repetitive sequences 
(based on sequence id) will be removed from the one in a MAG pair (no identify threshould) with lower socre. Finally, MAG with score lower than threshould 
will be removed from the final MAG set. MAG pairs will be produced from all MAG results from different single binner and the combination of multiple binner 
results. E.g., the combination of metabat2 and maxbin2 contains MAGs whose contigs are supported by both of the two binners.'''
# obtain current time
def time_current():
    current_time = f"[{str(datetime.now().replace(microsecond=0))}]"
    return current_time
def calculate_genome_length(dictionary, id_list):
    length = 0
    for each_id in id_list:
        if each_id in dictionary:
            length += dictionary[each_id]
    return length
def combine_two_dics(bin_dic1, bin_dic2, length_dic, threshold):
    result_dic = {} # key is 'genome_size' and 'contig_list'
    result_dic['genome_list'] = [] # bin id list
    result_dic['genome_size'] = {} # key is bin id, value is a number
    result_dic['contig_list'] = {} # key is bin id, value is a list containing contig ids
    n = 0
    refined_bin_name = ''
    for each_bin1 in bin_dic1['genome_list']:
        set1 = set(bin_dic1['contig_list'][each_bin1]) # transfer contig id list to set
        for each_bin2 in bin_dic2['genome_list']:
            set2 = set(bin_dic2['contig_list'][each_bin2])
            combine_set = set1 & set2 # obtain the intersections of two id sets
            if len(combine_set) > 0: # if some contigs were supported by more tha one binner, this contigs may belong to a potential refine bin
                # calculate refined genome size
                temp_length = 0
                for each_seq in combine_set:
                    if each_seq in length_dic:
                        temp_length += length_dic[each_seq]
                    else:
                        print(f'ERROR, {each_seq} was not found in the contig_length dictionary.')
                        exit(1)
                if temp_length >= threshold:
                    n += 1
                    refined_bin_name = f'{prefix}' + str(n)
                    if refined_bin_name not in result_dic['genome_size']:
                        result_dic['genome_size'][refined_bin_name] = 0
                    if refined_bin_name not in result_dic['contig_list']:
                        result_dic['contig_list'][refined_bin_name] = []
                    result_dic['contig_list'][refined_bin_name] = sorted(combine_set)
                    result_dic['genome_size'][refined_bin_name] = temp_length
                    result_dic['genome_list'].append(refined_bin_name)
    return result_dic
def summary_checkm(checkm_results, com_factor=1, con_factor=1, N50_factor=1e-10):
    output_dic = {}
    output_dic['score'] = {}
    output_dic['Completeness'] = {}
    output_dic['Contamination'] = {}
    output_dic['N50'] = {}
    df = pd.read_table(checkm_results, header=0)
    for rows in df.iterrows():
        Contamination = rows[1]['Contamination']
        try:
            Completeness = rows[1]['Completeness']
        except KeyError:
            Completeness_General = rows[1]['Completeness_General']
            Completeness_Specific = rows[1]['Completeness_Specific']
            Completeness_Model_Used = rows[1]['Completeness_Model_Used']
            if Completeness_Model_Used == 'Gradient Boost (General Model)':
                Completeness = Completeness_General
            elif Completeness_Model_Used == 'Neural Network (Specific Model)':
                Completeness = Completeness_Specific
        try:
            N50_value = rows[1]['N50']
        except KeyError:
            try:
                N50_value = rows[1]['Contig_N50']
            except KeyError:
                N50_value = rows[1]['N50 (contigs)']
        score = com_factor * float(Completeness) - con_factor * float(Contamination) + N50_factor * float(N50_value)
        try:
            bin_name = rows[1]['Name']
        except KeyError:
            bin_name = rows[1]['Bin Id']
        output_dic['score'][bin_name] = score
        output_dic['Completeness'][bin_name] = float(Completeness)
        output_dic['Contamination'][bin_name] = float(Contamination)
        output_dic['N50'][bin_name] = float(N50_value)
    return output_dic
def compare_two_genomes(list1, list2, com_threshold=50, con_threshold=10):
    com1, con1, score1 = list1
    genome1_ok = False
    com2, con2, score2 = list2
    genome2_ok = False
    if (com1 >= com_threshold) and (con1 <= con_threshold): genome1_ok = True
    if (com2 >= com_threshold) and (con2 <= con_threshold): genome2_ok = True
    if (genome1_ok == True) and (genome2_ok == True):
        if score1 >= score2: return 'left'
        else: return 'right'
    elif (genome1_ok == True) and (genome2_ok == False): return 'left'
    elif (genome1_ok == False) and (genome2_ok == True): return 'right'
    elif (genome1_ok == False) and (genome2_ok == False):
        if score1 >= score2: return 'left'
        else: return 'right'

##obtain options
parser = argparse.ArgumentParser(description='''This script obtain a non-redundant MAG set from mutiple binning results. Biopython, pandas, checkm1/checkm2 are needed.''', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', help="Input bin sets. At least two floders (bin sets) should be specified and floder name cannot have duplicates.", type=str, required=True, default=False, nargs='+')
parser.add_argument('-n', help="MAGs supported by at least 2 to this number of tools will be considered. Default is 0 which means pipeline will use the number of input floders.", type=int, required=False, default=0)
parser.add_argument('-t', '--threads', help="Threads for checkm/checkm2.", type=int, required=False, default=8)
parser.add_argument('--p1', help="The comleteness factor.", type=float, required=False, default=1.0)
parser.add_argument('--p2', help="The contamination factor.", type=float, required=False, default=1.0)
parser.add_argument('--p3', help="The N50 factor.", type=float, required=False, default=1e-10)
parser.add_argument('--min_size', help="The minimum value of refined genome (bp).", type=int, required=False, default=50000)
parser.add_argument('--checkm1', help="Use checkm1 rather than checkm2. Default will use checkm2.", required=False, default=False, action='store_true')
parser.add_argument('--percentage', help="Two refined genomes with identity higher than this percentage (%%) will be consideried as duplications.", type=float, required=False, default=95)
parser.add_argument('-c', help="Min comleteness percentage in the final quality check.", type=int, required=False, default=50)
parser.add_argument('-x', help="Max contamination percentage in the final quality check.", type=int, required=False, default=10)
parser.add_argument('--prefix', help="Prefix for refined bins. Default is 'bin.', final bin wiil be 'bin.1.fa'.", type=str, required=False, default='bin.')
parser.add_argument('-o', help="Output dir. It must be different from the path containing query binning results.", type=str, required=True, default='./')
parser.add_argument('--remove_tmp', help="Remove tmp files after finish this pipeline.", required=False, default=False, action='store_true')
args = parser.parse_args()
remove_tmp = args.remove_tmp
prefix = args.prefix
number = args.n
outdir = args.o
percentage = args.percentage
threads = args.threads
min_genome = int(args.min_size)
## Set up the logger
if not os.path.exists(outdir):
    os.mkdir(outdir)
log_file = os.path.join(outdir,'Bin_refinement.log')
if os.path.exists(log_file):
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

### begin pipeline ###
logger.info(f'{time_current()} Copy initial bin sets to output dir and calculating the number of input binning results.')
input_list = args.i
max_number = len(input_list)
logger.info(f'{time_current()} There are {max_number} bin sets were given.')
if max_number == 1:
    logger.info(f'{time_current()} Only one bin set was given, no need to perform bin refinement, exit.')
    exit(1)
if number > max_number:
    logger.info(f'{time_current()} Given number is greater than the number of input bin sets, pipeline will used the number of inital binning results: {max_number}.')
elif number <= max_number:
    if number == 1:
        logger.info(f'{time_current()} Number threshold must higher than 1, pipeline will use the default value.')
    elif number > 1:
        max_number = number
    elif number == 0:
        pass

binning_results = [] # Initial binning reslts list
refined_bins = {} # it is a dic, key is name of bin/refined_bin set, value is a dic whose key is bin id and value is contig ids in this bin
id_to_seq = {} # a dic, key is sequence ids (non redundant), value is sequences itself
id_to_length = {} # a dic, key is sequence ids (non redundant), value is the length of it
tool_count = 0
code_to_source = {}
for folder_dir in input_list:
    tool = os.path.basename(folder_dir).replace('_bins', '')
    tool_count += 1
    tool_code = f'tool{tool_count}'
    code_to_source[tool_code] = tool
    logger.info(f'{time_current()} Processing binning results of {tool}...')
    binning_results.append(tool_code) # initial folders name list
    each_dir = os.path.join(outdir, tool_code) # make initial folders in the output dir for each binning results 
    if os.path.exists(each_dir):
        os.system(f'rm -r {each_dir}')
        os.mkdir(each_dir)
    else:
        os.mkdir(each_dir)
    if tool_code not in refined_bins:
        refined_bins[tool_code] = {}
        refined_bins[tool_code]['genome_list'] = []
        refined_bins[tool_code]['genome_size'] = {} # key is bin id, value is genome size (a number)
        refined_bins[tool_code]['contig_list'] = {} # key is bin id, value is contigs in the genome (a list containing contig ids)
    # read MAGs in each folder
    initial_count = 0
    for root, dirs, files in os.walk(folder_dir):
        if root != folder_dir:
            break
        for file in files:
            initial_count += 1
            refine_name = f'{prefix}' + str(initial_count)
            bin_path = os.path.join(folder_dir, file)
            exit_value = os.system(f'cp {bin_path} {each_dir}/{refine_name}.fa')
            if exit_value != 0:
                logger.info(f'{time_current()} Something wrong when copying {bin_path} to folder {each_dir}.'); exit(1)
            if refine_name not in refined_bins[tool_code]['contig_list']:
                refined_bins[tool_code]['contig_list'][refine_name] = []
            if refine_name not in refined_bins[tool_code]['genome_size']:
                refined_bins[tool_code]['genome_size'][refine_name] = 0
            refined_bins[tool_code]['genome_list'].append(refine_name)
            mag_fasta = list(SeqIO.parse(os.path.join(folder_dir, file), 'fasta'))
            for records in mag_fasta:
                refined_bins[tool_code]['contig_list'][refine_name].append(records.id)
                refined_bins[tool_code]['genome_size'][refine_name] += len(records.seq)
                if records.id not in id_to_seq:
                    id_to_seq[records.id] = records.seq
                if records.id not in id_to_length:
                    id_to_length[records.id] = len(records.seq)
    logger.info(f'{time_current()} Found {initial_count} initial bins in {folder_dir}.')

logger.info(f'{time_current()} Begin bin refinement process.')
# generate posible combination of input bin sets
list_of_old_tuples = []
for i in range(1, max_number + 1, 1):
    combination_tuple_list = [] # inital the temp set list before each loop
    for combination_tuple in itertools.combinations(binning_results, i): # 'combination_tuple' is a tuple contains the combination of input folder names
    # if number is 5, there are tuples like these: (1)...(5); (1,2)...(4,5); (1,2,3)...(3,4,5); (1,2,3,4)...(2,3,4,5); (1,2,3,4,5)
        combination_set = set(combination_tuple)
        if i == 1:
            logger.info(f'{time_current()} Storage single combination_set {combination_tuple[0]}.')
            combination_tuple_list.append(combination_tuple) # storage set combination_tuples in the first loop
            continue
        # if i > 1, next lines will be performed
        logger.info(f"{time_current()} Combining following binning results: {combination_tuple}.")
        combination_tuple_list.append(combination_tuple)
        new_name = '_'.join(combination_tuple) # new name is a string, e.g., if a set is {1,2}, new name will be "1_2"
        for each_old_tuple in list_of_old_tuples:
            if set(each_old_tuple) < combination_set: # if a old set is the subset of a combination_set, e.g., old set is {1,2} and combination_set is {1,2,3}
                # then left_set_id is "1_2" and right_set_id is "3"
                left_dic_id = '_'.join(each_old_tuple)
                right_dic_id = ''.join(str(a) for a in (combination_set ^ set(each_old_tuple))) # both "combination_set ^ each_old_set" or "combination_set - each_old_set" are ok here
                break # exit loop to avoid repetition, e.g., both {1,2}, {1,3} and {2,3} are the subset of {1,2,3},
                      # if we want to get refined result "1_2_3", we only need combine "1_2" and "3".
        if new_name not in refined_bins:
            refined_bins[new_name] = {}
            refined_bins[new_name]['genome_list'] = []
            refined_bins[new_name]['genome_size'] = {}
            refined_bins[new_name]['contig_list'] = {}
        temp_combine_list = {}
        if new_name not in temp_combine_list:
            temp_combine_list[new_name] = {} # the value is two dics, 'contig_list' (value is a list) and 'genome_size' (value is a number)
        # refined_bins[left_dic_id]['contig_list'], key is genome id, value is a list contiaining contig ids
        temp_combine_list[new_name] = combine_two_dics(refined_bins[left_dic_id], refined_bins[right_dic_id], id_to_length, threshold=min_genome)
        refined_bins[new_name]['contig_list'] = temp_combine_list[new_name]['contig_list'] # dic
        refined_bins[new_name]['genome_size'] = temp_combine_list[new_name]['genome_size'] # dic
        refined_bins[new_name]['genome_list'] = temp_combine_list[new_name]['genome_list'] # list
        # write refned bins to folder
        refined_bin_ids = refined_bins[new_name]['genome_list']
        if len(refined_bin_ids) > 0:
            new_dir = os.path.join(outdir, new_name)
            if os.path.exists(new_dir):
                logger.info(f'{time_current()} Output dir for {new_dir} already exists, it will be removed.')
                os.system(f'rm -r {new_dir}')
                os.mkdir(new_dir)
            else:
                os.mkdir(new_dir)
            for each_bin in refined_bin_ids:
                write_path = os.path.join(new_dir, f'{each_bin}.fa')
                temp_contigs = refined_bins[new_name]['contig_list'][each_bin]
                with open(write_path, 'w') as temp_out:
                    for contig in temp_contigs:
                        temp_out.writelines(f'>{contig}\n{id_to_seq[contig]}\n')
            logger.info(f"{time_current()} {len(refined_bin_ids)} refined bins were generated from following binning results: {combination_tuple}.")
        elif len(refined_bin_ids) == 0:
            logger.info(f"{time_current()} No refined bins were generated from following binning results: {combination_tuple}.")
            del refined_bins[new_name]
    list_of_old_tuples = combination_tuple_list # storage set combination_sets in this loop, at the same time, old set combination_sets will be removed. It is a list contains multiple sets
logger.info(f'{time_current()} Obatin bin quality of each bin refinement.')
final_refined_combinations = refined_bins.keys() # refinement list, e.g., ['A', 'B_C']...
for each_refinement in final_refined_combinations:
    logger.info(f'{time_current()} Running checkm for {each_refinement}...')
    checkm_in = os.path.join(outdir, each_refinement)
    checkm_out = os.path.join(outdir, f'{each_refinement}.checkm')
    checkm_tmp = os.path.join(outdir, f'{each_refinement}.tmp')
    checkm_file = os.path.join(checkm_out, 'quality_report.tsv')
    if os.path.exists(checkm_file):
        logger.info(f'{time_current()} Checkm step for {each_refinement} was already done, run next step.')
    else:
        if not os.path.exists(checkm_tmp): os.mkdir(checkm_tmp)
        if os.path.exists(checkm_out): os.system(f'rm -r {checkm_out}')
        if args.checkm1 == False:
            exit_value =  os.system(f'checkm2 predict --quiet --input {checkm_in} --output-directory {checkm_out} -x fa \
                                    --threads {threads} --remove_intermediates --tmpdir {checkm_tmp}')
        elif args.checkm1 == True:
            exit_value =  os.system(f'checkm lineage_wf --quiet -x fa --pplacer_threads {threads} --threads {threads} --tmpdir {checkm_tmp} {checkm_in} {checkm_out} && \
                                    checkm qa -o 2 -f {checkm_file} --tab_table -t {threads} {checkm_out}/lineage.ms {checkm_out}')
        if exit_value != 0:
            logger.info(f'{time_current()} Something wrong when running checkm2 for combination {each_refinement}.'); exit(1)

score_dic = {}
for each_refinement in final_refined_combinations:
    checkm_out = os.path.join(outdir, f'{each_refinement}.checkm')
    checkm_file = os.path.join(checkm_out, 'quality_report.tsv')
    logger.info(f'{time_current()} Reading checkm results of {each_refinement}...')
    if each_refinement not in score_dic:
        score_dic[each_refinement] = {}
    score_dic[each_refinement] = summary_checkm(checkm_file, com_factor=float(args.p1), con_factor=float(args.p2), N50_factor=float(args.p3)) # a dic, key is refinement, value is for dics including score, com, con, N50
logger.info(f'{time_current()} Selecting best bin refinement results, performing deduplication...')
logger.info(f'{time_current()} Calculating the overlaping between bin pairs...')
temp_final_refinements = {} # contains the best version of bin paris, but there are still some duplications among bins
i = 0
for each_refinement in final_refined_combinations:
    if len(temp_final_refinements) == 0:
        logger.info(f'{time_current()} Add the first refinement to the analysis data set.')
        # initial the temp final results. The old version will be replaced if a better version was found. genomes in the right refinement and
        # not found in left refinement will also be added to the temp final result
        temp_final_refinements['combination'] = {}
        temp_final_refinements['genome_list'] = refined_bins[each_refinement]['genome_list'] # genome list, it will only increase or remain unchanged, not decrease
        temp_final_refinements['contig_list'] = refined_bins[each_refinement]['contig_list'] # dic
        temp_final_refinements['genome_size'] = refined_bins[each_refinement]['genome_size'] # dic
        temp_final_refinements['score'] = score_dic[each_refinement]['score'] # dic
        temp_final_refinements['Completeness'] = score_dic[each_refinement]['Completeness'] # dic
        temp_final_refinements['Contamination'] = score_dic[each_refinement]['Contamination'] # dic
        temp_final_refinements['N50'] = score_dic[each_refinement]['N50'] # dic
        for each_left_bin in temp_final_refinements['genome_list']:
            temp_final_refinements['combination'][each_left_bin] = str(each_refinement)
        logger.info(f'{time_current()} Start binning deduplication.')
    elif len(temp_final_refinements) > 0:
        temp_duplication = []
        left_bins = temp_final_refinements['genome_list']
        right_bins = refined_bins[each_refinement]['genome_list']
        for each_left_bin in left_bins:
            left_contigs = temp_final_refinements['contig_list'][each_left_bin]
            left_size = temp_final_refinements['genome_size'][each_left_bin]
            left_score = temp_final_refinements['score'][each_left_bin]
            left_N50 = temp_final_refinements['N50'][each_left_bin]
            left_completeness = temp_final_refinements['Completeness'][each_left_bin]
            left_contimination = temp_final_refinements['Contamination'][each_left_bin]
            left_combination = temp_final_refinements['combination'][each_left_bin]
            for each_right_bin in right_bins:
                i += 1
                right_contigs = refined_bins[each_refinement]['contig_list'][each_right_bin]
                right_size = refined_bins[each_refinement]['genome_size'][each_right_bin]
                right_score = score_dic[each_refinement]['score'][each_right_bin]
                right_completeness = score_dic[each_refinement]['Completeness'][each_right_bin]
                right_contimination = score_dic[each_refinement]['Contamination'][each_right_bin]
                right_N50 = score_dic[each_refinement]['N50'][each_right_bin]
                right_combination = str(each_refinement)
                alignment = set(left_contigs) & set(right_contigs)
                overlap_length = calculate_genome_length(id_to_length, list(alignment))
                left_ratio = 100 * (overlap_length / left_size)
                right_ratio = 100 * (overlap_length / right_size)
                overlap_ratio = max([left_ratio, right_ratio])
                if overlap_ratio >= percentage: # determine whether the right genome is the duplication of left genome
                    # whatever better version or not, we need add it to the list 'temp_duplication'.
                    if each_right_bin not in temp_duplication:
                        temp_duplication.append(each_right_bin)
                    better_genome = compare_two_genomes([left_completeness, left_contimination, left_score], 
                                                    [right_completeness, right_contimination, right_score], 
                                                    com_threshold=50, con_threshold=10)
                    if better_genome == 'right': # determine whether the right genome is better than the left genome
                        # replace the old version (left) information using better version (right) result, but do not change the genome id
                        # print(f'{each_right_bin} in {right_combination} with score {right_score} is better than {each_left_bin} in {left_combination} with {left_score}')
                        temp_final_refinements['contig_list'][each_left_bin] = right_contigs
                        temp_final_refinements['genome_size'][each_left_bin] = right_size
                        temp_final_refinements['combination'][each_left_bin] = right_combination
                        temp_final_refinements['score'][each_left_bin] = right_score
                        temp_final_refinements['Completeness'][each_left_bin] = right_completeness
                        temp_final_refinements['Contamination'][each_left_bin] = right_contimination
                        temp_final_refinements['N50'][each_left_bin] = right_N50
        if len(temp_duplication) < len(right_bins):# if Ture, it means that some genomes in right refinement was not found a dulication in left refinement.
            for each_right_bin in right_bins:
                if each_right_bin not in temp_duplication:
                    count = len(temp_final_refinements['genome_list']) + 1 # add non-duplication genomes in right refinement to the end of the temp final results
                    remain_name = f'{prefix}' + str(count)
                    right_contigs = refined_bins[each_refinement]['contig_list'][each_right_bin]
                    right_size = refined_bins[each_refinement]['genome_size'][each_right_bin]
                    right_score = score_dic[each_refinement]['score'][each_right_bin]
                    right_completeness = score_dic[each_refinement]['Completeness'][each_right_bin]
                    right_contimination = score_dic[each_refinement]['Contamination'][each_right_bin]
                    right_N50 = score_dic[each_refinement]['N50'][each_right_bin]
                    # print(f'{each_right_bin} in {right_combination} is not a duplication and was storage with new name {remain_name}')
                    temp_final_refinements['contig_list'][remain_name] = right_contigs
                    temp_final_refinements['genome_size'][remain_name] = right_size
                    temp_final_refinements['combination'][remain_name] = str(each_refinement)
                    temp_final_refinements['score'][remain_name] = right_score
                    temp_final_refinements['Completeness'][remain_name] = right_completeness
                    temp_final_refinements['Contamination'][remain_name] = right_contimination
                    temp_final_refinements['N50'][remain_name] = right_N50
                    temp_final_refinements['genome_list'].append(remain_name)
logger.info(f'{time_current()} Finished MAG deduplication, processed {i} pairs of genomes.')
logger.info(f'{time_current()} Removing repeat sequences from refined bins, sequence will be removed from the bins with lower socre...')
final_refinement = {}
i = 0
for combination_tuple in itertools.combinations(temp_final_refinements['genome_list'], 2):
    left_genome = combination_tuple[0]
    left_N50 = temp_final_refinements['N50'][left_genome]
    left_completeness = temp_final_refinements['Completeness'][left_genome]
    left_contimination = temp_final_refinements['Contamination'][left_genome]
    left_score = temp_final_refinements['score'][left_genome]
    left_contigs = temp_final_refinements['contig_list'][left_genome]
    right_genome = combination_tuple[1]
    right_N50 = temp_final_refinements['N50'][right_genome]
    right_completeness = temp_final_refinements['Completeness'][right_genome]
    right_contimination = temp_final_refinements['Contamination'][right_genome]
    right_score = temp_final_refinements['score'][right_genome]
    right_contigs = temp_final_refinements['contig_list'][right_genome]
    i += 1
    # if two genomes have intersection contigs, remove these contigs from the genome with lower score
    intersection_contigs = set(left_contigs) & set(right_contigs)
    if len(intersection_contigs) > 0:
        better_genome = compare_two_genomes([left_completeness, left_contimination, left_score], 
                                    [right_completeness, right_contimination, right_score], 
                                    com_threshold=50, con_threshold=10)
        if better_genome == 'right': # if right genome is better, remove repeat contig s from left genome
            left_contigs = list(filter(lambda x:x not in intersection_contigs, left_contigs))
            temp_final_refinements['contig_list'][left_genome] = left_contigs
        elif better_genome == 'left':
            right_contigs = list(filter(lambda x:x not in intersection_contigs, right_contigs))
            temp_final_refinements['contig_list'][right_genome] = right_contigs
logger.info(f"{time_current()} Finished removing duplicate sequences, processed {i} pairs of genomes.")
logger.info(f"{time_current()} Write clean genomes to folder, calculating newest genome size, genomes with no contigs will be discared.")
clean_bins_dir = os.path.join(outdir, 'refined_bins_temp')
if not os.path.exists(clean_bins_dir):
    os.mkdir(clean_bins_dir)
else:
    os.system(f'rm -r {clean_bins_dir}')
    os.mkdir(clean_bins_dir)
temp_final_refinements['genome_size'] = {}
temp_final_refinements['score'] = {}
for each_genome in temp_final_refinements['genome_list']:
    if len(temp_final_refinements['contig_list'][each_genome]) == 0: continue
    with open(f'{clean_bins_dir}/{each_genome}.fa', 'w') as tempout:
        for each_seq in temp_final_refinements['contig_list'][each_genome]:
            tempout.writelines(f'>{each_seq}\n{id_to_seq[each_seq]}\n')
logger.info(f"{time_current()} Rerun checkm for clean bins...")
checkm_in_final = clean_bins_dir
checkm_out_final = os.path.join(outdir, 'refined_bins_temp.checkm')
checkm_tmp_final = os.path.join(outdir, 'refined_bins_temp.tmp')
checkm_file_final = os.path.join(checkm_out_final, 'quality_report.tsv')

if os.path.exists(checkm_file_final):
    logger.info(f'{time_current()} The checkm step was already done for the final clean bins.')
else:
    if not os.path.exists(checkm_tmp_final): os.mkdir(checkm_tmp_final)
    if os.path.exists(checkm_out_final): os.system(f'rm -r {checkm_out_final}')
    if args.checkm1 == False:
        exit_value =  os.system(f'checkm2 predict --quiet --input {checkm_in_final} --output-directory {checkm_out_final} -x fa \
                                --threads {threads} --remove_intermediates --tmpdir {checkm_tmp_final}')
    elif args.checkm1 == True:
        exit_value =  os.system(f'checkm lineage_wf --quiet -x fa --pplacer_threads {threads} --threads {threads} --tmpdir {checkm_tmp_final} {checkm_in_final} {checkm_out_final} && \
                                checkm qa -o 2 -f {checkm_file_final} --tab_table -t {threads} {checkm_out_final}/lineage.ms {checkm_out_final}')
    if exit_value != 0:
        logger.info(f'{time_current()} Something wrong when running checkm2 for clean bins.'); exit(1)
    os.system(f'cp {checkm_file_final} {outdir}/refined_bins_temp.txt')
logger.info(f"{time_current()} Check quality of clean bins, remove bad bins.")
logger.info(f"{time_current()} Please NOTE, genomes with qualified compleness and contamination will be retained even their genome size is lower that threshold.")
Final_quality = os.path.join(outdir, 'Final_quality_of_refined_bins.txt')
Final_quality_write = open(Final_quality, 'w')
final_counts = 0
final_genomes = {}
clean_bins_dir_final = os.path.join(outdir, 'refined_bins_final')
if not os.path.exists(clean_bins_dir_final):
    os.mkdir(clean_bins_dir_final)
else:
    os.system(f'rm -r {clean_bins_dir_final}')
    os.mkdir(clean_bins_dir_final)
with open(checkm_file_final, 'r') as tempin:
    lines = tempin.read().splitlines()
for line in lines:
    if 'Contamination' in line:
        Final_quality_write.writelines(f'{line}\n') # print header of checkm result
    else:
        cut = re.split('\t', line)
        genome_id = cut[0]
        contamination = cut[2]
        if cut[4] == 'Neural Network (Specific Model)':
            completeness = cut[3]
        else:
            completeness = cut[1]
        if (float(contamination) <= args.x) and (float(completeness) >= args.c):
            final_counts += 1
            new_id = f'{prefix}' + str(final_counts)
            os.system(f'cp {clean_bins_dir}/{genome_id}.fa {clean_bins_dir_final}/{new_id}.fa')
            new_line = re.sub(genome_id, new_id, line)
            Final_quality_write.writelines(f'{new_line}\n')
            final_genomes[genome_id] = new_id
Final_quality_write.close()
if remove_tmp == True:
    for each_refinement in final_refined_combinations:
        checkm_in = os.path.join(outdir, each_refinement)
        checkm_out = os.path.join(outdir, f'{each_refinement}.checkm')
        checkm_tmp = os.path.join(outdir, f'{each_refinement}.tmp')
        os.system(f'rm -r {checkm_in}')
        os.system(f'rm -r {checkm_tmp}')
        os.system(f'rm -r {checkm_out}')
    os.system(f'rm -r {checkm_tmp_final}')
    os.system(f'rm -r {checkm_out_final}')
logger.info(f"{time_current()} Write source of clean bins.")
with open(f'{outdir}/Source_of_final_refined_bins.txt', 'w') as tempout:
    tempout.writelines('Bin Id\tCombination\n')
    for genome_id in sorted(final_genomes.keys()):
        temp_code = temp_final_refinements['combination'][genome_id]
        temp_code = re.sub('_', '__', temp_code)
        for each_code in code_to_source.keys():
            temp_code = re.sub(each_code, code_to_source[each_code], temp_code)
        tempout.writelines(f"{final_genomes[genome_id]}\t{temp_code}\n")

logger.info(f"{time_current()} Finished all steps, a total of {final_counts} refined bins were generated.")
logger.info(f"{time_current()} Source of refined bins were written to {outdir}/Source_of_final_refined_bins.txt.")
logger.info(f"{time_current()} Quality of refined bins were written to {Final_quality}.")
logger.info(f"{time_current()} Raw refined genomes before checkm check were retained to {outdir}/refined_bins_temp.")
logger.info(f"{time_current()} Quality of raw refined genomes were retained to {outdir}/refined_bins_temp.txt.")
