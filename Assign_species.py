#!/usr/bin/env python3
import argparse
import subprocess
import os, re
import pandas as pd
from datetime import datetime
# obtain current time
def time_current():
    current_time = f"[{str(datetime.now().replace(microsecond=0))}]"
    return current_time
##obtain options
parser = argparse.ArgumentParser(description="This scripts assign species ids based on the gtdb-tk results and dRep.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-d', help="Genomes dir.", type=str, required=True)
parser.add_argument('-x', help="Suffix of genomes.", type=str, required=False, default='fa')
parser.add_argument('-o', help="Output dir.", type=str, required=True)
parser.add_argument('-t', help="Threads for dRep.", type=int, required=False, default=8)
parser.add_argument('--gtdb', help="Input gtdb-tk results.", type=str, required=True)
parser.add_argument('--genome_column', help="The genome id column name in the gtdb-tk result.", type=str, required=False, default='user_genome')
parser.add_argument('--taxa_column', help="The taxonomy column name in the gtdb-tk result.", type=str, required=False, default='classification')
parser.add_argument('--type', help="Choose gtdb or dRep based when there are conflicts.", required=False, default='dRep', choices=['gtdb', 'dRep'])
args = parser.parse_args()
threads = args.t
# run dRep
genomes_dir = args.d
suffix = args.x
if suffix[0] == '.':
    suffix = suffix[1:]
outdir = args.o
if not os.path.exists(outdir):
    os.mkdir(outdir)

taxa_dict = {'d__':'domain', 'p__':'phylum', 'c__':'class', 'o__':'order', 'f__':'family', 'g__':'genus', 's__':'species'}
short_levels = list(taxa_dict.keys())
long_levels = list(taxa_dict.values())

# read gtdb results
print(f'{time_current()} Read the gtdb-tk result.')
genome_column = args.genome_column
taxa_colum =  args.taxa_column
gtdb_df = pd.read_table(args.gtdb, header=0, usecols=[genome_column, taxa_colum])
gtdb_df = gtdb_df.sort_values(by=taxa_colum, ascending=True)# it will cause incomplete taxa to be arranged above the complete taxa
archaea_df = gtdb_df[gtdb_df[taxa_colum].str.contains('d__Archaea', na=False)]
bacteria_df = gtdb_df[gtdb_df[taxa_colum].str.contains('d__Bacteria', na=False)]
print(f'{time_current()} There are {len(archaea_df)} archaeal and {len(bacteria_df)} bacterial MAGs, respectively.')
if len(archaea_df) == 0 or len(bacteria_df) == 0:
    print(f'{time_current()} No available genomes were found in the gtdb-tk results, exit.'); exit(1)
gtdb_df = gtdb_df.set_index(genome_column)

print(f'{time_current()} Perform species level cluster based on ANI.')
# get absulote path of bacteria and archaea genome files, make two list file for dRep
print(f'{time_current()}   Make file list for dRep.')
all_genomes = [] # since genome list in gtdb rsult may be redundant, we need get non-redundant genome list whih also be used in dRep.
archaea_list = os.path.join(outdir, '00_archaea_list.txt')
bacteria_list = os.path.join(outdir, '00_bacteria_list.txt')
if len(archaea_df) > 0:
    with open(archaea_list, 'w') as aro:
        for each in archaea_df[genome_column].tolist():
            bin_path = f"{os.path.join(genomes_dir, each + '.' + suffix)}"
            if os.path.exists(bin_path):
                aro.writelines(f"{os.path.join(genomes_dir, each + '.' + suffix)}\n")
                all_genomes.append(each)
            else:
                print(f'{time_current()}   Genome {each}.{suffix} was exists in the gtdb-tk result but not found in the MAG folder, skip it.')
if len(bacteria_df) > 0:
    with open(bacteria_list, 'w') as bao:
        for each in bacteria_df[genome_column].tolist():
            bin_path = f"{os.path.join(genomes_dir, each + '.' + suffix)}"
            if os.path.exists(bin_path):
                bao.writelines(f"{os.path.join(genomes_dir, each + '.' + suffix)}\n")
                all_genomes.append(each)
            else:
                print(f'{time_current()}   Genome {each}.{suffix} was exists in the gtdb-tk result but not found in the MAG folder, skip it.')

# run dRep, get dRep results
def get_name(text_in):
    return os.path.splitext(os.path.basename(text_in))[0]
def run_dRep(list_file, dRep_dir):
    global threads
    cdb_result = os.path.join(dRep_dir, 'data_tables', 'Cdb.csv')
    ndb_result = os.path.join(dRep_dir, 'data_tables', 'Ndb.csv')
    if not os.path.exists(dRep_dir):
        os.mkdir(dRep_dir)
    cmd = f"dRep dereplicate -g {list_file} -p {threads} --S_algorithm skani -sa 0.95 -pa 0.9 -comp 0 -con 100 {dRep_dir}"
    if (not os.path.exists(cdb_result)) or (not os.path.exists(ndb_result)):
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, shell=True)
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                print(output.strip())
        rc = process.wait()
        if rc != 0:
            raise subprocess.CalledProcessError(rc, cmd)
        else:
            print(f'{time_current()}     dRep runs successfully.')
    print(f'{time_current()}     Read dRep results...')
    dRep_df = pd.read_csv(cdb_result, header=0, usecols=['genome', 'secondary_cluster'])
    dRep_df['genome'] = dRep_df['genome'].apply(get_name)
    ani_df = pd.read_csv(ndb_result, header=0, usecols=['reference', 'querry', 'ani'])
    ani_df['reference'] = ani_df['reference'].apply(get_name)
    ani_df['querry'] = ani_df['querry'].apply(get_name)
    ani_df = ani_df[ani_df['querry'] != ani_df['reference']] # remove the self comparison
    return [dRep_df, ani_df]
def add_prefix(input_text, prefix):
    return str(prefix) + str(input_text)
total_dRep_df = pd.DataFrame()
total_ani_df = pd.DataFrame()
if len(archaea_df) > 0:
    print(f'{time_current()}   Running dRep for archaeal MAGs...')
    dRep_dir = os.path.join(outdir, '00_dRep_dir_archaea')
    archaea_drep_df, archaea_ani_df = run_dRep(archaea_list, dRep_dir)
    archaea_drep_df['secondary_cluster'] = archaea_drep_df['secondary_cluster'].apply(add_prefix, args=['arc_sp.'])
    archaea_species = archaea_drep_df['secondary_cluster']
    total_dRep_df = archaea_drep_df
    total_ani_df = archaea_ani_df
if len(bacteria_df) > 0:
    print(f'{time_current()}   Running dRep for bacterial MAGs...')
    dRep_dir = os.path.join(outdir, '00_dRep_dir_bacteria')
    bacteria_drep_df, bacteria_ani_df = run_dRep(bacteria_list, dRep_dir)
    bacteria_drep_df['secondary_cluster'] = bacteria_drep_df['secondary_cluster'].apply(add_prefix, args=['bac_sp.'])
    bacteria_species = bacteria_drep_df['secondary_cluster']
    if len(total_dRep_df) > 0:
        total_dRep_df = pd.concat([total_dRep_df, bacteria_drep_df], ignore_index=True)
        total_ani_df = pd.concat([total_ani_df, bacteria_ani_df], ignore_index=True)
    elif len(total_dRep_df) == 0:
        total_dRep_df = bacteria_drep_df
        total_ani_df = bacteria_ani_df

## assign taxonomy for genomes
'''
If a genome have unclear gtdb result, we will find genomes in the same species cluster of query genome, the taxonomy of target genomes will be assigned to query genome. 
If multiple target taxa are found, the taxa of the genome share the highest ANI will be used
    first, assign lineage for genomes with uncomplete lineage, e.g., gtdb taxonomy of 'bin.69.17' is: "d__Bacteria;p__JADFOP01;c__JADFOP01;o__JADFOP01;f__JADFOP01;g__;s__", 
while it share >95% ANI with "d__Bacteria;p__JADFOP01;c__JADFOP01;o__JADFOP01;f__JADFOP01;g__JADFOP01;s__JADFOP01 sp022562795" members, in other
word, 'bin.69.17' is not a new genus, or even not a new species.
    second, assign genus and species ID for genomes.
Note: bin.GT1.13 was identified as g__JACZQL01, but it shared higher ANI with g__JACZXR01 members. It cause a problem, that some members with incomplete lineage
in g__JACZXR01 share higheast ANI with bin.GT1.13, and these members may be wrongly assigned as g__JAZQL01.
'''
def compare_lineage_length(input_lineage_list:list):
    global short_levels
    length_dict = {}
    for each_lineage in input_lineage_list:
        n = 0
        for each in re.split(';', each_lineage):
            if each not in short_levels: n += 1
        if n not in length_dict:
            length_dict[n] = [each_lineage]
        else:
            length_dict[n] = length_dict[n] + [each_lineage]
    max_length = max(length_dict.keys())
    return length_dict[max_length]
print(f'{time_current()} Complete the lineage of MAGs based on the other genomes within the species cluster.')
taxa_change_log_df = pd.DataFrame(columns=['original_taxonomy', 'new_taxonomy', 'original_dRep_id', 'species_ids', 'lineage_source', 'Modified_lineage', 'Uncharacterized_species', 'New_species', 'NOTE'], index=all_genomes, data='')
# taxa_change_log_df['NOTE'] = taxa_change_log_df['NOTE'].astype(str)
cluster_list = sorted(list(set(total_dRep_df['secondary_cluster'].tolist())))
for each_cluster in cluster_list: # scan each species cluster in all clusters
    genomes_dRep_df = total_dRep_df[total_dRep_df['secondary_cluster'] == each_cluster]
    genomes = genomes_dRep_df['genome'].tolist() # genomes in the species cluster
    genomes_taxa_df = gtdb_df[gtdb_df.index.isin(genomes)] # gtdb results df of the species cluster
    genomes_ani_df = total_ani_df[(total_ani_df['reference'].isin(genomes)) | (total_ani_df['querry'].isin(genomes))]
    # check the taxa of genomes in the species cluster, if all have complete lineage,
    if len(genomes) == 1:
        taxa_change_log_df.loc[genomes[0],'new_taxonomy'] = taxa_change_log_df.loc[genomes[0],'original_taxonomy'] = gtdb_df.loc[genomes[0], taxa_colum]
        taxa_change_log_df.loc[genomes[0],'original_dRep_id'] = each_cluster
        taxa_change_log_df.loc[genomes[0],'NOTE'] += f'Only one genome was found in the species cluster. '
        taxa_change_log_df.loc[genomes[0],'lineage_source'] = 'original'
    else:
        for genome in genomes:
            taxa_change_log_df.loc[genome,'original_dRep_id'] = each_cluster
        ## found the reference complete lineage for incomplete lineages.remove('s__')
        total_lineage = list(set(genomes_taxa_df[taxa_colum].tolist())) # non-redundant lineage list
        if len(total_lineage) == 1:
            # all genomes have the most complete lineages
            for genome in genomes:
                taxa_change_log_df.loc[genome,'new_taxonomy'] = taxa_change_log_df.loc[genome,'original_taxonomy'] = genomes_taxa_df.loc[genome, taxa_colum]
                taxa_change_log_df.loc[genome,'NOTE'] += f'The genome already have the most complete lineage in the species. '
                taxa_change_log_df.loc[genome,'lineage_source'] = 'original'
        else: # found different lineage in the species, 1) there may be some uncomplete lineage (same or different), and complete lineages may be the same or different; 2) different taxa with the same length were clustered together
            # found the most complete lineage
            complete_lineages = compare_lineage_length(total_lineage)
            complete_lineage_df = genomes_taxa_df[genomes_taxa_df[taxa_colum].isin(complete_lineages)]
            complete_lineage_genomes = complete_lineage_df.index.tolist()
            for genome in genomes:
                old = taxa_change_log_df.loc[genome,'original_taxonomy'] = genomes_taxa_df.loc[genome, taxa_colum]
                if len(complete_lineages) == 1:
                    # only one kind of most complete lineage and some incomplete lineages were found
                    complete_lineage = complete_lineages[0]
                    new = taxa_change_log_df.loc[genome,'new_taxonomy'] = complete_lineage
                    if old == new:
                        taxa_change_log_df.loc[genome,'NOTE'] += f'The genome already have the most complete lineage in the species. '
                        taxa_change_log_df.loc[genome,'lineage_source'] = 'original'
                    else:
                        taxa_change_log_df.loc[genome,'NOTE'] += f'The lineage was modified. '
                        taxa_change_log_df.loc[genome,'lineage_source'] = 'The most complete lineage in the cluster.'
                elif len(complete_lineages) > 1:
                    # 1) some uncomplete lineage and genomes with different complete lineages:
                    # if all complete lineages have the same or different common ansestor with the incomplete lineage, only use ANI based selection
                    # if parts of complete lineages have the same common ansestor with the incomplete lineage,, first extract complete lineage with the same acsestor of incomplete lineage, then perform ANI if mutiple lineages were found or modified incomplete lineage if single lineage was found
                    # 2) only genomes with different complete lineages, in this situation, do not change original lineage
                    if old not in complete_lineages: # it means the old lineage should be the incomplete lineage
                        # extract ani df of genomes with most complete lineages
                        ref_ani_df = genomes_ani_df[((genomes_ani_df['reference'].isin(complete_lineage_genomes)) & (genomes_ani_df['querry'] == genome)) | ((genomes_ani_df['querry'].isin(complete_lineage_genomes)) & (genomes_ani_df['reference'] == genome))]
                        max_ani = ref_ani_df['ani'].max()
                        max_ani_df = ref_ani_df[ref_ani_df['ani'] == max_ani]
                        close_ref = list(set(max_ani_df['querry'].tolist() + max_ani_df['reference'].tolist()))
                        close_ref.remove(genome)
                        if len(close_ref) == 1:
                            close_ref = close_ref[0] # extract the single element in the list
                        elif len(close_ref) > 1:
                            # see if the taxa of mutiple close genomes is the same or not
                            close_taxa_df = genomes_taxa_df[genomes_taxa_df.index.isin(close_ref)] # gtdb result df of closeast genomes
                            temp_taxa_list = list(set(close_taxa_df[taxa_colum].tolist()))
                            if len(temp_taxa_list) == 1:
                                taxa_change_log_df.loc[genome,'NOTE'] += f'Reminder: This genome shares the same higheast ANI with multiple genomes with the same taxonomy. '
                            else:
                                taxa_change_log_df.loc[genome,'NOTE'] += f'WARNING: This genome shares the same higheast ANI with multiple genomes with different taxonomy, they are {temp_taxa_list} and we chose the first one. '
                            close_ref = close_ref[0]
                        else:
                            print(f'{time_current()} Some thing wrong when runing ANI summary.')
                            exit(1)
                        new = taxa_change_log_df.loc[genome,'new_taxonomy'] = genomes_taxa_df.loc[close_ref, taxa_colum]
                        taxa_change_log_df.loc[genome,'NOTE'] += f'The lineage was modified. '
                        taxa_change_log_df.loc[genome,'lineage_source'] = close_ref
                    else: # the lineage of the genome is one of the most complete lineages
                        taxa_change_log_df.loc[genome,'new_taxonomy'] = genomes_taxa_df.loc[genome, taxa_colum]
                        taxa_change_log_df.loc[genome,'NOTE'] += f'WARNING: This genome already have the most complete lineage in the species, but different species were found in the single cluster. '
                        taxa_change_log_df.loc[genome,'lineage_source'] = 'original'
                else:
                    print(f'{time_current()} Did not found complete lineages.')
                    exit(1)
## assign genus, change the raw id from dRep. e.g., arc_sp.2_1, means this genome belongs to the first species in the second genus in the Archaea domain
print(f'{time_current()} Summary the genus and species information, assgin the genus and species IDs.')
# scan the taxonomy assignment df, assign genus and species IDs for genomes
# 'species' id is the original dRep ids,
# for genus, if MAG have compelte lineage at genus level, genus is the genus level lineage, if not, use the primary cluster
# if MASH (90% ANI). e.g., CSP1-5 members have multiple primary cluster (like 19_1, 18_1...) but they still be assigned the same genus level id 'bac_sp.20_X'
# bin.SY51.8, bin.SY28.25, bin.GT1.8, bin.SY9.2, bin.SY3.11 have incomplete lineages (e.g., d__Bacteria;p__Pseudomonadota;c__;o__;f__;g__;s__), 
# and the phylogenetic position of them is different in the phylogenetic tree, and they have different primary cluster id (e.g., 27_1, 28-1), 
# therefore they should assigned with different genus ids
taxa_change_log_df = taxa_change_log_df.sort_values(by=['new_taxonomy', 'original_dRep_id'], ascending=[True, True])
old_genus = ''; old_species = ''; a_gn = 0; a_sn = 0; b_gn = 0; b_sn = 0
for index, row in taxa_change_log_df.iterrows():
    domain, phylumn, class_, order, family, genus, statu = re.split(';', row['new_taxonomy'])
    if statu == 's__':
        taxa_change_log_df.loc[index, 'New_species'] = 'Yes'
        taxa_change_log_df.loc[index, 'Uncharacterized_species'] = 'Yes'
    else:
        taxa_change_log_df.loc[index, 'New_species'] = 'No'
        if ' sp' in statu:
            temp_statu = re.split(' sp', statu, 1)[1]
            if str(temp_statu).isnumeric():
                taxa_change_log_df.loc[index, 'Uncharacterized_species'] = 'Yes'
            else:
                taxa_change_log_df.loc[index, 'Uncharacterized_species'] = 'No'
        else:
            taxa_change_log_df.loc[index, 'Uncharacterized_species'] = 'No'
    original_d = domain
    genus_level_lineage = ';'.join([domain, phylumn, class_, order, family, genus])
    species = row['original_dRep_id']
    if domain == 'd__Archaea':
        domain = 'arc_sp.'
        if old_genus != genus_level_lineage:
            old_genus = genus_level_lineage; a_gn += 1; a_sn = 1; old_species = species # initialize species id for each genus
        elif (old_genus == genus_level_lineage) and (old_species != species):
            if genus_level_lineage[-3:] == 'g__':
                a_gn += 1
            a_sn += 1; old_species = species
        taxa_change_log_df.loc[index,'species_ids'] = str(domain) + str(a_gn) + '_' + str(a_sn)
    elif domain == 'd__Bacteria':
        domain = 'bac_sp.'
        if old_genus != genus_level_lineage:
            b_gn += 1; b_sn = 1; old_genus = genus_level_lineage; old_species = species
        elif (old_genus == genus_level_lineage) and (old_species != species):
            if genus_level_lineage[-3:] == 'g__':
                b_gn += 1
            b_sn += 1; old_species = species
        taxa_change_log_df.loc[index,'species_ids'] = str(domain) + str(b_gn) + '_' + str(b_sn)
    new_species = 's__' + str(taxa_change_log_df.loc[index,'species_ids'])
    taxa_change_log_df.loc[index,'Modified_lineage'] = ';'.join([original_d, phylumn, class_, order, family, genus, new_species])
taxa_change_log_df = taxa_change_log_df.reset_index(drop=False).rename(columns={'index':'MAG_ids'})
final_results = os.path.join(outdir, 'Final_taxa_assignment.txt')
taxa_change_log_df.to_csv(final_results, sep='\t', header=True, index=False, mode='w', na_rep='NA')
print(f'{time_current()} Final results were wirtten to file: {final_results}.')
