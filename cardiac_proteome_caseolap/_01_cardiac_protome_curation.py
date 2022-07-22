import os
import sys
import re
import time # is this necessary?
import subprocess # is this necessary?
import random # is this necessary?
from itertools import combinations
import collections
import numpy as np
import pandas as pd
import urllib.parse
import urllib.request
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_theme()

 
'''
Extract proteomes =============================
'''
# Assume fasta format
def get_proteins_from_fasta(file_name):
    proteins = set()
    reviewed = set()
    unreviewed = set()
    for l in [l.strip("\n") for l in open(file_name,"r").readlines()]:
        if ">" in l:
            accession = l.split("|")[1]
            proteins.add(accession)
            if "tr" in l:
                unreviewed.add(accession)
            elif "sp" in l:
                reviewed.add(accession)
    print("%d proteins (%d reviwed, %d unreviewed) in %s" %(len(proteins),len(reviewed),len(unreviewed),file_name))
    return proteins,reviewed,unreviewed


def merge_lists(input_directory = "./protein_lists", 
                out_file = "./merged_list.txt",
                remove_isoforms = False,
                filter_with_uniprot_db_file = None,
               debug=False):
    '''
    This function reads all the .txt files in the input directory and merges the unique identifiers into merged_list.txt
    remove_isoforms replaces protein-isoforms with its base protein (e.g. ProteinA-2 -> ProteinA)
    filter_with_uniprot_db_file retains only protein identifiers found within the reference database 
        (e.g. filter against Human ref proteome)
    Return values: unique_proteins, filtered_id_to_proteins
        unique_proteins - a set of unique protein identifiers among all lists
        filtered_id_to_proteins - a mapping from file name to proteins found within
    '''
    input_files = []
    # get all file names
    for filename in os.listdir(input_directory):
        if filename.endswith(".txt"):
            input_files += [os.path.join(input_directory, filename)]
        else:
            continue
    if debug:
        print("%d files in directory %s" %(len(input_files),input_directory))

    # load uniprot database if filter_with_uniprot_db_file is specified
    if filter_with_uniprot_db_file != None:
        uniprot_all, uniprot_reviewed, uniprot_unreviewed = get_proteins_from_fasta(filter_with_uniprot_db_file)

    # read all files and extract proteins for each
    id_to_proteins = {}
    total_prot_count = 0
    for f in input_files:
        id_ = f.split("/")[-1].strip(".txt")
        proteins = set(filter(None, [l.strip("\n") for l in open(f,"r").readlines()]))
        id_to_proteins[id_] = proteins
        total_prot_count += len(proteins)
    if debug:
        print("Total number of proteins among all files: %s"%total_prot_count)

    # get unique protein list
    unique_proteins = set()
    for id_ in id_to_proteins.keys():
        proteins = id_to_proteins[id_]
        unique_proteins = unique_proteins.union(proteins)
    if debug:
        print("%d unique proteins identified" %(len(unique_proteins)))

    # remove isoforms
    if remove_isoforms:
        prots = list(unique_proteins)
        unique_proteins = set()
        num_isoforms = 0
        for p in prots:
            if "-" not in p:
                unique_proteins.add(p)
            else:
                unique_proteins.add(p.split("-")[0])
                num_isoforms += 1
        if debug:
            print("Number of isoforms replaced: %d" %num_isoforms)
            print("%d unique proteins identified after removing isoforms" %(len(unique_proteins)))

        
    # cross reference against uniprot database
    if (filter_with_uniprot_db_file != None):
        prots = list(unique_proteins)
        unique_proteins = set()

        # only keep proteins in the uniprot db
        for p in prots:
            if p in uniprot_all:
                unique_proteins.add(p)

        # separated by experiment, only keep proteins in uniprot db
        filtered_id_to_proteins = {}
        for id in id_to_proteins.keys():
            proteins = id_to_proteins[id]
            filtered_proteins = proteins.intersection(unique_proteins)
            if len(filtered_proteins) > 0:
                filtered_id_to_proteins[id] = filtered_proteins
        if debug:
            print("After filtering, %d experiments have relevant proteins" %(len(filtered_id_to_proteins)))
            print("Total number of proteins among all files, filtering with Uniprot: %s"%len(unique_proteins))

    # write file
    with open(out_file,"w") as out_f:
        out_f.write("\n".join(list(unique_proteins)))
        print("Written to file %s"%out_file)
        
    return unique_proteins, filtered_id_to_proteins


def extract_proteomes(species_to_folder_and_db, 
                      remove_isoforms=True, 
                      out_folder = "./protein_lists", 
                      debug=True):
    '''
    This function reads the protein lists from four species and merges into one list TODO
    '''
    # extract proteomes
    species_to_unique_proteins_and_id_to_proteins = {}
    for species, f_db in species_to_folder_and_db.items():
        input_dir, db_file = f_db
        out_f = out_folder + "/merged_%s_list.txt"%(species)
        unique_proteins, id_to_proteins = merge_lists(remove_isoforms=remove_isoforms, 
                                                          out_file=out_f,
                                                          input_directory = input_dir,
                                                          filter_with_uniprot_db_file = db_file)
        species_to_unique_proteins_and_id_to_proteins[species] = (unique_proteins,id_to_proteins)
        
    # merge all
    unique_proteins = set()
    id_to_proteins = dict()
    for species, up_idp in species_to_unique_proteins_and_id_to_proteins.items():
        up,idp = up_idp
        unique_proteins = unique_proteins.union(up)
        id_to_proteins = id_to_proteins | idp
        
        if debug:
            print("%s (%d):"%(species,len(idp)))
            for pxid, proteins in idp.items():
                print(pxid,len(proteins))
    species_to_unique_proteins_and_id_to_proteins['Combined'] = (unique_proteins,id_to_proteins)
    
    # write to file
    with open("./merged_lists.txt","w") as out_file:
        out_file.write("\n".join(unique_proteins))

    if debug:
        print("%d unique proteins"%(len(unique_proteins)))
        unique_experiments = set()
        for exp in id_to_proteins.keys():
            unique_experiments.add(exp.split("_")[0])
        print("%d experimental datasets"%(len(unique_experiments)))
    
    return species_to_unique_proteins_and_id_to_proteins


'''
Gene to proteins =============================
'''
def read_prot_to_gene_map(file_name):
    prot_to_gene = {}
    gene_set = set()
    lines =    [l.strip("\n") for l in open(file_name,"r").readlines()]
    for l in lines[1:]:
        prot,gene = l.split("\t")
        prot_to_gene[prot] = gene
        gene_set.add(gene)
    print("%d proteins mapped to %d genes" % (len(prot_to_gene),len(gene_set)))
    return prot_to_gene, gene_set


def separate_gene_mapping_by_species(species_to_unique_proteins_and_id_to_proteins, prot_to_gene, debug=False):
    species_to_gene_mapping_and_gene_set = {}
    
    for species, up_idp in species_to_unique_proteins_and_id_to_proteins.items():
        unique_proteins,_ = up_idp
        protein_to_gene = {k:v for k,v in prot_to_gene.items() if k in unique_proteins}
        gene_set = set(protein_to_gene.values())
        protein_set = set(protein_to_gene.items())
        species_to_gene_mapping_and_gene_set[species] = protein_to_gene, gene_set
        if debug:
            print("%s: %d genes %d proteins (out of %d; difference: %d)"%(species,
                                    len(gene_set), 
                                    len(protein_set),
                                    len(unique_proteins),
                                    len(unique_proteins)-len(set(protein_set))))
    return species_to_gene_mapping_and_gene_set

       
# need to reverse mapping for the histogram
def gene_to_proteins(protein_to_gene):
    gene_to_protein = {}
    for prot in protein_to_gene.keys():
        gene = protein_to_gene[prot]
        if gene not in gene_to_protein:
            gene_to_protein[gene] = set()
        gene_to_protein[gene].add(prot)
    return gene_to_protein


def get_combinations(my_list):
    '''
    This function gets every combination (in order) from my_list
    my_list is size n, this function returns combinations of sizes [2 ... n]
    '''
    combs = []
    for size in (range(2,len(my_list)+1)):
        combs += list(combinations(my_list,size))
    return combs

def overlap_gene_set(species_to_genes):
    '''
    This function calculates the data needed for the staircase plot.
    Reports:
            'All': the number of unique genes among all species
            species: the number of genes found in the species
            combinations: the number of genes found in between combination of species
    '''
    return_dict = {}

    # calculate overlap of all
    union_set = set()
    for genes in species_to_genes.values():
        union_set = union_set.union(genes)
    return_dict['All'] = union_set

    # next, input the individual lists to keep the ordering correct
    for species, genes in zip(species_to_genes.keys(), species_to_genes.values()):
        return_dict[species] = genes

    # finally, their combinations
    combs = get_combinations(species_to_genes.keys())
    for species_list in combs:
        comb_set = union_set
        for s in species_list:
            comb_set = comb_set.intersection(species_to_genes[s])
        return_dict[species_list] = comb_set
    
    return return_dict


def staircase_plot_horiz(overlap_sets):
    '''
    This function visualizes the overlap_sets from previous cell
    '''

    # extract x and y data for bar chart
    x = []
    y = []
    for label, genes in zip(overlap_sets.keys(),overlap_sets.values()):
        if type(label) == tuple:
            label_italicized = ['$\it{%s}$ $\it{%s}$'%(l.split(" ")[0], l.split(" ")[1]) for l in label] # italicize each label
            x_ = " ∩ ".join(label_italicized)
        elif label != 'All':
            # if you don't split by spaces, it skips the space...
            label_ital = '$\it{%s}$ $\it{%s}$'%(label.split(" ")[0], label.split(" ")[1])
            x_ = label_ital # italicized
        else:
            x_ = label
        x += [x_]
        y += [len(genes)]
    # reverse values
    x = x[::-1]
    y = y[::-1]

    # labels formatted with comma
    y_labels = [f"{y_:,}" for y_ in y]
    
    # plot
    fig, ax = plt.subplots()
    plt.barh(x, y, log=False) 
    plt.xlabel('Number of genes') 
    
    # Text on the top of each bar
    for i in range(len(x)):
            plt.text(y = i-0.2 , x = y[i], s = y_labels[i], size = 8)
    
    # remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    
def ital_species(s):
    ''' italicizes species name, need to skip the space'''
    return '$\it{%s}$ $\it{%s}$'%(s.split(" ")[0], s.split(" ")[1])

def staircase_string_formatter(overlap_sets):
    ''' converts overlap_sets tuples to string to display on plot'''
    
    relabeled_overlap_sets = dict()
    for label, genes in zip(overlap_sets.keys(),overlap_sets.values()):
        if type(label) == tuple:
            label_italicized = [ital_species(s) for s in label] # italicize each label
            new_label = " ∩ ".join(label_italicized)
        elif label != 'All':
            new_label = ital_species(label) # italicized
        else:
            # don't italicize 'All'
            new_label = label
        relabeled_overlap_sets[new_label] = genes
    return relabeled_overlap_sets

def staircase_plot_horiz_1(new_overlap_sets):
    '''
    This function visualizes the overlap_sets from previous cell
    '''

    # extract x and y data for bar chart
    x = []
    y = []
    for label, genes in zip(new_overlap_sets.keys(),new_overlap_sets.values()):
        x += [label]
        y += [len(genes)]

    # reverse values
    x = x[::-1]
    y = y[::-1]

    # labels formatted with comma
    y_labels = [f"{y_:,}" for y_ in y]
    
    # plot
    fig, ax = plt.subplots()
    plt.barh(x, y, log=False) 
    plt.xlabel('Number of genes') 
    
    # Text on the top of each bar
    for i in range(len(x)):
            plt.text(y = i-0.2 , x = y[i], s = y_labels[i], size = 8)
    
    # remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    
    # background white
    ax.set_facecolor('white')
    
    

def gene_histogram(gene_to_protein_set, tag="", out_file = "./figures/gene_histogram.pdf"):
    '''
    This function takes a gene to protein set mapping and makes a histogram
    showing the number of proteins that correspond to each gene. 
    '''
    
    # Make plot area and background to see on notebook better
    fig, ax = plt.subplots(figsize = [10,4])
    fig.patch.set_facecolor('white')

    # collect frequencies
    freqs = []
    for gene in gene_to_protein_set:
        freqs += [len(gene_to_protein_set[gene])]
    x = np.array(range(max(freqs)))+1
    
    # gather histogram data and plot
    coll = collections.Counter(freqs)
    y = [coll[x_] for x_ in x]
    plt.bar(x, y, log=True)

    # Text on the top of each bar
    for i in range(len(x)):
        if y[i] > 0:
            plt.text(x = x[i] , y = y[i]+0.5, s = str(y[i]), size = 12, va='bottom', ha = 'center')
            
    # formatting 
    plt.xlabel('Number proteins isoforms per gene') 
    plt.ylabel('Frequency') 
    plt.title("%s Cardiac Protein-Gene Mapping"%tag)
    ax.tick_params(axis='both', labelsize=14)
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # set x-axis
    if max(x) <= 10:
        plt.xlim([0,10+1])
    else:
        plt.xlim([0,max(x)+1])
    
    # # save file if file name is given
    if len(out_file) > 0:
            fig.savefig(out_file,bbox_inches="tight")
    plt.show()
    
    
def protein_experiment_histogram(unique_proteins, id_to_proteins, tag="", out_file = "./figures/dataset_histogram.pdf"):
    '''
    This function plots a histogram of number of proteins found within each experiment. 
    '''
    
    # Make plot area and background to see on notebook better
    fig, ax = plt.subplots(figsize = [5,3.5])
    fig.patch.set_facecolor('white')

    # frequency for each protein
    protein_to_datasets = {}
    freqs = []
    for p in unique_proteins:
        datasets = set()
        for dset in id_to_proteins.keys():
            if p in id_to_proteins[dset]:
                datasets.add(dset)
        freqs += [len(datasets)]
    
    x = np.array(range(max(freqs)))+1
    coll = collections.Counter(freqs)
    y = [coll[x_] for x_ in x]
    plt.bar(x, y, log=True)

   # Text on the top of each bar
    for i in range(len(x)):
        if y[i] > 0:
            plt.text(x = x[i] , y = y[i]+0.5, s = str(y[i]), size = 12, va='bottom', ha = 'center')
            
    # formatting 
    plt.xlabel('Number of %s Cardiac Datasets'%(tag)) 
    plt.ylabel('Protein Observation Frequency') 
    plt.title("%s Cardiac Proteins"%tag)
    ax.tick_params(axis='both', labelsize=14)
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # set x-axis
#     plt.xlim([0,12.5])
    
    # # save file if file name is given
    if len(out_file) > 0:
            fig.savefig(out_file,bbox_inches="tight")
    plt.show()
    
    
def protein_experiment_heatmap(id_to_proteins, tag="", out_file = "./figures/dataset_histogram.pdf"):
    '''
    This function reads all the .txt files in the input directory and merges the unique identifiers into merged_list.txt
    '''

    fig = plt.gcf()
    fig.set_size_inches(10, 7)
    
    # get overlap matrix
    df = pd.DataFrame(columns=id_to_proteins.keys())
    for id1 in id_to_proteins.keys():
        
        row = {}
        for id2 in id_to_proteins.keys():
            p1 = id_to_proteins[id1]
            p2 = id_to_proteins[id2]
            overlap = p1.intersection(p2)
            row[id2] = len(overlap)
        df = df.append(row,ignore_index=True)
    
    # replace rows and column labels
    new_headers = list(id_to_proteins.keys())
    # rename one label
    if 'mayer_pubs' in new_headers: # special case for this experiment which has multiple identifiers
        new_headers[new_headers.index("mayer_pubs")] = "PXD003918, PXD003919, PXD003926"
    new_headers = [h.split("_")[0] for h in new_headers]
    # new_headers[new_headers.index("mayer_pubs")] = "Mayr Lab"
    df.columns = new_headers
    df.index = new_headers
    df = df.astype('int32')
    
    # grey out the lower right quadrant
    mask = []
    for i in range(len(id_to_proteins)):
        m = []
        for j in range(len(id_to_proteins)):
            m += [(j <= i)]
        mask += [m]
    mask = pd.DataFrame(mask)
    mask.columns=new_headers
    mask.index=new_headers

    ax = sns.heatmap(df, annot=True, fmt="d", mask=mask, 
                                     vmin=0, cmap='YlGnBu',
                                    cbar_kws={'label':'Color Code'})
    
    # set title text
    ax.set_title("%s Cardiac Datasets"%(tag))
    
    # rotate x axis labels
    ax.set_xticklabels(
        ax.get_xticklabels(), 
        rotation=45, 
        horizontalalignment='right'
    )
    
    # set tick marks
    ax.tick_params(left=True,bottom=True)
        
    # # save file if file name is given
    if len(out_file) > 0:
            fig.savefig(out_file,bbox_inches="tight")
    plt.show()
    
    
'''
Remove mutations based on sequence similarity =============================
'''
def parse_fasta(file_name):
    prot_to_seq = {}
    lines = [l.strip("\n") for l in open(file_name,"r").readlines()]
    accession = ""
    sequence = ""
    for l in lines:
        if ">" in l:
            if len(sequence) > 0:
                # keep track of previous accession and sequence extraction
                prot_to_seq[accession] = sequence
            accession = l.split("|")[1]
            sequence = ""
        else:
            sequence += l
    prot_to_seq[accession] = sequence
    print("%d proteins parsed in %s" %(len(prot_to_seq),file_name))
    return prot_to_seq



def score_seqs(seq1,seq2):
    '''
    This function reports the number of mismatches and insertions between seq1 and seq2.
    '''
    alignment_score = pairwise2.align.globalxx(seq1, seq2, score_only=True)
    score = max(len(seq1),len(seq2))-alignment_score
    return score

     

def filter_protein_set(protein_list, prot_to_seq,
                                        reviewed_list = None, filter_threshold = 3, debug = False):
    '''
    This function filters the gene_to_protein_set by aligning pairs of sequences
    and removing arbitrarily proteins that are too similar. Two proteins are too 
    similar if they differ in fewer than or equal to filter_threshold amino acids.
    Proteins that are in the reviewed_list (if included) are retained in priority.
    Pairs of reviewed proteins and pairs of unreviewed proteins are chosen 
    arbitrarily.
    '''

    # don't remove anything if only one protein
    if len(protein_list) <= 1 :
        return protein_list

    # get every pairwise combination from the list
    pairs = list(combinations(protein_list, 2))
    
    # keep track of their scores
    scores = {}
    for s1, s2 in pairs:
        seq1 = prot_to_seq[s1]
        seq2 = prot_to_seq[s2]
        scores[(s1,s2)] = score_seqs(seq1,seq2)
    # if(debug):
    #     print(scores)

    # remove proteins based on scores 
    to_keep = set(protein_list)
    for s1, s2 in pairs:
        score = scores[(s1,s2)]

        # remove one if score less than equal to threshold
        if score <= filter_threshold:

            keep_s1 = True
            if reviewed_list != None:
                # break ties based on reviewed list if helpful
                if s1 in reviewed_list and s2 not in reviewed_list:
                    keep_s1 = True
                elif s2 in reviewed_list and s1 not in reviewed_list:
                    keep_s1 = False
                else:
                    keep_s1 = np.random.rand() > .5    # 50-50 chance
            else:
                keep_s1 = np.random.rand() > .5    # 50-50 chance
            
            if(debug):
                print("%s and %s have a score of %s"%(s1,s2,score))
            if keep_s1:
                if(debug):
                    print("Keeping %s"%(s1))
                to_keep.add(s1)
                to_remove = s2
            else:
                if(debug):
                    print("Keeping %s"%(s2))
                to_keep.add(s2)
                to_remove = s1
            if to_remove in to_keep:
                if(debug):
                    print("Removing %s"%(to_remove))
                to_keep.remove(to_remove)

    if(debug and len(to_keep) < len(protein_list)):
        print("Removed %d out of %d sequences."%(len(protein_list)-len(to_keep),len(protein_list)))
        print("%d remaining." %(len(to_keep)))
    return to_keep


def filter_gene_to_protein_set(gene_to_protein_set, prot_to_seq,
                                        reviewed_list = None, filter_threshold = 3, debug=False):
    '''
    This function filters the gene_to_protein_set by aligning pairs of sequences
    and removing arbitrarily proteins that are too similar. Two proteins are too 
    similar if they differ in fewer than or equal to filter_threshold amino acids.
    Proteins that are in the reviewed_list (if included) are retained in priority.
    Pairs of reviewed proteins and pairs of unreviewed proteins are chosen 
    arbitrarily.
    '''

    filtered_gene_to_protein_set = {}
    for gene in gene_to_protein_set.keys():
        proteins = gene_to_protein_set[gene]
        proteins_to_keep = filter_protein_set(proteins, prot_to_seq,
                                                            reviewed_list = reviewed_list, 
                                                            filter_threshold = filter_threshold, debug=debug)
        filtered_gene_to_protein_set[gene] = proteins_to_keep
    return filtered_gene_to_protein_set


'''
I/O =============================
'''
# Output the proteins from the filtered results
def output_cardiac_proteome(gene_to_protein_set, out_file_name = "./cardiac_proteome.txt"):
    genes = set()
    proteins = set()

    for gene, protein_set in zip(gene_to_protein_set.keys(), gene_to_protein_set.values()):
        genes.add(gene)
        for p in protein_set:
            proteins.add(p)

    print("%d genes and %d proteins"%(len(genes),len(proteins)))
    
    with open(out_file_name,"w") as out_file:
        out_file.write("\n".join(list(proteins)))
        print("%d proteins written to %s" %(len(proteins),out_file_name))
        
    genes_out_file = out_file_name[:-4] + "_genes"+out_file_name[-4:]
    with open(genes_out_file,"w") as out_file:
        out_file.write("\n".join(list(genes)))
        print("%d genes written to %s" %(len(genes),genes_out_file))
        

def load_cardiac_proteome(species_to_gene_and_protein_lists, gene_to_protein_mapping_file):
    # read the gene mapping (protein -> gene)
    protein_to_gene, _ = read_prot_to_gene_map(proteins_to_gene_map)
                                               
    # parse the cardiac proteomes for each species
    species_to_cardiac_proteome = {}
    for species, gene_protein_list_files in species_to_gene_and_protein_lists.items():
        # read the gene list
        gene_list = set([l.strip("\n") for l in open(gene_protein_list_files[1],"r").readlines()])
        # read the protein list
        protein_list = set([l.strip("\n") for l in open(gene_protein_list_files[0],"r").readlines()])
        
        # keep only gene to protein mappings found within gene and protein lists
        cardiac_proteome = {}
        for protein in protein_list:
            if protein in protein_to_gene:
                gene = protein_to_gene[protein]
                if gene not in cardiac_proteome:
                    cardiac_proteome[gene] = set()
                cardiac_proteome[gene].add(protein)
                
        genes_parsed = cardiac_proteome.keys()
        proteins_parsed = set()
        for g,ps in cardiac_proteome.items():
             proteins_parsed = proteins_parsed.union(ps)
        print("%d (out of %d) genes and %d (out of %d) proteins parsed for %s"%(len(genes_parsed),
                                                        len(gene_list),
                                                        len(proteins_parsed),
                                                        len(protein_list),
                                                        species))
        species_to_cardiac_proteome[species] = cardiac_proteome
    return species_to_cardiac_proteome

    
    
def convert_string_ids_to_uniprot(query_list, out_file="./mappings/uniref_90_mapping_table.tsv"):
    '''
    This function maps STRING ID's to Uniprot accession.
    Returns a list of UniProt proteins and a mapping function
    from STRING ID's to Uniprot accession.
    '''

    # Use UniProt API to convert STRING ID's to UniProt
    url = 'https://www.uniprot.org/uploadlists/'

    params = {
    'from': 'ACC', # from UniProt
    'to': 'NF90', # to UniRef90
    'format': 'tab',
    'query': " ".join(query_list)
    }

    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
        response = f.read()
    
    # parse the results into a map
    uniprot_to_uniref = {}
    uniref_to_uniprot = {}
    raw_string = response.decode('utf-8')
    # split into table but leave out header and last line
    mapping_table = raw_string.split("\n")[1:-1]
    for l in mapping_table:
        uniprot_id,uniref_id = l.split("\t")
        # uniprot to uniref is a 1-to-1 mapping
        uniprot_to_uniref[uniprot_id] = uniref_id
        # uniref90 is a 1-to-many mapping
        if uniref_id not in uniref_to_uniprot:
            uniref_to_uniprot[uniref_id] = set()
        uniref_to_uniprot[uniref_id].add(uniprot_id)
    print("%d out of %d UniProt ID's mapped to %d UniRef90 ID's"%(len(uniprot_to_uniref.keys()),
                                                                 len(query_list),
                                                                 len(uniref_to_uniprot.keys())))
    if(out_file):
        with open(out_file,"w") as out:
            count = 1
            out.write("From\tTo")
            for uniprot, uniref in uniprot_to_uniref.items():
                out.write("\n"+"\t".join([uniprot,uniref]))
                count+=1
        print("Wrote %d lines to file %s"%(count,out_file))
    return uniprot_to_uniref, uniref_to_uniprot


def convert_uniprot_to_gene_names(query_list, out_file="./mappings/merged_protein_lists_to_genes.tsv"):
    '''
    This function maps UniProt accession to Gene names.
    Returns a list of UniProt proteins and a mapping function
    from STRING ID's to Uniprot accession.
    '''

    # Use UniProt API to convert STRING ID's to UniProt
    url = 'https://www.uniprot.org/uploadlists/'

    params = {
    'from': 'ACC', # from UniProt
    'to': 'gene', # to UniRef90
    'format': 'tab',
    'query': " ".join(query_list)
    }

    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
        response = f.read()
    
    # parse the results into a map
    uniprot_to_gene = {}
    gene_to_uniprot = {}
    raw_string = response.decode('utf-8')
    # split into table but leave out header and last line
    mapping_table = raw_string.split("\n")[1:-1]
    for l in mapping_table:
        uniprot_id,gene = l.split("\t")
        uniprot_to_gene[uniprot_id] = gene
        if gene not in gene_to_uniprot:
            gene_to_uniprot[gene] = set()
        gene_to_uniprot[gene].add(uniprot_id)
    print("%d out of %d UniProt ID's mapped to %d UniRef90 ID's"%(len(uniprot_to_gene.keys()),
                                                                 len(query_list),
                                                                 len(gene_to_uniprot.keys())))
    if(out_file):
        with open(out_file,"w") as out:
            count = 1
            out.write("From\tTo")
            for uniprot, uniref in uniprot_to_uniref.items():
                out.write("\n"+"\t".join([uniprot,uniref]))
                count+=1
        print("Wrote %d lines to file %s"%(count,out_file))
    return uniprot_to_gene, gene_to_uniprot