from caseolap._08_count_synonyms import *
import sys, json, time, os
from elasticsearch import Elasticsearch
from elasticsearch_dsl import Search, Q
from multiprocessing import cpu_count, Process

'''
Parameters
'''
### Input
#entity_dict_path = 'data/filtered_entities100ENTS.txt'
entity_dict_path = 'data/casesensitive_entities.txt'  
textcube_pmid2category = 'data/textcube_pmid2category.json'

### Intermediary file (produced as output, used as input)
syn_pmid_count = 'data/syn_pmid_count.txt'

### Output 
pmid_syn_count_out = 'data/pmid_synonym_counts.json'   # PMID Syn|Count...Syn|Cnt
synfound_pmid2cat = 'data/synfound_pmid2category.txt'  # PMID--->CategoryNumber
logfile = 'log/synonymcount_log.txt'                   # #hits:Synonym

### Other parameters
index_name = 'pubmed' # Index name
key = 'abstract'      # Choose if searching the abstracts and titles
#key = 'full_text'    # Choose if searchines the abstracts, titles, and full text


'''
Main code
'''
if __name__ == '__main__':
    
    # Instantiate the object
    CS = CountSynonyms(entity_dict_path, textcube_pmid2category)
    
    # Search for the synonyms in the indexed text
    CS.synonym_search(key, logfile, syn_pmid_count, index_name) 
   
    # Finalize the output files
    CS.finish_synonym_search(logfile, syn_pmid_count,\
                             pmid_syn_count_out, synfound_pmid2cat)
