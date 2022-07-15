'''
The purpose of this file is so that the user can see why each entity scored
highly; the user can see the ranked entities.
Get the ranked entities.
Get the ranked synonyms.
'''
from caseolap._13_inspect_entity_scores import *



'''
Parameters
'''
### Input path 
#entities_path = 'data/filtered_entities100ENTS.txt'
id2syns_in = 'data/id2syns.json'                     # The case-varied entity dict 
scores_in = 'result/caseolap.csv'                    # The CaseOLAP scores
pmid_syn_count_in = 'data/pmid_synonym_counts.json'  # Counts of each synonym
remove_syns_in = 'data/remove_these_synonyms.txt'    # Syns that were not used
cat2pmids_in = 'data/metadata_category2pmids.json'   # Category->[PMID,...,PMID] 

### Output paths
ranked_syns_out = 'result/ranked synonyms/ranked_synonyms.txt' # Syns ranked by counts
ranked_ent_out = 'result/ranked entities/ranked_entities.txt'  # Ents ranked by score


'''
Main Code
'''
if __name__ == '__main__':
    # Initialize the class
    IES = InspectEntityScores(scores_in, id2syns_in, remove_syns_in,\
                              pmid_syn_count_in, cat2pmids_in)

    ### SAVE TO FILE: the ranked synonyms
    IES.get_ranked_synonyms_found(cat2pmids_in, pmid_syn_count_in, ranked_syns_out)
    
    ### Get the ranked entities
    # SORT: entities by CaseOLAP scores
    IES.sort_all_scores(scores_in)
    
    # DISPLAY: proportion entities found / entities searched for
    IES.prop_entities_found()
    
    # SAVE TO DICT: ranked entites, their synonyms, and their synonym counts
    IES.get_entity_syn_counts()
    
    # SAVE TO FILE & DISPLAY: ranked entities, their synonyms, and synonym counts
    IES.rank_each_category_score(ranked_ent_out)
    
    