from caseolap._10_make_entity_counts import *
'''
The purpose of this class is to make the entitycount file that maps 
PMIDs to entities to entity counts. This uses a "PMIDs to synonyms to synonym counts"
mapping and an "entity to synonym" mapping to do this. No ElasticSearch querying
is needed in this step.
'''


'''
Parameters
'''
### Input
remove_syns_infile = 'data/remove_these_synonyms.txt'
id2syns_path = 'data/id2syns.json'
pmid_syn_count_in = 'data/pmid_synonym_counts.json'

### Output
entitycount_outfile = 'data/entitycount.txt'

'''
Main code
'''
if __name__ == '__main__':

    # Instantiate and initialize class
    MEC = MakeEntityCounts(id2syns_path, remove_syns_infile, pmid_syn_count_in)
    
    # Makes "pmid->entity->count" from "pmid->syn->count" & "entity->syn" mappings
    MEC.entitycount(entitycount_outfile)