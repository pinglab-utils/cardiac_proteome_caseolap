'''
The purpose of this file is to create case-sensitive variations
of the synonyms. The synonyms will then be queried.
'''
from caseolap._07_vary_synonyms_cases import *


'''
Paths
'''
### Input
entity_dict_path = 'input/entities.txt'      # 1st entity dict: Entity_ID|syn1|...|synN
species = ['human', 'pig', 'mouse', 'rat']  # Species studied. Used for permuting syns.

### Output
case_entities_path = 'data/casesensitive_entities.txt' # Case sensitive entity dict



'''
Main code
'''
if __name__ == '__main__':
    # Instantiates the class
    VSC = VarySynonymsCases()
    
    # Loads entity dictionary mapping ID to synonyms
    VSC.load_id2syns_dict(entity_dict_path)
    
    # Adds some more synonyms
    VSC.add_species_syns(species)

    # Makes case-sensitive variations of the synonyms
    VSC.gets_final_syns(case_entities_path)