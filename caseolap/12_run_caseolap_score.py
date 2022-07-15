'''
The purpose of this file is to produce CaseOLAP scores for the entities
based on their hits in each document (pmid2pcount_path) and the documents'
category (category2pmids_path).
'''
import pandas as pd, numpy as np, matplotlib.pyplot as plt, seaborn as sns, json
from caseolap._12_caseolap_score import *



'''
Parameters
'''
### Input data directories
cat2pmids_path = './data/metadata_category2pmids.json' # {CategoryName:[PMID,...],...}
pmid2pcount_path = './data/metadata_pmid2pcount.json'  # {PMID:{Entity:Count,...},...}

### Output data path
result_dir = "result/"
logFilePath = "./log/caseolap_score_log.txt" # Logs #PMIDs for each category



'''
Main Code
'''
if __name__ == '__main__':

    logfile = open(logFilePath, "w") 
    
    cell2pmids = json.load(open(cat2pmids_path, 'r'))   

    pmid2pcount = json.load(open(pmid2pcount_path, 'r'))


    ### Test Run
    C = Caseolap(cell2pmids,pmid2pcount,result_dir,logfile)
    C.cell_pmids_collector(dump =True,verbose =True)
    #C.cell_pmids
    C.cell_pmid2pcount_collector()
    #C.cell_pmid2pcount
    C.all_protein_finder(dump =True,verbose = True)
    #C.all_proteins
    C.all_protein_finder()
    #C.all_proteins
    #C.cell_uniqp   
    C.cell_p2tf_finder()
    #C.cell_p2tf
    C.cell_tf_finder()
    #C.cell_tf
    C.cell_pop_finder(dump=True)
    #C.cell_pop
    C.cell_p2pmid_finder()
    #C.cell_p2pmid
    C.cell_ntf_finder()
    #C.cell_ntf
    C.cell_ndf_finder()
    #C.cell_ndf
    C.cell_rel_finder()
    #C.cell_rel
    C.cell_dist_finder(dump=True)
    #C.cell_dist
    C.cell_cseolap_finder(dump=True)
    #C.cell_caseolap
    logfile.close()