'''
The purpose of this file is to download the zipped files containing
the PubMed publications (i.e. documents). These will be mined later.
'''
import os, sys, json
from caseolap._01_download import *



'''
Parameters
'''
### Input
data_dir = './'
download_config_file_path = './config/download_config.json'# What to download
ftp_config_file_path = './config/ftp_config.json'          # Where to download from
baseline_files = './ftp.ncbi.nlm.nih.gov/pubmed/baseline/' # Documents from prior years
baseline_files = os.path.join(data_dir, baseline_files)    # " "
update_files = 'ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/'  # Documents from this year 
update_files = os.path.join(data_dir, update_files)        # " "
logfile_path = './log/download_log.txt'                    # Log download progress




'''
Main Code
'''
if __name__ == '__main__':
    
    ### Open files
    logfile = open(logfile_path, "w") 
    download_config = json.load(open(download_config_file_path, 'r'))
    ftp_config = json.load(open(ftp_config_file_path, 'r'))
    
    ### Check main directory
    if not os.path.isdir(data_dir):
        print("Directory not found:", data_dir) 
    
    ### Start download
    download_pubmed(data_dir, download_config, ftp_config, logfile)
    
    ### Verify download
    check_all_md5_in_dir(baseline_files, logfile)
    check_all_md5_in_dir(update_files, logfile)

    ### Extract downloaded files
    extract_all_gz_in_dir(baseline_files, logfile)
    extract_all_gz_in_dir(update_files, logfile)
    
    ### Close log file
    logfile.close()
    
    