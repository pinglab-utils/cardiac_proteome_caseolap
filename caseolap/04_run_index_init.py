'''
The purpose of this file is to initialize the ElasticSearch index.
'''
import json
from elasticsearch import Elasticsearch



'''
Parameters
'''
### Index parameters
index_name = "pubmed"      # Index name (match with 05 file)
type_name = "pubmed_meta"  # Index type name (match with 05 file)
number_shards = 1          # Set to 1 if no cluster
number_replicas = 0    
case_sensitive = True      # Index the text as case sensitive (True) or lower case (False)

### Input file
index_init_config_file = './config/index_init_config.json'



'''
Main Code
'''
if __name__ == '__main__':
    ### Load the indexing config file
    index_init_config = json.load(open(index_init_config_file,'r')) 
    
    ### Create an index request body
    request_body = {
        "settings": {
            "number_of_shards": number_shards,
            "number_of_replicas": number_replicas},
        "mappings": {
            type_name: {"properties": index_init_config}}}
    
    ### Indicate if the indexing will be case sensitive
    if case_sensitive == True:
        request_body["settings"]["analysis"] = {
            "analyzer":{
                "casesensitive_text": {
                    "type":"custom",
                    "tokenizer":"standard",
                    "filter": ["stop"]}}}
       
    
    
    ### Start elasticsearch
    es = Elasticsearch()
              
    ### Delete the old index if it exists
    if es.indices.exists(index_name):
        res = es.indices.delete(index = index_name)
        print("Deleting index %s , Response: %s" % (index_name, res))

    ### Create a new index   
    res = es.indices.create(index = index_name, body = request_body)
    print("Successfully crested index %s , Response: %s" % (index_name, res))
    
    