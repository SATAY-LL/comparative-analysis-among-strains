# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Python 3.8.10 ('satay-dev')
#     language: python
#     name: python3
# ---

# +
from goatools import obo_parser
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.test_data.genes_NCBI_10090_ProteinCoding import GENEID2NT as geneid2nt

import Bio.UniProt.GOA as GOA

import pandas as pd
import numpy as np



# +
# Load the Gene Ontology (GO) database
go_file = "../postprocessed-data/go-basic.obo"  # Path to the Gene Ontology OBO file
go = obo_parser.GODag(go_file)




# +
import gzip

## download from https://ftp.ebi.ac.uk/pub/databases/GO/goa/YEAST/goa_yeast.gaf.gz
# File is a gunzip file, so we need to open it in this way
with gzip.open("../postprocessed-data/goa_yeast.gaf.gz", 'rt') as yeast_gaf_fp:
    yeast_funcs = {}  # Initialise the dictionary of functions
    
    # Iterate on each function using Bio.UniProt.GOA library.
    for entry in GOA.gafiterator(yeast_gaf_fp):
        uniprot_id = entry.pop('DB_Object_ID')
        yeast_funcs[uniprot_id] = entry
# -

yeast_funcs_pd=pd.DataFrame.from_dict(yeast_funcs,orient='index')

# +
## Import the interactors of BEM1 
pop = yeast_funcs.keys()
assoc = {}
bem1PI=np.loadtxt("../postprocessed-data/positive_satay_genes_bem1.txt",dtype=str)
bem1NI=np.loadtxt("../postprocessed-data/negative_satay_genes_bem1.txt",dtype=str)
bem3PI=np.loadtxt("../postprocessed-data/positive_satay_genes_bem3.txt",dtype=str)
bem3NI=np.loadtxt("../postprocessed-data/negative_satay_genes_bem3.txt",dtype=str)
nrp1PI=np.loadtxt("../postprocessed-data/positive_satay_genes_nrp1.txt",dtype=str)
nrp1NI=np.loadtxt("../postprocessed-data/negative_satay_genes_nrp1.txt",dtype=str)
bem1bem3PI=np.loadtxt("../postprocessed-data/positive_satay_genes_bem1_bem3.txt",dtype=str)
bem1bem3NI=np.loadtxt("../postprocessed-data/negative_satay_genes_bem1_bem3.txt",dtype=str)

for x in yeast_funcs:
    if x not in assoc:
        assoc[x] = set()
    assoc[x].add(str(yeast_funcs[x]['GO_ID']))



# -

methods = ["bonferroni", "sidak", "holm", "fdr"]

# +

g = GOEnrichmentStudy(pop, assoc, go,
                         propagate_counts=True,
                         alpha=0.9,
                         methods=[methods[0]])



# +
study=yeast_funcs_pd[yeast_funcs_pd['DB_Object_Symbol'].isin(bem1bem3PI)].index ## enrichment of Bem1  positive interactors

g_res = g.run_study(study)

# +
s_fdr_term = []
s_fdr_p = []
for x in g_res:
    s_fdr_term.append(x.goterm.id)
    s_fdr_p.append(x.p_bonferroni)

aspect=["C"]
# -

yeast_funcs_pd.query('DB_Object_Symbol in @bem1bem3PI').query('GO_ID in @s_fdr_term').query('Aspect in @aspect').sort_values(by='DB_Object_Symbol')

yeast_funcs_pd[yeast_funcs_pd["GO_ID"]==s_fdr_term[2]]

from genes_ncbi_yeast_proteincoding import GENEID2NT as GeneID2nt_yeast

from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS

#run one time to initialize
obo_fname = download_go_basic_obo()
fin_gene2go = download_ncbi_associations()
obodag = GODag("go-basic.obo")

# +
#run one time to initialize
mapper = {}

for key in GeneID2nt_yeast:
    mapper[GeneID2nt_yeast[key].Symbol] = GeneID2nt_yeast[key].GeneID
    
inv_map = {v: k for k, v in mapper.items()}

# +
#run one time to initialize

# Read NCBI's gene2go. Store annotations in a list of namedtuples
objanno = Gene2GoReader(fin_gene2go, taxids=[559292])
# Get namespace2association where:
#    namespace is:
#        BP: biological_process               
#        MF: molecular_function
#        CC: cellular_component
#    assocation is a dict:
#        key: NCBI GeneID
#        value: A set of GO IDs associated with that gene
ns2assoc = objanno.get_ns2assc()
# -

#run one time to initialize
goeaobj = GOEnrichmentStudyNS(
        GeneID2nt_yeast.keys(), # List of mouse protein-coding genes
        ns2assoc, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.9, # default significance cut-off
        methods = ['holm']) # defult multipletest correction method

# +
#run one time to initialize
GO_items = []

temp = goeaobj.ns2objgoea['BP'].assoc
for item in temp:
    GO_items += temp[item]
    

temp = goeaobj.ns2objgoea['CC'].assoc
for item in temp:
    GO_items += temp[item]
    

temp = goeaobj.ns2objgoea['MF'].assoc
for item in temp:
    GO_items += temp[item]


# -

#pass list of gene symbols
def go_it(test_genes):
    print(f'input genes: {len(test_genes)}')
    
    mapped_genes = []
    for gene in test_genes:
        try:
            mapped_genes.append(mapper[gene])
        except:
            pass
    print(f'mapped genes: {len(mapped_genes)}')
    
    goea_results_all = goeaobj.run_study(mapped_genes)
    goea_results_sig = [r for r in goea_results_all if r.p_holm < 0.9]
    GO = pd.DataFrame(list(map(lambda x: [x.GO, x.goterm.name, x.goterm.namespace, x.p_uncorrected, x.p_holm,\
                   x.ratio_in_study[0], x.ratio_in_study[1], GO_items.count(x.GO), list(map(lambda y: inv_map[y], x.study_items)),\
                   ], goea_results_sig)), columns = ['GO', 'term', 'class', 'p', 'p_corr', 'n_genes',\
                                                    'n_study', 'n_go', 'study_genes'])

    GO = GO[GO.n_genes > 1]
    return GO


test=["CDC42","BEM1","BEM3","CDC24"]

df = go_it(test)

df

# +
## look here: https://github.com/mousepixels/sanbomics_scripts/blob/main/GO_in_python.ipynb
