from socket import gethostname

def get_pyensembl_cache_location():
    hostname = gethostname()
    if hostname == "nazo":
        return "/home/nazif/thesis/data"
    else:
        return "/truba/home/mtasbas/data"

def get_rnaduplex_location():
    hostname = gethostname()
    if hostname == "nazo":
        return "/usr/local/bin/RNAduplex"
    else:
        return "/truba/home/mtasbas/miniconda3/envs/venv/bin/RNAduplex"
    
PYENSEMBL_CACHE_DIR = get_pyensembl_cache_location()
RNADUPLEX_LOCATION = get_rnaduplex_location()

GRCH37_DIR = "data/fasta/grch37"
MIRNA_COORDS_DIR = "data/mirna_coordinates"
TA_SPS_CSV = "data/ta_sps/ta_sps.csv"
MIRNA_CSV = "data/mirna/mirna.csv"
XGB_MODEL = "misc/models/model_with_no_close_proximity.json"
NUCLEOTIDE_OFFSET = 30

AWK_SCRIPT_PATH = "scripts/rnaduplex_to_csv.awk"

QUANTILE_RANGE = 0.25


