from sqlalchemy import create_engine
import pandas as pd
import scripts.pyensembl_operations as po

def setup_notebook():
    """
    Set up the environment for Jupyter Notebooks.
    """
    # Show all columns in pandas
    pd.set_option('display.max_columns', None)

    # Create SQLite engine
    engine = create_engine('sqlite:///data/db/mirscribe.db')

    # Import PyEnsembl data for genome build 37
    g37 = po.import_pyensembl(37)

    # Load data from SQLite database
    transcripts = pd.read_sql_table("transcripts", engine)
    genes = pd.read_sql_table("genes", engine)
    mirnas = pd.read_sql_table("mirnas", engine)

    return engine, g37, transcripts, genes, mirnas
