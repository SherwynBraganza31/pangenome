import sqlite3
import pandas as pd
import os


def condense_to_db(dbname: str, path: str, trimmed_matrix: pd.DataFrame, gene_ont: pd.DataFrame):
    if os.path.exists(f"{path}{dbname}"):
        os.remove(f"{path}{dbname}")

    cnx = sqlite3.connect(":memory")

    genes_table = trimmed_matrix.loc[:, ['Gene ID', 'Gene Family']].set_index('Gene ID')
    genes_table.to_sql('genes', cnx, if_exists='fail')

    organism_link = trimmed_matrix.loc[
        trimmed_matrix['Gene ID'] == trimmed_matrix['Gene Family'], ['Organism', 'Genome ID']]
    organism_link = organism_link.drop_duplicates(subset=['Organism']).set_index('Organism')
    organism_link.to_sql('organisms', cnx, if_exists='fail')

    gene_families_table = trimmed_matrix.drop_duplicates(subset=['Gene Family']).set_index('Gene Family').loc[:,
                          ['Annotation', 'Partition', '# Isolates', '# Sequences', 'Avg seqs/isolate', 'UniProtKB']]
    gene_families_table.to_sql('gene_families', cnx, if_exists='fail')

    go_processes = gene_ont.set_index('Gene Ontology ID').loc[:,
                   ['Gene Ontology Classification', 'Classification Value', 'UniProtKB']]
    go_processes.to_sql('go_processes', cnx, if_exists='fail')

    uniprots = gene_ont.drop_duplicates(subset=['UniProtKB']).set_index('UniProtKB').loc[:,
               ['Protein Name', 'Count']]
    uniprots.to_sql('uniprot_proteins', cnx, if_exists='fail')

    cnx.execute("vacuum main into " + f"\'{path}{dbname}\'")
    cnx.close()


if __name__ == '__main__':
    parent_dir = input("Enter postprocessed data directory: ")
    dbname = input("Enter name of the database you want to store it as: ")

    parent_dir = parent_dir if parent_dir[-1] == '/' else parent_dir + '/'
    dbname = dbname if ".db" in dbname else dbname + ".db"
    
    try:
        trimmed_matrix = pd.read_csv(parent_dir + 'genes.csv')
    except FileNotFoundError:
        print("Please run matrix_trimmer.py first.")
        quit()

    try:
        gene_ont = pd.read_csv(parent_dir + 'go_ids.csv')
    except FileNotFoundError:
        print("Please run gene_ont_proc.py before running this")
        quit()

    condense_to_db(dbname, parent_dir, trimmed_matrix, gene_ont)
    