import json
import pandas as pd
import os
from tqdm.auto import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


class GeneOntology:
    def __init__(self, parent_dir: str, results_path: str = None):
        self.parent_dir = parent_dir if parent_dir[-1] == '/' else parent_dir + '/'
        self.results_path = results_path if results_path is not None else self.parent_dir + 'trimmed_matrix_files/'
        self.ecNumbers = None
        self.goIDs = None
        self.biological = None
        self.cellular = None
        self.molecular = None

    def process_go_ids(self):
        with open(self.results_path + 'go_results.json', 'r', encoding='utf-8') as ifile:
            results = json.load(ifile)

        with open(self.results_path + 'uniprot_freqs.json', 'r') as ifile:
            uniprot_freqs = json.load(ifile)

        df = pd.DataFrame(columns=['Gene Ontology ID', 'Protein Name', 'UniProtKB',
                                   'Gene Ontology Classification', 'Classification Value',
                                   'Count'])
        ecNumbers = []
        for entry in tqdm(results['results']):
            uniprot_id = entry['from']
            prot_name = entry['to']['proteinDescription']['recommendedName']['fullName']['value']
            try:
                temp_nums = entry['to']['proteinDescription']['recommendedName']['ecNumbers']
                for dictionary in temp_nums:
                    ecNumbers.append(dictionary['value'])
            except KeyError:
                pass
            for db_entry in entry['to']['uniProtKBCrossReferences']:
                if db_entry['database'] == 'GO':
                    go_id = db_entry['id']
                    extracted_val = db_entry['properties'][0]['value']
                    go_class_value = extracted_val[2:]
                    if extracted_val[0] == 'P':
                        go_class = 'Biological Process'
                    elif extracted_val[0] == 'C':
                        go_class = 'Cellular Function'
                    elif extracted_val[0] == 'F':
                        go_class = 'Molecular Function'
                    else:
                        go_class = 'Unknown Class'

                    try:
                        gene_count = int(uniprot_freqs[uniprot_id])
                    except KeyError:
                        print(f'Uniprot : {uniprot_id} for GO id : {go_id} does not exist?')
                        gene_count = None

                    df.loc[len(df.index)] = [go_id, prot_name, uniprot_id,
                                             go_class, go_class_value, gene_count]

        self.goIDs = df
        self.ecNumbers = ecNumbers
        self.goIDs.set_index("Gene Ontology ID", inplace=True)
        self.goIDs.to_csv(self.results_path + 'go_ids.csv')

        with open(self.results_path + 'ecNumbers.txt', "a") as f:
            for x in self.ecNumbers:
                f.write(x + '\n')
        return

    def go_class_splitter(self):
        # Read entries from the .csv that links GO ID(s) to geneID(s)
        # The aim is to group up all UniProt ID(s) with the same GO ID
        # and then count the number of Genes associated with a particular
        # GO subfunction for each of the GO Function Categories - Biological Processes,
        # Molecular Functions and Cellular Functions
        # We then want to calculate the scores for each of these subprocesses, defined
        # as - (# genes involved in a subprocess) / (total # of genes)

        if (self.goIDs is None) and ("go_ids.csv" not in os.listdir(self.results_path)):
            raise TypeError("No GO IDs present")
        elif self.goIDs is None:
            self.goIDs = pd.read_csv(self.results_path + "go_ids.csv")
        else:
            pass

        if self.goIDs.index.name is None:
            self.goIDs.set_index("Gene Ontology ID", inplace=True)
        self.biological = self.goIDs.loc[self.goIDs["Gene Ontology Classification"].isin(["Biological Process"])]
        self.cellular = self.goIDs.loc[self.goIDs["Gene Ontology Classification"].isin(["Cellular Function"])]
        self.molecular = self.goIDs.loc[self.goIDs["Gene Ontology Classification"].isin(["Molecular Function"])]

        self.calculateScores()

        self.biological.to_csv(self.results_path + 'biological.csv')
        self.cellular.to_csv(self.results_path + 'cellular.csv')
        self.molecular.to_csv(self.results_path + 'molecular.csv')

    def calculateScores(self):
        dataframes = [self.biological, self.cellular, self.molecular]
        for idx, dataframe in enumerate(dataframes):
            sum_genes = np.sum(dataframe.loc[:, "Count"])
            dataframe.loc[:, "Score (% of genes involved)"] = (
                    dataframe.loc[:, "Count"].to_numpy() * 100/sum_genes)
            dataframe.sort_values(by="Score (% of genes involved)", ascending=False, inplace=True)

    def plotGeneOntology(self, fig_save_loc=None):
        fig, axs = plt.subplots(3, 1, figsize=(8, 20), sharex=True)
        plt.subplots_adjust(hspace=.0)
        data_labels = ['Biological Process', 'Cellular Function', 'Molecular Function']
        color_list = ['#cbe2e9', '#f1ffc4', '#ffcaaf']

        fig_save_loc = self.parent_dir + 'gene_ontology_distribution.jpeg'

        dataframes = [self.biological, self.cellular, self.molecular]
        for i, data in enumerate(dataframes):
            sns.scatterplot(ax=axs[i], data=data[0:5], x="Score (% of genes involved)", y="Classification Value",
                            size="Count", legend=False, sizes=(20, 1000), color='#5A5A5A')
            axs[i].set_ylabel(data_labels[i])
            axs[i].set_facecolor(color_list[i])
            axs[i].grid(visible=True)

        plt.rcParams.update({'font.size': 10})
        plt.savefig(fig_save_loc, dpi=512, bbox_inches='tight', pad_inches=0.2)


if __name__ == '__main__':
    parent_dir = input('Enter parent directory: ')
    temp_obj = GeneOntology(parent_dir=parent_dir)
    temp_obj.process_go_ids()
    temp_obj.go_class_splitter()
    temp_obj.plotGeneOntology()

