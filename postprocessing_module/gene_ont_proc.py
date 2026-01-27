import json
import pandas as pd
import os
from tqdm.auto import tqdm
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns


class GeneOntology:
    def __init__(self, source_dir: str, results_path: str = None):
        self.source_dir = source_dir if source_dir[-1] == '/' else source_dir + '/'
        self.results_path = self.source_dir + 'postprocessing_results/'
        self.ecNumbers = None
        self.goIDs = None
        self.biological = None
        self.cellular = None
        self.molecular = None
        self.goFreqs = None

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
                ecNumbers.append(entry['to']['proteinDescription']['recommendedName']['ecNumbers'][0]['value'])
                # for dictionary in temp_nums:
                #     ecNumbers.append(dictionary['value'])
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

        go_uniques = pd.unique(self.goIDs.loc[:, 'Gene Ontology ID'])
        freq = np.zeros(go_uniques.shape[0])
        for idx, gid in enumerate(go_uniques):
            freq[idx] = self.goIDs.loc[self.goIDs.loc[:, 'Gene Ontology ID'].isin([gid])].shape[0]
        self.goFreqs = self.goIDs.drop_duplicates(['Gene Ontology ID'])
        self.goFreqs.loc[:, 'Count'] = freq
        self.goFreqs.set_index("Gene Ontology ID", inplace=True)
        self.goFreqs.to_csv(self.results_path + 'go_freqs.csv')

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

        if (self.goFreqs is None) and ("go_freqs.csv" not in os.listdir(self.results_path)):
            raise TypeError("No GO IDs present")
        elif self.goFreqs is None:
            self.goFreqs = pd.read_csv(self.results_path + "go_freqs.csv")
        else:
            pass

        if self.goFreqs.index.name is None:
            self.goFreqs.set_index("Gene Ontology ID", inplace=True)
        self.biological = self.goFreqs.loc[self.goFreqs["Gene Ontology Classification"].isin(["Biological Process"])]
        self.cellular = self.goFreqs.loc[self.goFreqs["Gene Ontology Classification"].isin(["Cellular Function"])]
        self.molecular = self.goFreqs.loc[self.goFreqs["Gene Ontology Classification"].isin(["Molecular Function"])]

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
        plt.rcParams.update({'font.size': 24 })
        plt.subplots_adjust(hspace=.0)
        data_labels = ['Biological Process', 'Cellular Function', 'Molecular Function']
        color_list = ['#67cceb', '#abd12f', '#d1845e']

        fig_save_loc = self.results_path + 'gene_ontology_distribution.jpeg'

        dataframes = [self.biological, self.cellular, self.molecular]
        for i, data in enumerate(dataframes):
            sns.scatterplot(ax=axs[i], data=data[0:10], x="Score (% of genes involved)", y="Classification Value",
                            size="Count", legend=False, sizes=(20, 1000), color=color_list[i])
            axs[i].set_ylabel(data_labels[i])
            axs[i].grid(visible=False)
        plt.savefig(fig_save_loc, dpi=512, bbox_inches='tight', pad_inches=0.2)


if __name__ == '__main__':
    source_dir = input('Enter parent directory: ')
    temp_obj = GeneOntology(source_dir=source_dir)
    temp_obj.process_go_ids()
    temp_obj.go_class_splitter()
    temp_obj.plotGeneOntology()

