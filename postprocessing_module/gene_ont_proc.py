import json
import pandas as pd
import os

class GeneOntology:
    def __init__(self, parent_dir: str, results_path: str):
        self.parent_dir = parent_dir
        self.results_path = results_path
        self.ecNumbers = None
        self.goIDs = None

    def process_GO_IDs(self):
        with open(self.results_path + '/go_results.json', 'r', encoding='utf-8') as ifile:
            results = json.load(ifile)
        df = pd.DataFrame(columns=['Gene Ontology ID', 'Protein Name', 'UniProtKB',
                                   'Gene Ontology Classification', 'Classification Value'])
        ecNumbers = []
        for entry in results['results']:
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

                    df.loc[len(df.index)] = [go_id, prot_name, uniprot_id,
                                             go_class, go_class_value]

        self.goIDs = df
        self.ecNumbers = ecNumbers

        if "go_ids.csv" in os.listdir(self.results_path):
            df = pd.read_csv(self.results_path + "/go_ids.csv")
            self.goIDs = pd.concat([df, self.goIDs])

        self.goIDs.to_csv(self.results_path + '/go_ids.csv')

        with open(self.results_path + '/ecNumbers.txt', "a") as f:
            for x in self.ecNumbers:
                f.write(x + '\n')

        return

    def goClassSplitter(self):
        # Read entries from the .csv that links GO ID(s) to geneID(s)
        # The aim is to group up all UniProt ID(s) with the same GO ID
        # and then count the number of Genes associated with a particular
        # GO subfunction for each of the GO Function Categories - Biological Processes,
        # Molecular Functions and Cellular Functions
        # We then want to calculate the scores for each of these subprocesses, defined
        # as - (# genes involved in a subprocess) / (total # of genes)

        if not self.goIDs and "go_ids.csv" not in os.listdir(self.results_path):
            raise TypeError("No GO IDs present")
        elif self.goIDs:
            pass
        else:
            self.goIDs = pd.read_csv(self.results_path + "/go_ids.csv")

        with open(self.results_path + "/uniprot_freqs.json", "r", encoding="utf-8") as file:
            uniprot_freqs = json.load(file)

        # seperate each gene based on GO Process Class into a dictionary
        bio_dict, mol_dict, cel_dict = {}, {}, {}
        uniques_dict = [bio_dict, mol_dict, cel_dict]
        class_match = ['Biological Process', 'Molecular Function', 'Cellular Function']
        max_vals = [0] * 3

        for index, row in self.goIDs.iterrows():
            go_id = row.loc['Gene Ontology ID']
            go_class = class_match.index(row.loc['Gene Ontology Classification'])

            if go_id in uniques_dict[go_class].keys():
                uniques_dict[go_class][go_id][1] += uniprot_freqs[row.loc['UniProtKB']]
                max_vals[go_class] += uniprot_freqs[row.loc['UniProtKB']]
            else:
                uniques_dict[go_class].update({go_id: [row.loc['Classification Value'],
                                                       uniprot_freqs[row.loc['UniProtKB']],
                                                       0]})
                max_vals[go_class] += uniprot_freqs[row.loc['UniProtKB']]

        # Fill in Scores
        for go_class in range(len(class_match)):
            for key in uniques_dict[go_class]:
                uniques_dict[go_class][key][2] = round(uniques_dict[go_class][key][1] / max_vals[go_class] * 100, 3)

        for k in range(len(uniques_dict)):
            uniques_dict[k] = dict(sorted(uniques_dict[k].items(),
                                          key=lambda kv: kv[1][1], reverse=True))
            uniques_df = pd.DataFrame.from_dict(uniques_dict[k], orient='index',
                                                columns=['Metabolic Pathway',
                                                         'Gene Count', 'Score'])

            uniques_df.to_csv('/content/drive/MyDrive/geobacillus_pangenome/Gene_Ontology/' +
                              class_match[k].replace(' ', '_').lower() + '_go_count.csv')



