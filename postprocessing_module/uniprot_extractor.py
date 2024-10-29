import os, gffutils, pandas as pd, json
from tqdm.auto import tqdm
import numpy as np


class UniProtExtractor:
    def __init__(self, source_dir):
        self.source_dir = source_dir if source_dir[-1] == '/' else source_dir + '/'
        self.gff_path = self.source_dir + 'ppan_dataset/GFF/'
        self.results_path = self.source_dir + 'postprocessing_results/'
        if not os.path.exists(self.results_path):
            os.makedirs(self.results_path)

    def parseGFF(self):
        file_list = os.listdir(self.gff_path)
        gene_count = {}
        out_df = pd.DataFrame(columns=['Gene ID', 'Gene Symbol', 'UniProtKB',
                                       'Annotation', 'Organism'])
        for file in file_list:
            count = 0
            if '.gff' not in file:
                continue
            # create a sqlite3 db of the file in memory
            print(f"Processing file: {file}")
            curr_file = gffutils.create_db(file, dbfn=':memory:',
                                           force=True, keep_order=True, merge_strategy='merge',
                                           sort_attribute_values=True)
            row_count = curr_file.execute("SELECT COUNT(*) FROM features").fetchone()[0]
            pbar = tqdm(total=row_count)
            # sql stmt to select attributes of all rows
            for row in curr_file.execute('SELECT attributes FROM features'):
                # parse string attributes into dict
                attribs = json.loads(row['attributes'])
                res = self.updateDatabase(out_df, attribs, file.split('.')[0])
                count += 1
                pbar.update(1)
            gene_count.update({file: count})

        uniprot_freqs = np.unique(out_df.loc[~out_df['UniProtKB'].isna(), ['UniProtKB']].to_numpy(), return_counts=True)
        uniprot_freqs = dict(zip(uniprot_freqs[0].astype(str), uniprot_freqs[1].astype(str)))
        with open(self.results_path + 'uniprot_freqs.json', 'w') as ofile:
            json.dump(uniprot_freqs, ofile)

        out_df.set_index('Gene ID', inplace=True)
        out_df.to_csv(self.results_path + 'gff_sequencing.csv')
        print(gene_count)
    def updateDatabase(db: pd.DataFrame, attributes: dict, organism: str) -> bool:
        """
          Takes in an attributes dictionary and adds data from it as a row
          in the db DataFrame.
          If the sub-attributes that its looking for doesn't exist,
          it skips the update and returns a False; else True.

          Attribute Checks:
          - Checks to see the dict contains data for the db cols (col_keys)
          - checks to see that inference is a list of 2 items -- The UniProtKB item
          is only present if inferece as the second element:

          All elements of the dictionary are lists of strings.

          @params
              db: pandas.DataFrame
                  DataFrame that represents the main 'database' of cherry
                  picked data
              attributes: dict
                  Dictionary corresponding to the attributes of a row of the
                  .gff file
              organims: str
                  The organism from which the gene came from

            @returns
              True:  if attributes dictionary passes all the checks and updates
              False: if it fails the checks

        """
        col_keys = ['ID', 'Name', 'inference', 'product']

        if 'ID' in attributes.keys():
            ID, Name, UniProtKB, product = None, None, None, None

            # split the string by :, the last element contains the UniProtKB
            for substring in attributes['inference']:
                if 'UniProtKB' in substring:
                    UniProtKB = substring.split(':')[-1]
                    break

            try:
                ID = attributes['ID'][0]
            except KeyError:
                ID = attributes['locus-tag']

            try:
                Name = attributes['Name'][0]
            except KeyError:
                pass

            try:
                product = attributes['product'][0]
            except KeyError:
                pass

            db.loc[len(db.index)] = [ID, Name, UniProtKB, product, organism]
            return True
        else:
            return False
