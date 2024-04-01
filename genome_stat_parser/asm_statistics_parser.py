import pandas as pd
import os
from io import StringIO

class InvalidOutputFormat(Exception):
    "Raised when the save file output is unknown"
    pass

class AssembStatParser:
    def __init__(self, directory: str = '.'):
        self.dir = directory
        self.dir_list = os.listdir(self.dir)
        self.meta_data = self.create_metadata_container()

        for x in self.dir_list:
            if '.txt' in x:
                self.meta_data = pd.concat([self.meta_data, self.extract_data(x, self.meta_data.columns)],
                                           ignore_index=True)

    def create_metadata_container(self) -> pd.DataFrame:
        col_names = ['Organism name',
                     'Strain',
                     'Assembly name',
                     'Assembly Level',
                     'Size(Mb)',
                     'Taxid',
                     'RefSeq assembly accession',
                     'Genome Representation',
                     'total-length',
                     'Coverage',
                     'Scaffolds',
                     'Scaffold N50',
                     'Scaffold L50',
                     'Contigs',
                     'Total Gap Length',
                     'Spanned Gaps',
                     'GC%',
                     'Study: cultured/metagenomics',
                     '16sRNA',
                     'F',
                     'M',
                     'CDS',
                     'Gene']

        return pd.DataFrame(columns=col_names)

    def extract_data(self, filename: str, column_names: list[str]) -> pd.DataFrame:
        f = open(filename, "r")
        core_data = pd.DataFrame(columns=['unit-name', 'molecule-name', 'molecule-type/loc',
                                          'sequence-type', 'statistic', 'value'])
        temp_meta_data = pd.DataFrame(columns=column_names)
        temp_meta_data.loc[0] = ['']*len(column_names)
        for x in f:
            if x[0] != '#':
                line_data = pd.read_csv(StringIO(x), sep='\t', header=None)
                line_data.columns = core_data.columns
                core_data = pd.concat([core_data, line_data], ignore_index=True)
            else:
                for col in column_names[0:8]:
                    if col.casefold() in x.casefold():
                        start_idx = x.index('=') + 1 if col == 'Strain' else x.index(':') + 2
                        temp_meta_data.loc[0, col] = x[start_idx: -1]

        return temp_meta_data

    def save_dataframe(self, filename:str):
        if '.csv' in filename:
            self.meta_data.to_csv(filename)
        elif '.json' in filename:
            self.meta_data.to_json(filename)
        elif '.xlsx' in filename:
            self.meta_data.to_excel(filename)
        else:
            raise InvalidOutputFormat


if __name__ == '__main__' :
    temp = AssembStatParser()
    temp.save_dataframe('statistics.csv')
