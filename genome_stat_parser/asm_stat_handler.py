import json
import pandas as pd
import numpy as np
import matplotlib
import os
import shutil


class ASMStatParser:
    """
    Class that extracts data from assembly_data_report.jsonl and fills data in
    data_summary.tsv

    Attributes
    ----------
    source_dir: str
        zoned in directory that contains all subdirectories of all the assemblies
    jsonlist: list[dict]
        list of dictionaries containing assembly report data for each accession


    Methods
    -------
    extract_checkM_stats()
        Extract the checkM scores from assembly_report.jsonl

    generate_ranker_columns()
        Generate arbitrary ranker columns for data_summary2.tsv

    generate_plots()

    convertJSONlines()
    """
    def __init__(self, source_dir: str):
        # source directory error checking
        self.source_dir = source_dir + '/' if source_dir[-1] != '/' else source_dir
        # TODO: check if source directory has a the correct dataset folder tree
        # source_dir -> ncbi_dataset -> data -> (assembly_data_report.jsonl, dataset_summary.tsv, dataset_catalog.json)


        self.jsonlist = list(open(self.source_dir + 'assembly_data_report.jsonl', 'r', encoding='utf-8'))
        self.curated_genome_frame = None

    def extract_checkM_stats(self):
        """
            Extract the checkM scores from assembly_data_report.jsonl and fill it into
            data_summary_v2.csv
            :return:
        """
        dataframe = pd.read_csv(self.source_dir + 'dataset_summary_v2.csv')
        dataframe.set_index('Assembly Accession', inplace=True)
        dataframe['CheckM Score'] = [''] * dataframe.shape[0]
        dataframe_cols = list(dataframe.columns)

        # dive into json to find checkM percent
        for idx in range(len(self.jsonlist)):
            json_dict = json.loads(self.jsonlist[idx])
            accession = json_dict['accession']
            try:
                dataframe.loc[accession, "CheckM Score"] = json_dict['checkmInfo']['completeness']
            except KeyError:
                print(f'No completeness % for {accession}')

        # move columns around
        checkM_idx = dataframe_cols.index('Level')
        dataframe_cols.insert(checkM_idx + 1, dataframe_cols.pop(len(dataframe_cols) - 1))
        dataframe.reindex(columns=dataframe_cols)
        dataframe.to_csv(self.source_dir + 'modded_data_summary.csv')

    def generate_ranker_columns(self):
        data = pd.read_csv(f'{self.source_dir}modded_data_summary.csv')
        data['Ranker CheckM-BUSCO'] = [0] * data.shape[0]
        data['Ranker BUSCO-CheckM'] = [0] * data.shape[0]
        data['Ranker Equal'] = [0] * data.shape[0]
        data_columns = list(data.columns)
        data = data.to_numpy()
        level_encoding = {'Complete Genome': 3, 'Chromosome': 2, 'Scaffold': 1, 'Contig': 0}
        # encode_ranker = np.vectorize(lambda a, b: level_encoding[a] * 1000 + b)
        data[:, -3] = data[:, data_columns.index('CheckM Score')] * 1000 + data[:, data_columns.index('BUSCO')]
        data[:, -2] = data[:, data_columns.index('CheckM Score')] + data[:, data_columns.index('BUSCO')] * 1000
        data[:, -1] = data[:, data_columns.index('CheckM Score')] + data[:, data_columns.index('BUSCO')]
        sorted_indices = [np.argsort(data[:, -1])[-1::-1], np.argsort(data[:, -2])[-1::-1],
                          np.argsort(data[:, -3][-1::-1])]
        filenames = ['data_summ_CheckM.tsv', 'data_summ_BUSCO.tsv', 'data_summ_equal.tsv']
        for i in range(3):
            sorted_index = np.argsort(data[:, -(3-i)])[-1::-1]
            data = data[sorted_index, :]
            df = pd.DataFrame(data, columns=data_columns)
            df.set_index('Assembly Accession', inplace=True)
            df.to_csv(f'{self.source_dir}{filenames[i]}', sep='\t')
        return

    def curate_genomes(self, busco_cutoff: float = 95, checkM_cutoff: float = 99.5):
        df = pd.read_csv(self.source_dir + 'modded_data_summary.csv')
        self.curated_genome_frame = df[((df['Level'].isin(['Complete Genome', 'Chromosome']))
                                        | (df['CheckM Score'] > checkM_cutoff) & (df['BUSCO'] > busco_cutoff))]

        if 'assembly_data_report.json' not in os.listdir(self.source_dir):
            self.convert_json_lines()

        linking_dict = dict(zip(df['Assembly Accession'].values.tolist(),
                                np.zeros(self.curated_genome_frame.shape[0])))

        def fill_linking_dict():
            """
            Generates a mapping between the accession id and the name of the organism

            Args:
                linking_dict: Mapping between accession id and org name
                source: The source folder that contains the accessions

            Returns:
                None
            """
            with open(self.source_dir + 'assembly_data_report.json') as ifile:
                jsonlist = json.load(ifile)

            for x in jsonlist:
                accession = x['accession']
                if accession not in linking_dict.keys():
                    pass

                orgname = str.split(x['organism']['organismName'], sep=' ')

                # Strain doesn't always exist? In that case, empty strain
                try:
                    strain = x['organism']['infraspecificNames']['strain']
                except KeyError:
                    strain = ''

                # Sometimes the strain is present in the original name.
                # Remove if persists to prevent duplication
                try:
                    orgname.remove(strain)
                except ValueError:
                    pass

                # handle case if there is no species
                try:
                    species = orgname[1][:-1] if orgname[1][-1] == '.' else orgname[1]
                except IndexError:
                    species = 'sp'

                genus = orgname[0]
                fullname = genus + '_' + species + '_' + strain

                while fullname[-1] == '_':
                    fullname = fullname[:-1]

                linking_dict.update({accession: fullname})

        fill_linking_dict()
        self.curated_genome_frame['Filename'] = \
            [linking_dict[x] for x in self.curated_genome_frame['Assembly Accession']]

        def drop_dupes():
            uniques = np.unique(self.curated_genome_frame['Filename'])
            for x in uniques:
                dupes = self.curated_genome_frame[self.curated_genome_frame['Filename'].isin([x])]
                best = dupes['BUSCO'].idxmax()
                self.curated_genome_frame = self.curated_genome_frame.drop(dupes.drop(best).index)

        drop_dupes()

        with open(self.source_dir + 'linking_dict.json', 'w') as ofile:
            json.dump(linking_dict, ofile)

        self.curated_genome_frame.to_csv(self.source_dir + 'curated_genome_frame.csv', index=False)

    def create_curated_dir(self):
        if 'curated_genome_list.json' not in os.listdir(self.source_dir):
            self.curate_genomes()

        source_root = ''
        for x in self.source_dir.split('/')[:-4]:
            source_root += x + '/'

    def generate_plots(self):
        data = pd.read_csv(f'{self.source_dir}data_summ_equal.tsv', sep='\t')
        data_columns = list(data.columns)
        data = data.to_numpy()

    def convert_json_lines(self):
        json_list = []
        for x in self.jsonlist:
            json_list.append(json.loads(x))
        with open(self.source_dir + 'assembly_data_report.json', 'w') as ofile:
            json.dump(json_list, ofile)
        

if __name__ == '__main__':
    temp = ASMStatParser(source_dir='/home/sbraganza/projects/171genomes/171genomes_ncbi/ncbi_dataset/data/')
    temp.extract_checkM_stats()
    temp.curate_genomes()

