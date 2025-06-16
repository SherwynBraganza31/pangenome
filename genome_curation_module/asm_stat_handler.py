import json
import pandas as pd
import numpy as np
import os
from busco_handler import BuscoHandler
import shutil
import matplotlib as mpl
import matplotlib.pyplot as plt


class ASMStatParser:
    """
    Class that extracts data from assembly_data_report.jsonl and fills data in
    data_summary.tsv

    Attributes
    ----------
    working_dir: str
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
        # minor working directory error checking
        self.source_dir = source_dir + '/' if [-1] != '/' else source_dir
        self.working_dir = self.source_dir + 'ncbi_dataset/data/'

        self.jsonlist = list(open(self.working_dir + 'assembly_data_report.jsonl', 'r', encoding='utf-8'))
        self.curated_genome_frame = None

    def extract_checkM_stats(self):
        """
            Extract the checkM scores from assembly_data_report.jsonl and fill it into
            data_summary.csv
            :return:
        """
        dataframe = pd.read_csv(self.working_dir + 'dataset_summary.tsv', sep='\t')
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
        dataframe.to_csv(self.working_dir + 'data_summary.tsv', sep='\t')

    def get_busco_scores(self):
        busco_handler = BuscoHandler(working_dir=self.working_dir)
        busco_handler.run_BUSCO()
        busco_handler.update_BUSCO_scores()

    def curate_genomes(self, busco_cutoff: float = 95, checkM_cutoff: float = 99.5):
        df = pd.read_csv(self.working_dir + 'modded_data_summary.csv')

        # Choose genomes that meet:
        # (CheckM Level == Complete & Chromosome) | (CheckM Score > cutoff & BUSCO > cutoff)
        self.curated_genome_frame = df[((df['Level'].isin(['Complete Genome', 'Chromosome']))
                                        | (df['CheckM Score'] > checkM_cutoff) & (df['BUSCO'] > busco_cutoff))]

        if 'assembly_data_report.json' not in os.listdir(self.working_dir):
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
            with open(self.working_dir + 'assembly_data_report.json') as ifile:
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

        with open(self.working_dir + 'linking_dict.json', 'w') as ofile:
            json.dump(linking_dict, ofile)

        self.curated_genome_frame.to_csv(self.working_dir + 'curated_genome_frame.csv', index=False)

    def create_curated_dir(self):

        # error checking for if curate_genomes has not been called
        if 'curated_genome_list.json' not in os.listdir(self.working_dir):
            self.curate_genomes()

        curated_dir_subpath = self.source_dir + 'curated_genomes/data/'
        os.makedirs(curated_dir_subpath)

        curated_list = self.curated_genome_frame.index.to_list()

        for x in os.listdir(self.working_dir):
            if (x in curated_list or x == "curated_genome_frame.csv" or x == "assembly_data_report.jsonl"
                    or x == "linking_dict.json"):
                shutil.copytree(self.working_dir + "/" + x, curated_dir_subpath + '/' + x)

    def convert_json_lines(self):
        json_list = []
        for x in self.jsonlist:
            json_list.append(json.loads(x))
        with open(self.working_dir + 'assembly_data_report.json', 'w') as ofile:
            json.dump(json_list, ofile)

    def plot_genome_stats(self):
        with open(self.working_dir + 'linking_dict.json') as file:
            uncurated_filenames = np.asarray(list(json.load(file).values()))

        curated_genome_frame = pd.read_csv(self.source_dir + '58genomes/ncbi_dataset/data/curated_genome_frame.csv')
        curated_filenames = curated_genome_frame['Filename'].to_numpy()

        # splits the filenames into subspecies and returns assembly counts for each subspecies
        def parse_filenames_to_species(genome_names):
            subspecies = np.zeros(genome_names.shape).astype(str)
            grab_species = np.vectorize(lambda x: x.split('_')[1] if (
                        x.split('_')[1] != 'sp' and x.split('_')[1] != 'subsp') else 'no subspecies')
            species_distr, counts = np.unique(grab_species(genome_names), return_counts=True)
            counts = np.append(counts, [1, 1])

            return species_distr, counts

        uncurated_species, uncurated_counts = parse_filenames_to_species(uncurated_filenames)
        curated_species, curated_counts = parse_filenames_to_species(curated_filenames)
        mask = ~np.isin(uncurated_species, curated_species)
        curated_species, curated_counts = (
            np.hstack((curated_species, uncurated_species[mask])),
            np.hstack((curated_counts, np.zeros(
                shape=uncurated_species.shape[0] - curated_species.shape[0]))))

        sort_idx = np.argsort(curated_species)
        curated_species = curated_species[sort_idx]
        curated_counts = curated_counts[sort_idx]

        sort_idx = np.argsort(uncurated_species)
        uncurated_species = uncurated_species[sort_idx]
        uncurated_counts = uncurated_counts[sort_idx]

        def get_random_colors(n):
            colors = np.array(list(mpl.colors.CSS4_COLORS.values()))
            return np.random.choice(colors, size=n, replace=False)

        x = 6 * np.arange(len(uncurated_species))  # the label locations
        bar_width = 4  # the width of the bars

        fig, ax = plt.subplots(layout='constrained')

        rects = ax.bar(x, curated_counts, bar_width, label='Curated Genomes', color='#440154FF')
        ax.bar_label(rects, padding=3)
        rects = ax.bar(x + bar_width, uncurated_counts, bar_width, label='Uncurated Genomes', color='#FDE725FF')
        ax.bar_label(rects, padding=3)

        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax.set_ylabel('# of Assemblies')
        ax.set_title('Genome Asseblies Distribution')
        ax.set_xticks(x + bar_width / 2, uncurated_species, rotation=45)

        fig.savefig(self.source_dir + 'genome_distribution_sidebyside.jpeg', dpi=600)


if __name__ == '__main__':
    temp = ASMStatParser(source_dir='/home/sbraganza/projects/171genomes/171genomes_ncbi/ncbi_dataset/data/')
    temp.extract_checkM_stats()
    temp.curate_genomes()

