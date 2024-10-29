import os
import subprocess
import pandas as pd
import json


class BuscoHandler:
    """
    Class that handles running of BUSCO and updates the data_summary.tsv file with the results from the
    BUSCO run.

    Attributes
    ----------
    working_dir: str
        zoned in directory that contains all subdirectories of all the assemblies
    busco_dir: str
        new directory to be created to store BUSCO run results


    Methods
    -------
    run_BUSCO()
        handles running of BUSCO as a subprocess and storing/hashing the generated data

    update_BUSCO_scores()
        handles modifying the data_summary.tsv file with the BUSCO scores
    """
    def __init__(self, working_dir: str):
        """
        :param working_dir: str
            The directory in which all the accession fasta files are located
        :param busco_dir:
            The new directory to be created to store busco run results
        """
        self.working_dir = working_dir if working_dir[-1] == '/' else working_dir + '/'
        self.score_dict = {} # dictionary that scores busco

    def run_BUSCO(self):
        """
        Run BUSCO as a subprocess for all the assembly accessions in the source directory.

        :return: None
        """
        running_dir = os.getcwd()
        os.chdir(self.working_dir)

        cpus = input("Enter number of cpus to use (default=8): ")
        if cpus == '' or cpus is None:
            cpus = '8'

        lineage = input("Enter lineage (default is bacteria_odb10): ")
        if lineage == '' or lineage is None:
            lineage = 'bacteria_odb10'

        for x in os.listdir(self.working_dir):  # parse through all assembly folders and files in dir
            if os.path.isdir(self.working_dir + x):  # filter to only get assembly directories
                for fasta_file in os.listdir(self.working_dir + x):
                    if '.fna' in fasta_file:  # only pick the fasta file in assembly dir
                        # run BUSCO analysis
                        subprocess.run(f'busco -i {self.working_dir}{x}/{fasta_file} --mode genome --cpu {cpus}' +
                                       f' --lineage_dataset {lineage} -o BUSCO/{x} -f',
                                       shell=True)
                        for json_file in os.listdir(f'BUSCO/{x}'):
                            if '.json' in json_file:
                                with open(f'BUSCO/{x}/{json_file}', encoding="utf-8") as f:
                                    json_dict = json.load(f)
                                    self.score_dict.update({x:json_dict['results']['Complete percentage']})

        with open(self.working_dir + 'busco_scores.json', 'w') as file:
            json.dump(self.score_dict, file)
        os.chdir(running_dir)  # return to the running directory
        return

    def update_BUSCO_scores(self):
        """
        Modify data_summary.tsv to include a column for BUSCO scores

        :return: None
        """
        df = pd.read_csv(f'{self.working_dir}/data_summary.tsv', sep='\t')
        df['BUSCO'] = [0] * df.shape[0]  # add an extra column of 0s
        columns = list(df.columns)
        df = df.to_numpy()

        # read in score dict
        # empty dictionarys get converted to False bools
        if not bool(self.score_dict):  # if dictionary is empty
            try:
                with open(self.working_dir + 'busco_scores.json', 'r') as file:
                    self.score_dict = json.load(file)
            except FileNotFoundError:
                print('Empty Busco Scores Dictionary and cannot find file to load from')
                return

        for key in self.score_dict.keys():
            df[df[:, columns.index('Assembly Accession')] == key, -1] = self.score_dict[key]

        df = pd.DataFrame(df, columns=columns)
        # shuffle column around so its after checkM
        checkM_idx = columns.index('Level')
        columns.insert(checkM_idx, columns.pop(len(columns) - 1))
        df.reindex(columns=columns)
        df.set_index('Assembly Accession', inplace=True)
        df.to_csv(f'{self.working_dir}/data_sumamry.tsv', sep='\t')
        return


# if __name__ == '__main__':
#     busco_handler = BuscoHandler(working_dir='/home/sbraganza/projects/171genomes/171genomes_ncbi/ncbi_dataset/data')
#     placeholder = input("Hit Enter to start BUSCO analysis...\n")
#     busco_handler.run_BUSCO()
#     busco_handler.update_BUSCO_scores()
