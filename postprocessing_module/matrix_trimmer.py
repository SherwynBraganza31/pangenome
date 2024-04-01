import pandas as pd
import json
import numpy as np
from tqdm import tqdm
import multiprocess as mp

class GeneParser:
    """
        Extract Genes from the matrix file and map gene families,
        UniProtKB(s)

        Attributes:
        -----------
        parsed_genedata_path: str
            path to where generated data files will be stored
        ppan_dir: str
            path to the ppanggolin generated dataset
        uniprot_path

        Storage Tree:

                                parent_directory
                                       |
                                       |
         PPaNGGOLIN_RUN (ppan_dir)    --  trimmed_matrix_files   --    prokka_run_data  -- ......
                    /                           \
                   /                             \
                  /                               \
                matrix                      genes.csv ; hypothetical.csv ; recognized.csv ; gff_seqeuencing.csv



    """
    def __init__(self, parent_path: str, ppan_dir: str, uniprot_path: str):
        self.parent_path = parent_path if parent_path[-1] == '/' else parent_path + '/'
        self.ppan_dir = parent_path + "PPaNGGOLiN_RUN/"
        self.matrix_path = self.ppan_dir + '/matrix/matrix.csv'
        self.parsed_genedata_path = self.parent_path + '/trimmed_matrix_files/'
        self.uniprot_path = self.parsed_genedata_path

    def parseGenes(self):
        """
            Grab a reduced version of the PPanGGOLiN's matrix file with the following fields -
            'Gene ID', 'Organism', 'Gene Family', 'Genome ID', 'UniProtKB' and 'Partition'
            Each row in a matrix file starts with the name of a Gene Family. The latter columns
            correspond to all the organisms under consideration and contain the names of the genes
            present in that organism that belong to the same gene family.

            The aim of this function is to split these genes up as its own entity and have the
            gene family as an attribute.

            In addition, it links up UniProtIDs for each gene listed.

            Saves the generated table as genes.csv.

            Returns: None

        """
        matrix = pd.read_csv(self.matrix_path)
        self.genes = pd.DataFrame(columns=['Gene ID', 'Organism', 'Gene Family', 'Genome ID', 'UniProtKB', 'Partition'])

        organism_idx = list(
            range(matrix.columns.to_list().index('Avg group size nuc') + 1, len(matrix.columns.to_list())))
        matrix_cols = matrix.columns.to_numpy()
        matrix = matrix.to_numpy()
        pbar = tqdm(total=matrix.shape[0] * matrix.shape[1])    # progress bar

        for row in range(0, matrix.shape[0]):
            for col in organism_idx:
                gene_family = matrix[row, 0]
                genome = gene_family.split('_')[0]
                organism = matrix_cols[col]
                partition = matrix[row, 1]
                uniprot = None
                # skip if the organism doesn't share a gene with family
                if pd.isna(matrix[row, col]):
                    continue
                # it may share a list of genes, split them and process each one
                gene_list = str(matrix[row, col]).replace('\"', "").split(' ')
                for gene in gene_list:
                    self.genes.loc[len(self.genes.index)] = [gene, organism, gene_family, genome, uniprot, partition]
                pbar.update(1)
        pbar.close()

        self.linkUniProt()
        self.genes.set_index('Gene ID', inplace=True)
        self.genes.to_csv(self.parsed_genedata_path + 'genes.csv')

    def parseGenes_mp(self):
        """
        Multiprocess version of the above code.

        TODO Not complete. Experimental version
        """
        whole_matrix = pd.read_csv(self.matrix_path)
        output = pd.DataFrame(columns=['Gene ID', 'Organism', 'Gene Family', 'Genome ID', 'UniProtKB', 'Partition'])

        def splitGenes(idx, matrix):
            out_df = pd.DataFrame(columns=['Gene ID', 'Organism', 'Gene Family', 'Genome ID', 'UniProtKB', 'Partition'])
            pbar = tqdm(total=matrix.shape[0] * matrix.shape[1])
            for row in range(0, matrix.shape[0]):
                for col in organism_idx:
                    gene_family = matrix[row, 0]
                    genome = gene_family.split('_')[0]
                    organism = matrix_cols[col]
                    partition = matrix[row, 1]
                    uniprot = ""

                    # skip if the organism doesn't share a gene with family
                    if pd.isna(matrix[row, col]):
                        continue
                    # it may share a list of genes, split them and process each one
                    gene_list = str(matrix[row, col]).replace('\"', "").split(' ')
                    for gene in gene_list:
                        out_df.loc[len(out_df.index)] = [gene, organism, gene_family,
                                                         genome, uniprot, partition]
                    pbar.update(1)

            pbar.close()
            # output.set_index('Gene ID', inplace=True)
            # output.to_csv(self.parsed_genedata_path + f'filtered_matrix_{idx}.csv')
            return out_df

        organism_idx = list(
            range(whole_matrix.columns.to_list().index('Avg group size nuc') + 1, len(whole_matrix.columns.to_list())))
        matrix_cols = whole_matrix.columns.to_numpy()
        whole_matrix = whole_matrix.to_numpy()
        max_rows = whole_matrix.shape[0]
        cpus = mp.cpu_count()
        input(f"Found {cpus} CPU Cores. Hit Enter to continue...")
        splits = [x * max_rows // cpus for x in range(cpus)] + [max_rows]

        with mp.Pool() as p:
            return_val = p.starmap(splitGenes, [(x, whole_matrix[splits[x]:splits[x + 1]]) for x in range(0, 5)])

        output = pd.concat([output]+return_val, axis=0, ignore_index=True)
        self.genes = output
        self.linkUniProt()
        self.genes.set_index('Gene ID', inplace=True)
        self.genes.to_csv(self.parsed_genedata_path + 'genes.csv')

    def linkUniProt(self):
        """
        Matches the UniProt values with the corresponding genes. Matchings lie in the .csv file
        generated by the gff_parser.

        Returns: None
        """
        save_flag = False
        if self.genes is None:
            self.genes = pd.read_csv(self.parsed_genedata_path + 'genes.csv')
            save_flag = True
        columns, genes_np = self.genes.columns.to_list(), self.genes.to_numpy()
        pbar = tqdm(total=self.genes.shape[0], unit='matches')

        uniprot_db = pd.read_csv(self.uniprot_path + '/gff_sequencing.csv').to_numpy()
        # create a map from gene_id to UniProtKB built from gff sequencing
        uniprot_db = dict(zip(uniprot_db[:, 0], uniprot_db[:, 2]))

        def uniprot_match(x, progress_bar):
            element_match = ''
            try:
                element_match = uniprot_db[x]
            except KeyError:
                pass
            progress_bar.update(1)
            return str(element_match)
        vect_uni_match = np.vectorize(uniprot_match)
        genes_np[:, 4] = vect_uni_match(genes_np[:, 0], pbar)

        self.genes = pd.DataFrame(data=genes_np, columns=self.genes.columns.to_list())
        if save_flag:
            self.genes.set_index('Gene ID', inplace=True)
            self.genes.to_csv(self.parsed_genedata_path + 'genes.csv')

    def splitProteins(self):
        """
        Splits the self.genes dataframe into two new DataFrames - hypothetical and recognized.
        Hypotheticals - Genes that do not have a linked UniProt Id
        Recognized - Genes that are assigned a UniProt ID.

        Returns:
            None
        """
        self.hypothetical = pd.DataFrame(columns=self.genes.columns)
        self.recognized = pd.DataFrame(columns=self.genes.columns)
        pbar = tqdm(total=self.genes.shape[0])
        uniprot_freqs = {}

        for idx, rows in self.genes.iterrows():
            if pd.isna(rows.loc['UniProtKB']) or rows.loc['UniProtKB'] == (None or 'nan'):
                self.hypothetical = pd.concat([self.hypothetical, rows], axis=0, ignore_index=True)
            else:
                self.recognized = pd.concat([self.recognized, rows], axis=0, ignore_index=True)
                uniprotKB = rows.loc['UniProtKB']
                if uniprotKB in uniprot_freqs.keys():
                    uniprot_freqs[uniprotKB] += 1
                else:
                    uniprot_freqs.update({uniprotKB : 1})
            pbar.update(1)

        pbar.close()
        self.hypothetical.set_index('Gene ID', inplace=True)
        self.hypothetical.to_csv(self.parsed_genedata_path + 'hypotheticals.csv')
        self.recognized.set_index('Gene ID', inplace=True)
        self.recognized.to_csv(self.parsed_genedata_path + 'recognized.csv')

        with open(self.parent_path + "/Gene_Ontology/uniprot_freqs.json", "w", encoding='utf-8') as file:
            json.dump(uniprot_freqs, file)

    def generateUniProtFreqs(self):
        """

        columns = ['Gene ID', 'Organism', 'Gene Family', 'Genome ID', 'UniProtKB', 'Partition']

        Returns:

        """

        matrix = self.recognized.to_numpy()
        uniprots = {}
        for row in range(matrix.shape[0]):
            if matrix[row, 4] not in uniprots.keys():
                uniprots.update({matrix[row, 4]: 1})
            else:
                uniprots[matrix[row, 4]] += 1

        with open(self.parent_dir + "/Gene_Ontology/uniprot_freqs.json", "w") as file:
            json.dump(uniprots, file)


if __name__ == '__main__':
    temp = GeneParser(parent_path="../geobacillus_thermodenitrificans_pangenome",
                      ppan_dir="../geobacillus_thermodenitrificans_pangenome/PPaNGGOLiN_RUN",
                      uniprot_path="../geobacillus_thermodenitrificans_pangenome/trimmed_matrix_files")
    temp.parseGenes_mp()
    temp.linkUniProt()
    # temp.splitProteins()
