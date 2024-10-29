from fasta_folder_gen import generateFASTA
from prokka_handler import ProkkaHandler

class AnnotationController:
    def __init__(self):
        self.source_dir = self.get_source_dir()
        self.outdir = self.source_dir + "prepped_fasta/"
        self.check_source_dir()
        self.assembly_path = self.source_dir + "ncbi_dataset/data/"

        generateFASTA(self.assembly_path, self.outdir)
        prokka_handler = ProkkaHandler(parent_dir=self.source_dir,
                                       fasta_dir=self.outdir)

        return

    def get_source_dir(self):
        #TODO
        # get input directory from user
        # it should be of the file heirarchy as follow
        # | source_dir
        # |     |-- ncbi_dataset
        # |         |-- data
        # |             assembly_report.jsonl
        # |             data_summary.tsv
        # |             dataset_catalog.json
        # |             |-- GCF_xxxxxxx
        # |                 xxx.fasta
        # |             |-- GCF_xxxxxxx
        return

    def check_source_dir(self):
        return