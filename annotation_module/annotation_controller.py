from annotation_module.fasta_folder_gen import generateFASTA
from annotation_module.prokka_handler import ProkkaHandler
import os


class AnnotationController:
    def __init__(self):
        self.source_dir = self.get_source_dir()
        if self.source_dir is None:
            print('Ending program ......')
            return

        self.prepped_fasta_dir = self.source_dir + "prepped_fasta/"
        self.assembly_path = self.source_dir + "ncbi_dataset/data/"

        generateFASTA(self.assembly_path, self.prepped_fasta_dir)
        prokka_handler = ProkkaHandler(parent_dir=self.source_dir, fasta_dir=self.prepped_fasta_dir)
        prokka_handler.createRunDataFolder()
        prokka_handler.executeProkkaCalls_mp()
        print('\nCreating PPaNGGOLiN Dataset ....\n\n')
        prokka_handler.createPPaNdataset()
        return

    def get_source_dir(self):

        source_dir = input('Enter in the source directory that contains the dataset: ')

        if '~' in source_dir.split('/'):
            source_dir = source_dir.replace('~', os.path.expanduser('~'))

        source_dir = source_dir if source_dir[-1] == '/' else source_dir + '/'

        return check_source_dir(source_dir)

def check_source_dir(dir):
    # resolve ~ if in string
    if "ncbi_dataset" not in os.listdir(dir):
        print('Source Directory does not have the correct file hierarchy')
        return None
    else:
        return dir
