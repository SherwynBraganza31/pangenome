"""
    Genomes downloaded from NCBI get zipped into a hierarchy of folders. Moreover,
    each assembly is named by its accession/assembly ID. The aim of this module is
    to remove the element of human interaction needed to rename the accessions to
    the species name.

    The folder hierarchy is as follows:
    ncbi_dataset / data / <assembly_ID> / <assembly_ID>.fna

    The relations between assembly_ID and organism are stored in a JSON List file -
    assemlby_data_report.jsonl. The data in this file will be used to rename the
    .fna files according to the requirements of the following annotation/PPaNGGOLiN
    modules.

"""

import json
import os
import shutil


def generateFASTA(source: str, dest: str):
    """
    Parsing function that reads through the downloaded source directory tree and extracts the fasta files
    to format them according the requirements of the following modules.

    Input (unzipped NCBI dataset folder) ----> Output (Ideally named prepped_fasta)

    Args:
        source: The source folder in the above described format to retrieve the
                fasta files from
        dest: The location where the output folder should be created

    Returns:

    """

    def createDSET(linking_dict, source, dest):
        dir_list = os.listdir(source)
        for x in linking_dict.keys():
            if x in dir_list:
                try:
                    subdir_list = os.listdir(source + x)
                    fasta_filename = [x for x in subdir_list if '.fna' in x][0]
                    shutil.copy(source + x + '/' + fasta_filename,
                                dest + linking_dict[x] + '.fna')
                # If source and destination are same
                except shutil.SameFileError:
                    print("Source and destination represents the same file.")

                # If there is any permission issue
                except PermissionError:
                    print("Permission denied.")

    if not os.path.exists(source):
        print('Source folder does not exist.')
        raise FileNotFoundError

    if not os.path.exists(dest):
        os.mkdir(dest)

    with open(source + 'linking_dict.json') as ifile:
        linking_dict = json.load(ifile)
    createDSET(linking_dict, source, dest)


if __name__ == '__main__':
    print('This module is designed to automatically rename and process the FASTA files'
          ' needed for further modules. \nA few requirements need to be met first:\n'
          '1) The directory containing the datasets from NCBI must be unzipped.\n'
          '2) The name of the folder to be processed must be entered.\n')

    directory_name = input('Enter the directory name to be processed: ')
    directory_name = directory_name if "/" == directory_name[-1] else directory_name + "/"
    output_folder = input('Enter output directory: ')
    output_folder = output_folder if "/" == output_folder[-1] else output_folder + "/"

    assembly_path = directory_name + 'ncbi_dataset/data/'

    generateFASTA(assembly_path, output_folder)
    print('Number of genomes processed : {}'.format(len(os.listdir(output_folder))))
