import os
import shutil

def runAntiSmash(source_dir: str):
    antismash_dir = source_dir + 'antiSMASH_dset/'

    

def create_antiSMASH_dset(source_dir: str):
    """
        Sift through Prokka generated data and grab .gbf files, rename them and create a new directory
        for antismash analysis.

        Returns:
            None

    """
    # do a clean create
    antismash_dir = source_dir + 'antiSMASH_dset/'
    prokka_run_dir = source_dir + 'prokka_run_data/'
    if os.path.exists(antismash_dir):
        shutil.rmtree(antismash_dir)

    os.mkdir(antismash_dir)

    for folder in os.listdir(prokka_run_dir):
        if not os.path.isdir(prokka_run_dir + folder):
            continue
        for file in os.listdir(prokka_run_dir + folder + '/'):
            if '.gbf' in file:
                shutil.copy(src=prokka_run_dir + folder + '/' + file,
                            dst=antismash_dir)
                # rename
                shutil.move(antismash_dir + file,
                            dst=antismash_dir + folder + '.gbk')