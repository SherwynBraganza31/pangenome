import os
import shutil
def create_dbCAN_dset(source_dir: str):
    """
        Sift through Prokka generated data and grab .gbf files, rename them and create a new directory
        for dbCAN analysis.

        Returns:
            None

    """
    # do a clean create
    dbCAN_dir = source_dir + 'dbCAN_dset/'
    prokka_run_dir = source_dir + 'prokka_run_data/'
    if os.path.exists(dbCAN_dir):
        shutil.rmtree(dbCAN_dir)

    os.mkdir(dbCAN_dir)

    for folder in os.listdir(prokka_run_dir):
        if not os.path.isdir(prokka_run_dir + folder):
            continue
        for file in os.listdir(prokka_run_dir + folder + '/'):
            if '.faa' in file:
                shutil.copy(src=prokka_run_dir + folder + '/' + file,
                            dst=dbCAN_dir)
                # rename
                shutil.move(dbCAN_dir + file,
                         dst=dbCAN_dir + folder + '.faa')


def run_dbCAN()