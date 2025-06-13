import time
import os
from prokka_handler import ProkkaHandler
import shutil
import numpy as np
import json
from sklearn import model_selection

parent_dir = '/home/sbraganza/projects/annotation_testset/'

fasta_filenames = list(os.listdir(parent_dir))
time_dict = {}


def bootstrap_run(iterations, set_size):
    bootstrap = model_selection.ShuffleSplit(n_splits=iterations, test_size=set_size, train_size=None)
    sp_runtime = []
    mp_runtime = []

    # get a bootstrap subset  and run prokka on it
    for count, (waste, file_idx) in enumerate(bootstrap.split(fasta_filenames)):
        print(f'Working on subset {count}')
        os.mkdir(parent_dir + 'subset/')

        for file in np.array(fasta_filenames)[file_idx.astype(int)]:
            shutil.copy(src=parent_dir + file, dst=parent_dir + 'subset/' + file)

        temp = ProkkaHandler(parent_dir=parent_dir + 'subset/', fasta_dir=parent_dir + 'subset/')
        temp.createRunDataFolder()

        # run and record multiprocess time
        start = time.time()
        temp.executeProkkaCalls_mp()
        mp_runtime.append(time.time()-start)

        # run and record singleprocess time
        start = time.time()
        temp.executeProkkaCalls()
        sp_runtime.append(time.time()-start)

        # remove subset directory created for new batch
        shutil.rmtree(parent_dir + 'subset/')

    return {"sp": sp_runtime, "mp": mp_runtime}


for bootstrap_size in range(40, 130, 10):
    print(f'Working on Bootstrap Size {bootstrap_size}')
    time_dict.update({bootstrap_size: bootstrap_run(5, bootstrap_size)})
    with open(parent_dir + 'time_dict.json', 'w') as file:
        json.dump(time_dict, file)

