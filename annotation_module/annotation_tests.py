import time
from torch.utils.data import DataLoader, Dataset
import os
from prokka_handler import ProkkaHandler
import shutil
import numpy as np
import json

parent_dir = '/home/sbraganza/projects/annotation_testset/'


class MyDataset(Dataset):
    def __init__(self, data):
        self.data = data

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        return self.data[idx]


fasta_filenames = list(os.listdir(parent_dir))

dataloader_10batch = DataLoader(dataset=MyDataset(fasta_filenames), batch_size=10, shuffle=True, drop_last=False)
dataloader_50batch = DataLoader(dataset=MyDataset(fasta_filenames), batch_size=50, shuffle=True, drop_last=False)
dataloader_full = DataLoader(dataset=MyDataset(fasta_filenames), batch_size=150, shuffle=True, drop_last=False)


def subset_run(dataloader):
    sp_runtime = []
    mp_runtime = []

    # create a subset of batch size and run prokka on it
    for count, data in enumerate(dataloader):
        print(f'Working on subset {count}')
        os.mkdir(parent_dir + 'subset/')

        for file in data:
            shutil.copy(src=parent_dir + file, dst=parent_dir + 'subset/' + file)

        temp = ProkkaHandler(parent_dir=parent_dir + 'subset/',
                             fasta_dir=parent_dir + 'subset/')
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


time_dict = {
    '10batch': subset_run(dataloader_10batch),
    '50batch': subset_run(dataloader_50batch),
    'full': subset_run(dataloader_full)
}

with open(parent_dir + 'time_dict.json', 'w') as file:
    json.dump(time_dict, file)
