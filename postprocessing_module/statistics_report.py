import pandas as pd
import numpy as np
import tables
import os
import json


class PostProcessing:
    def __init__(self, source_dir: str):
        self.source_dir = source_dir if source_dir[-1] == '/' else source_dir + '/'
        self.pangenome_dir = source_dir + 'ppanggolin_run/'
        self.post_proc_results = self.source_dir + 'postprocessing_results/'
        try:
            self.pangenome_filename = [x for x in os.listdir(self.pangenome_dir) if '.h5' in x][0]
        except IndexError:
            print('Directory provided does not contain a pangenome hdf5 file.')
            return


    def grab_pangenome_stats(self):
        def create_info_dict(info_group: tables.group.Group):
            """
            Read the pangenome content

            :param info_group: group in pangenome HDF5 file containing information about pangenome
            """
            attributes = info_group._v_attrs._f_list()

            info_dict = {"Genes": int(info_group._v_attrs['numberOfGenes'])}

            if "numberOfGenomes" in attributes:
                info_dict["Genomes"] = int(info_group._v_attrs['numberOfGenomes'])

            if "numberOfClusters" in attributes:
                info_dict["Families"] = int(info_group._v_attrs['numberOfClusters'])

            if "numberOfEdges" in attributes:
                info_dict["Edges"] = int(info_group._v_attrs['numberOfEdges'])

            if 'numberOfCloud' in attributes:  # then all the others are there

                persistent_stat = {"Family_count": int(info_group._v_attrs['numberOfPersistent'])}
                persistent_stat.update(info_group._v_attrs['persistentStats'])
                info_dict["Persistent"] = persistent_stat

                shell_stat = {"Family_count": int(info_group._v_attrs['numberOfShell'])}
                shell_stat.update(info_group._v_attrs['shellStats'])
                info_dict["Shell"] = shell_stat

                cloud_stat = {"Family_count": int(info_group._v_attrs['numberOfCloud'])}
                cloud_stat.update(info_group._v_attrs['cloudStats'])
                info_dict["Cloud"] = cloud_stat

                info_dict["Number_of_partitions"] = int(info_group._v_attrs['numberOfPartitions'])

                if info_group._v_attrs['numberOfPartitions'] != 3:
                    subpartition_stat = {f"Shell_{key}": int(val) for key, val in
                                         info_group._v_attrs['numberOfSubpartitions'].items()}
                    info_dict.update(subpartition_stat)

            if 'genomes_fluidity' in attributes:
                info_dict["Genomes_fluidity"] = {key: round(val, 3) for key, val in
                                                 info_group._v_attrs['genomes_fluidity'].items()}

            if 'family_fluidity' in attributes:
                info_dict["Family_fluidity"] = info_group._v_attrs['family_fluidity']

            if 'numberOfRGP' in attributes:
                info_dict["RGP"] = int(info_group._v_attrs['numberOfRGP'])

            if 'numberOfSpots' in attributes:
                info_dict["Spots"] = int(info_group._v_attrs['numberOfSpots'])

            if 'numberOfModules' in attributes:
                info_dict["Modules"] = {
                    'Number_of_modules': int(info_group._v_attrs['numberOfModules']),
                    'Families_in_Modules': int(info_group._v_attrs['numberOfFamiliesInModules']),
                    'Partition_composition': {
                        "Persistent": info_group._v_attrs['PersistentSpecInModules']['percent'],
                        "Shell": info_group._v_attrs['ShellSpecInModules']['percent'],
                        "Cloud": info_group._v_attrs['CloudSpecInModules']['percent']
                    }
                }
            return info_dict

        pangenome_data = tables.open_file(self.pangenome_dir + self.pangenome_filename, driver="H5FD_CORE")
        pangenome_info = create_info_dict(pangenome_data.root.info)
        with open(self.post_proc_results + 'pangenome_stats.json', 'w') as json_ofile:
            json.dump(pangenome_info, json_ofile)

        for x in pangenome_info.keys():
            print(f' {x}: {pangenome_info[x]}')

        pangenome_data.close()

    def get_gene_families(self):
        pangenome_data = tables.open_file(self.source_dir + self.pangenome_filename, driver="H5FD_CORE")
        gene_families = pangenome_data.get_node('/geneFamilies').read()




