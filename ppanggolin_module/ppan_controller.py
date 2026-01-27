import subprocess as sub
import os

class PpanController:
    def __init__(self, source_dir):

        if '~' in source_dir.split('/'):
            source_dir = source_dir.replace('~', os.path.expanduser('~'))


        source_dir = source_dir if source_dir[-1] == '/' else source_dir + '/'
        self.source_dir = source_dir
        os.makedirs(self.source_dir + 'ppanggolin_run', exist_ok=True)

        self.base_filename = "clustering_results"
        self.cpu = str(os.cpu_count() - 3)

        self.base_call()
        self.pangenome_data_gen()
        self.curves_matrix()
        self.rgp_gen()
        self.spots_of_insertion()
        self.rarefaction_curves()
        self.partitions_projections_families()
        self.fasta_outputs()
        self.phylogeny()

        print('PPanGGOLiN Module has done generating its files.\n\n')

    def base_call(self):
        # Change directory & create a clustered .h5 file
        print('Running base clustering commands....')
        os.chdir(self.source_dir + 'ppan_dataset')
        sub.run('ppanggolin annotate --anno organisms.gff.list --fasta organisms.fasta.list --output ' +
                       self.source_dir + 'ppanggolin_run --basename ' + self.base_filename + ' -f --cpu ' + self.cpu, shell=True)
        os.chdir(self.source_dir + 'ppanggolin_run')
        sub.run('ppanggolin cluster -p ' + self.base_filename + '.h5 --identity 0.5 --coverage 0.8 -f --cpu ' + self.cpu, shell=True)
        print('Finished clustering.\n')

    def pangenome_data_gen(self):
        # Generate Pangenome Data
        print('Running commands to generated pangenome data....')
        os.chdir(self.source_dir + 'ppanggolin_run')
        sub.run('ppanggolin graph -p ' + self.base_filename + '.h5 -f', shell=True)
        sub.run('ppanggolin partition -p ' + self.base_filename + '.h5 -f', shell=True)
        sub.run('ppanggolin info -p ' + self.base_filename + '.h5 --content', shell=True)
        print('Finished pangenome data generation.\n')

    def curves_matrix(self):
        # Generate Curves and Matrix Files
        print('Running commandst to generated Curves and Matrix File....')
        os.chdir(self.source_dir + 'ppanggolin_run')
        sub.run('ppanggolin draw -p ' + self.base_filename + '.h5 --ucurve --tile_plot --output '
                + self.base_filename + '_withcloud -f', shell=True)
        sub.run('ppanggolin draw -p ' + self.base_filename + '.h5 --tile_plot --nocloud --output '
                + self.base_filename + '_nocloud -f', shell=True)
        sub.run('ppanggolin write_pangenome -p ' + self.base_filename + '.h5 --csv --output matrix -f', shell=True)
        sub.run('ppanggolin write_pangenome -p ' + self.base_filename + '.h5 --stats --output stats -f', shell=True)
        sub.run('ppanggolin write_pangenome -p ' + self.base_filename + '.h5 --light_gexf --output light_gexf -f',
                shell=True)
        sub.run('ppanggolin write_pangenome -p ' + self.base_filename + '.h5 --gexf --output gexf -f',
                shell=True)
        sub.run('ppanggolin write_pangenome -p ' + self.base_filename + '.h5 --light_gexf --output light_gexf -f',
                shell=True)
        sub.run('ppanggolin write_pangenome -p ' + self.base_filename + '.h5 --partitions --output partitions -f',
                shell=True)
        print('Outputted matrix file and curves.\n')

    def rgp_gen(self):
        # Generate Regions of genome plasticity
        print('Running commands to generate regions of genome plasticity files....')
        os.chdir(self.source_dir + 'ppanggolin_run')
        sub.run('ppanggolin rgp -p ' + self.base_filename + '.h5', shell=True)
        sub.run('ppanggolin write -p ' + self.base_filename + '.h5 --regions --output RGP-Regions -f', shell=True)
        print('RGP files outputted.\n')

    def spots_of_insertion(self):
        print('Running commands to generate Spots Of Insertion File....')
        os.chdir(self.source_dir + 'ppanggolin_run')
        # Generate Spot Figures
        sub.run('ppanggolin spot -p ' + self.base_filename + '.h5 -f', shell=True)
        sub.run('ppanggolin write -p ' + self.base_filename + '.h5 --spots --output spots -f', shell=True)
        sub.run('ppanggolin draw -p ' + self.base_filename + '.h5 --spots all --output spot-figures -f', shell=True)
        print('Outputted spot figures.\n')

    def rarefaction_curves(self):
        print('Running commands to generate Rarefaction Curves Graph File....')
        os.chdir(self.source_dir + 'ppanggolin_run')
        # Generate Regions of genome plasticity
        sub.run('ppanggolin rarefaction -p ' + self.base_filename + '.h5 --output rarefaction -f', shell=True)
        print('Outputted rarefaction curves.\n')

    def partitions_projections_families(self):
        print('Running commands to generate partitions, projections, and gene family files....')
        os.chdir(self.source_dir + 'ppanggolin_run')
        # Generate Partitions, Projections, Gene Families
        sub.run('ppanggolin write -p ' + self.base_filename + '.h5 --partitions --output partitions -f', shell=True)
        sub.run('ppanggolin write -p ' + self.base_filename + '.h5 --projection --output projection -f', shell=True)
        sub.run('ppanggolin write -p ' + self.base_filename + '.h5 --families_tsv --output families-tsv -f', shell=True)
        print('Done generating the above files.\n')

    def fasta_outputs(self):
        print('Generating required FASTA outputs....')
        os.chdir(self.source_dir + 'ppanggolin_run')
        # Output FASTA files of Gene Data
        sub.run('ppanggolin fasta -p ' + self.base_filename + '.h5 --output MY_GENES --genes all -f', shell=True)
        sub.run('ppanggolin fasta -p ' + self.base_filename + '.h5 --output MY_GENES/persistent --genes persistent -f', shell=True)
        sub.run('ppanggolin fasta -p ' + self.base_filename + '.h5 --output MY_GENES/shell --genes shell -f', shell=True)
        sub.run('ppanggolin fasta -p ' + self.base_filename + '.h5 --output MY_GENES/cloud --genes cloud -f', shell=True)
        sub.run('ppanggolin fasta -p ' + self.base_filename + '.h5 --output MY_GENES/core --genes core -f', shell=True)
        sub.run('ppanggolin fasta -p ' + self.base_filename + '.h5 --output MY_GENES/softcore --genes softcore -f', shell=True)
        sub.run('ppanggolin fasta -p ' + self.base_filename + '.h5 --output MY_GENES/rgp --genes rgp -f', shell=True)

        # Output FASTA files of Protein Families
        sub.run('ppanggolin fasta -p ' + self.base_filename + '.h5 --output MY_PROT --prot_families all -f', shell=True)
        sub.run('ppanggolin fasta -p ' + self.base_filename + '.h5 --output MY_PROT/persistent --prot_families persistent -f', shell=True)
        sub.run('ppanggolin fasta -p ' + self.base_filename + '.h5 --output MY_PROT/shell --prot_families shell -f', shell=True)
        sub.run('ppanggolin fasta -p ' + self.base_filename + '.h5 --output MY_PROT/cloud --prot_families cloud -f', shell=True)
        sub.run('ppanggolin fasta -p ' + self.base_filename + '.h5 --output MY_PROT/core --prot_families core -f', shell=True)
        sub.run('ppanggolin fasta -p ' + self.base_filename + '.h5 --output MY_PROT/softcore --prot_families softcore -f', shell=True)
        sub.run('ppanggolin fasta -p ' + self.base_filename + '.h5 --output MY_PROT/rgp --prot_families rgp -f', shell=True)
        print('Done generating FASTA outputs.\n')

    def phylogeny(self):
        print('Generating phylogeny tree files....')
        os.chdir(self.source_dir + 'ppanggolin_run')
        sub.run('ppanggolin msa -p ' + self.base_filename + '.h5 --phylo -o "./Phylo" -f', shell=True)
        sub.run('ppanggolin msa -p ' + self.base_filename + '.h5 -o "./msa" -f', shell=True)
        print('Done generating phylogeny trees and its files.\n')

