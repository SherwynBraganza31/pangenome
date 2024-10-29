import subprocess as sub
import os

class PpanController:
    def __init__(self, source_dir):
        self.source_dir = source_dir
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

    
    def base_call(self):
        # Change directory & create a clustered .h5 file
        os.chdir(self.source_dir + '/ppan_dataset')
        sub.run('ppanggolin annotate --anno organisms.gff.list --fasta organisms.fasta.list --output ' +
                       self.source_dir + '/ppanggolin_run --basename ' + self.base_filename + ' -f --self.cpu ' + self.cpu, shell=True)
        os.chdir(self.source_dir + '/ppanggolin_run')
        sub.run('ppanggolin cluster -p ' + self.base_filename + '.h5 --identity 0.5 --coverage 0.8 -f --cpu ' + self.cpu, shell=True)


    def pangenome_data_gen(self):
        # Generate Pangenome Data
        os.chdir(self.source_dir + '/ppanggolin_run')
        sub.run('ppanggolin graph -p ' + self.base_filename + '.h5 -f', shell=True)
        sub.run('ppanggolin partition -p ' + self.base_filename + '.h5 -f', shell=True)
        sub.run('ppanggolin info -p ' + self.base_filename + '.h5 --content', shell=True)


    def curves_matrix(self):
        # Generate Curves and Matrix Files
        os.chdir(self.source_dir + '/ppanggolin_run')
        sub.run('ppanggolin draw -p ' + self.base_filename + '.h5 --ucurve --tile_plot --output ' + self.base_filename + '_withcloud -f', shell=True)
        sub.run('ppanggolin draw -p ' + self.base_filename + '.h5 --tile_plot --nocloud --output ' + self.base_filename + '_nocloud -f', shell=True)
        sub.run('ppanggolin write -p ' + self.base_filename + '.h5 --gexf --csv --output matrix -f', shell=True)
        sub.run('ppanggolin write -p ' + self.base_filename + '.h5 --light_gexf --output light_gexf -f', shell=True)


    def rgp_gen(self):
        # Generate Regions of genome plasticity
        os.chdir(self.source_dir + '/ppanggolin_run')
        sub.run('ppanggolin rgp -p ' + self.base_filename + '.h5', shell=True)
        sub.run('ppanggolin write -p ' + self.base_filename + '.h5 --regions --output RGP-Regions -f', shell=True)


    def spots_of_insertion(self):
        os.chdir(self.source_dir + '/ppanggolin_run')
        # Generate Spot Figures
        sub.run('ppanggolin spot -p ' + self.base_filename + '.h5 -f', shell=True)
        sub.run('ppanggolin write -p ' + self.base_filename + '.h5 --spots --output spots -f', shell=True)
        sub.run('ppanggolin draw -p ' + self.base_filename + '.h5 --spots all --output spot-figures -f', shell=True)


    def rarefaction_curves(self):
        os.chdir(self.source_dir + '/ppanggolin_run')
        # Generate Regions of genome plasticity
        sub.run('ppanggolin rarefaction -p ' + self.base_filename + '.h5 --output rarefaction -f', shell=True)


    def partitions_projections_families(self):
        os.chdir(self.source_dir + '/ppanggolin_run')
        # Generate Partitions, Projections, Gene Families
        sub.run('ppanggolin write -p ' + self.base_filename + '.h5 --partitions --output partitions -f', shell=True)
        sub.run('ppanggolin write -p ' + self.base_filename + '.h5 --projection --output projection -f', shell=True)
        sub.run('ppanggolin write -p ' + self.base_filename + '.h5 --families_tsv --output families-tsv -f', shell=True)


    def fasta_outputs(self):
        os.chdir(self.source_dir + '/ppanggolin_run')
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


    def phylogeny(self):
        os.chdir(self.source_dir + '/ppanggolin_run')
        sub.run('ppanggolin msa -p ' + self.base_filename + '.h5 --phylo -o "./Phylo" -f', shell=True)
        sub.run('ppanggolin msa -p ' + self.base_filename + '.h5 -o "./msa" -f', shell=True)

