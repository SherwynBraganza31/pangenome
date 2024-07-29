import subprocess as sub
import os

def runCommand():
    # pangenome_path = "/home/sbraganza/projects/10genomes"  # temp
    # cpu = "20"                                             # temp
    base_filename = "tengenomes"

    # Obtain from user
    pangenome_path = input("Enter parent directory: ")
    cpu = input("Enter cpu: ")

    # Check for files before running
    try:
        os.chdir(pangenome_path + '\ppan_dataset')
    except FileNotFoundError as e:
        log_error(f"Error changing directory: {e}")
        return


    def base_call():
        # Change directory & create a clustered .h5 file
        os.chdir(pangenome_path + '/ppan_dataset')
        sub.run('ppanggolin annotate --anno organisms.gff.list --fasta organisms.fasta.list --output ' +
                       pangenome_path + '/ppanggolin_run --basename ' + base_filename + ' -f --cpu ' + cpu, shell=True)
        os.chdir(pangenome_path + '/ppanggolin_run')
        sub.run('ppanggolin cluster -p ' + base_filename + '.h5 --identity 0.5 --coverage 0.8 -f --cpu ' + cpu, shell=True)


    def pangenome_data_gen():
        # Generate Pangenome Data
        os.chdir(pangenome_path + '/ppanggolin_run')
        sub.run('ppanggolin graph -p ' + base_filename + '.h5 -f', shell=True)
        sub.run('ppanggolin partition -p ' + base_filename + '.h5 -f', shell=True)
        sub.run('ppanggolin info -p ' + base_filename + '.h5 --content', shell=True)


    def curves_matrix():
        # Generate Curves and Matrix Files
        os.chdir(pangenome_path + '/ppanggolin_run')
        sub.run('ppanggolin draw -p ' + base_filename + '.h5 --ucurve --tile_plot --output ' + base_filename + '_withcloud -f', shell=True)
        sub.run('ppanggolin draw -p ' + base_filename + '.h5 --tile_plot --nocloud --output ' + base_filename + '_nocloud -f', shell=True)
        sub.run('ppanggolin write -p ' + base_filename + '.h5 --gexf --csv --output matrix -f', shell=True)
        sub.run('ppanggolin write -p ' + base_filename + '.h5 --light_gexf --output light_gexf -f', shell=True)


    def rgp_gen():
        # Generate Regions of genome plasticity
        os.chdir(pangenome_path + '/ppanggolin_run')
        sub.run('ppanggolin rgp -p ' + base_filename + '.h5', shell=True)
        sub.run('ppanggolin write -p ' + base_filename + '.h5 --regions --output RGP-Regions -f', shell=True)


    def spots_of_insertion():
        os.chdir(pangenome_path + '/ppanggolin_run')
        # Generate Spot Figures
        sub.run('ppanggolin spot -p ' + base_filename + '.h5 -f', shell=True)
        sub.run('ppanggolin write -p ' + base_filename + '.h5 --spots --output spots -f', shell=True)
        sub.run('ppanggolin draw -p ' + base_filename + '.h5 --spots all --output spot-figures -f', shell=True)


    def rarefaction_curves():
        os.chdir(pangenome_path + '/ppanggolin_run')
        # Generate Regions of genome plasticity
        sub.run('ppanggolin rarefaction -p ' + base_filename + '.h5 --output rarefaction -f', shell=True)


    def partitions_projections_families():
        os.chdir(pangenome_path + '/ppanggolin_run')
        # Generate Partitions, Projections, Gene Families
        sub.run('ppanggolin write -p ' + base_filename + '.h5 --partitions --output partitions -f', shell=True)
        sub.run('ppanggolin write -p ' + base_filename + '.h5 --projection --output projection -f', shell=True)
        sub.run('ppanggolin write -p ' + base_filename + '.h5 --families_tsv --output families-tsv -f', shell=True)


    def fasta_outputs():
        os.chdir(pangenome_path + '/ppanggolin_run')
        # Output FASTA files of Gene Data
        sub.run('ppanggolin fasta -p ' + base_filename + '.h5 --output MY_GENES --genes all -f', shell=True)
        sub.run('ppanggolin fasta -p ' + base_filename + '.h5 --output MY_GENES/persistent --genes persistent -f', shell=True)
        sub.run('ppanggolin fasta -p ' + base_filename + '.h5 --output MY_GENES/shell --genes shell -f', shell=True)
        sub.run('ppanggolin fasta -p ' + base_filename + '.h5 --output MY_GENES/cloud --genes cloud -f', shell=True)
        sub.run('ppanggolin fasta -p ' + base_filename + '.h5 --output MY_GENES/core --genes core -f', shell=True)
        sub.run('ppanggolin fasta -p ' + base_filename + '.h5 --output MY_GENES/softcore --genes softcore -f', shell=True)
        sub.run('ppanggolin fasta -p ' + base_filename + '.h5 --output MY_GENES/rgp --genes rgp -f', shell=True)

        # Output FASTA files of Protein Families
        sub.run('ppanggolin fasta -p ' + base_filename + '.h5 --output MY_PROT --prot_families all -f', shell=True)
        sub.run('ppanggolin fasta -p ' + base_filename + '.h5 --output MY_PROT/persistent --prot_families persistent -f', shell=True)
        sub.run('ppanggolin fasta -p ' + base_filename + '.h5 --output MY_PROT/shell --prot_families shell -f', shell=True)
        sub.run('ppanggolin fasta -p ' + base_filename + '.h5 --output MY_PROT/cloud --prot_families cloud -f', shell=True)
        sub.run('ppanggolin fasta -p ' + base_filename + '.h5 --output MY_PROT/core --prot_families core -f', shell=True)
        sub.run('ppanggolin fasta -p ' + base_filename + '.h5 --output MY_PROT/softcore --prot_families softcore -f', shell=True)
        sub.run('ppanggolin fasta -p ' + base_filename + '.h5 --output MY_PROT/rgp --prot_families rgp -f', shell=True)


    def phylogeny():
        os.chdir(pangenome_path + '/ppanggolin_run')
        sub.run('ppanggolin msa -p ' + base_filename + '.h5 --phylo -o "./Phylo" -f', shell=True)
        sub.run('ppanggolin msa -p ' + base_filename + '.h5 -o "./msa" -f', shell=True)


    # Call Functions
    base_call()
    pangenome_data_gen()
    curves_matrix()
    rgp_gen()
    spots_of_insertion()
    rarefaction_curves()
    partitions_projections_families()
    fasta_outputs()
    phylogeny()


if __name__ == "__main__":
    runCommand()