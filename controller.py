import sys
import os

# Import Files
import genome_stat_parser.asm_stat_handler
import genome_stat_parser.busco_handler


if __name__ == '__main__':
    # Obtain parent directory from user
    pangenome_path = input("Enter parent directory: ")
    # pangenome_path = "/home/sbraganza/projects/10genomes"  # temp

    # pass into other modules
    reduced_dataset = pangenome_path + '/curated_dataset'

    # asm_stat_handler.py
    # temp = ASMStatParser(source_dir='/home/sbraganza/projects/171genomes/ncbi_dataset/data/')
    temp = ASMStatParser(source_dir = pangenome_path + '/ncbi_dataset/data/')
    temp.extract_checkM_stats()
    temp.curate_genomes()
    temp.create_curated_dir()

    # busco_handler.py
    # busco_handler = BuscoHandler(source_dir='/home/sbraganza/projects/171genomes/ncbi_dataset/data')
    busco_handler = BuscoHandler(source_dir = pangenome_path + '/ncbi_dataset/data/')
    busco_handler.run_BUSCO()
    busco_handler.update_BUSCO_scores()

    # annotation_module
    #  -> fasta_folder_gen.py
    #  -> prokka_handler.py