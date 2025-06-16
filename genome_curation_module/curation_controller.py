from asm_stat_handler import ASMStatParser
import os


class CurationController:
    def __init__(self):
        self.source_dir = self.get_source_dir()
        self.check_source_dir()

        curation_object = ASMStatParser(source_dir=self.source_dir)
        curation_object.get_busco_scores()
        curation_object.extract_checkM_stats()
        curation_object.curate_genomes()
        curation_object.create_curated_dir()
        curation_object.plot_genome_stats()
        return

    def get_source_dir(self):

        # Obtain parent directory from user
        # Test path: "/home/sbraganza/Documents/School/projects/pangenome"
        source_dir = input("Enter parent directory: ")

        return source_dir

    def check_source_dir(self):
        # Directory/files to check if exists.
        required_subdir = "ncbi_dataset/data"
        required_files = ["assembly_data_report.jsonl", "data_summary.tsv", "dataset_catalog.json"]
        missing_files = []

        while True:
            # Check if ~/source_dir/ncbi_dataset/data exists
            combined_path = os.path.join(self.source_dir, required_subdir) # os.path.join deals with '/' automatically

            if os.path.exists(combined_path):
                # Check if the required files exist in ~/source_dir/ncbi_dataset/data/
                for file_name in required_files:
                    file_path = os.path.join(combined_path, file_name)
                    if not os.path.isfile(file_path):
                        missing_files.append(file_name)

                if not missing_files:
                    print("Success.")
                    break
                else:
                    print(f"The required file(s) {', '.join(missing_files)} do not exist.")
            else:
                print(f"The directory '{self.source_dir}' does not contain the " +
                      "required directory tree. Please try again.")

            # Ask for source_dir until valid path with required files are obtained
            self.source_dir = self.get_source_dir()


if __name__ == '__main__':
    temp = CurationController()