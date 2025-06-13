from postprocessing_module.uniprot_extractor import UniProtExtractor
from postprocessing_module.matrix_trimmer import GeneParser
from postprocessing_module.statistics_report import PostProcessing
from postprocessing_module.uniprotREST import UniProtRESTHandler
from postprocessing_module.gene_ont_proc import GeneOntology
import json


class PostprocessingController:
    def __init__(self, source_dir):
        self.source_dir = source_dir if source_dir[-1] == '/' else source_dir + '/'
        self.ppan_dir = self.source_dir + 'ppanggolin_run/'
        self.post_proc_dir = self.source_dir + 'post_processing/'
        self.grabStatReport()
        self.callUniProtExtractor()
        self.callMatrixTrimmer()
        return

    def grabStatReport(self):
        postproc_results = PostProcessing(self.source_dir)
        postproc_results.grab_pangenome_stats()
        return

    def callUniProtExtractor(self):
        uniprot_extractor = UniProtExtractor(self.source_dir)
        uniprot_extractor.parseGFF_mp()
        return

    def callMatrixTrimmer(self):
        matrix_trimmer = GeneParser(self.source_dir)
        matrix_trimmer.parseGenes_mp()

    def callUniProtREST(self):
        with open(self.post_proc_dir + "uniprot_freqs.json", "r", encoding='utf-8') as file:
            uniprots = json.load(file)

        result_handler = UniProtRESTHandler(polling_interval=30, num_retries=30)
        result_handler.submit_id_mapping(from_db="UniProtKB_AC-ID", to_db="UniProtKB",
                                         ids=uniprots.keys())

        if result_handler.check_id_mapping_results_ready():
            link = result_handler.get_id_mapping_results_link()
            results = result_handler.get_id_mapping_results_search(link)

        try:
            print(f"{len(results['results'])} results were found. {len(results['failedIds'])} failed")
        except KeyError:
            print(f"{len(results['results'])} results were found. 0 failed")

        with open(self.post_proc_dir + "go_results.json", 'w', encoding='utf8') as ofile:
            json.dump(results, ofile)


    def callGeneOntProc(self):
        ont_proc = GeneOntology(source_dir = self.source_dir)
        ont_proc.process_go_ids()
        ont_proc.go_class_splitter()
        ont_proc.plotGeneOntology()