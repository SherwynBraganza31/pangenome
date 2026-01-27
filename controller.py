import sys
import os

from annotation_module.annotation_controller import  AnnotationController
from genome_curation_module.curation_controller import CurationController
from ppanggolin_module.ppan_controller import PpanController
from postprocessing_module.postprocessing_controller import PostprocessingController


if __name__ == '__main__':
    curation_controller = CurationController()
    annotation_controller = AnnotationController(source_dir=curation_controller.source_dir)
    ppan_controller = PpanController(source_dir=curation_controller.source_dir)
    postprocessing_controller = PostprocessingController(source_dir=curation_controller.source_dir)

