version development

import "https://api.firecloud.org/ga4gh/v1/tools/pinellolab:BBMapDemultiplexOrchestratorWorkflow/versions/6/plain-WDL/descriptor" as demultiplex
import "https://api.firecloud.org/ga4gh/v1/tools/pinellolab:CrisprSelfEditMappingOrchestratorWorkflow/versions/4/plain-WDL/descriptor" as mapping

workflow CrisprSensorPreprocessing_Workflow {
    input {
        # DEMULTIPLEX INPUT SAMPLES
        Map[String, Array[Pair[AnnotatedSample, Array[String]]]] output_screenIdToSampleMap

        #
        #   FOR COUNTING
        #
        # NEW contains the screenId to reporter TSV map
        File? input_whitelistGuideReporterTsv
        Map[String, File]? input_screenIdToWhitelistGuideReporterTsv
        Map[String, File]? input_screenIdToGuideAnnotationsTsv # NOTE: This can be based on the IGVF guide specification or any flexible schema

        # NEW contains annotations of demultiplex to annotations table
        String? input_umiToolsHeaderBarcodeRegex
        String? input_umiToolsUmiPatternRegex
        
        Int? input_surrogateHammingThresholdStrict
        Int? input_barcodeHammingThresholdStrict
        Int? input_protospacerHammingThresholdStrict

        # RUNTIME PARAMS
        String guideMappingTaskDockerImage = "pinellolab/crispr_selfedit_mapping:release-0.0.140"
        Int guideMappingTaskPreemptible = 1
        Int guideMappingTaskDiskGB = 10
        Int guideMappingTaskMemoryGB = 2
        Int guideMappingTaskMaxRetries = 0
        String guideMappingTaskDiskType = "HDD"
        Int guideMappingTaskCpus = 1
    }

    call mapping.CrisprSelfEditMappingOrchestratorWorkflow as mappingWorkflow {
        input:
            input_screenIdToSampleMap=output_screenIdToSampleMap,

            input_whitelistGuideReporterTsv=input_whitelistGuideReporterTsv,
            input_screenIdToWhitelistGuideReporterTsv=input_screenIdToWhitelistGuideReporterTsv,
            input_screenIdToGuideAnnotationsTsv=input_screenIdToGuideAnnotationsTsv,

            input_umiToolsHeaderBarcodeRegex=input_umiToolsHeaderBarcodeRegex,
            input_umiToolsUmiPatternRegex=input_umiToolsUmiPatternRegex,
        
            input_surrogateHammingThresholdStrict=input_surrogateHammingThresholdStrict,
            input_barcodeHammingThresholdStrict=input_barcodeHammingThresholdStrict,
            input_protospacerHammingThresholdStrict=input_protospacerHammingThresholdStrict,

            dockerImage=guideMappingTaskDockerImage,
            preemptible=guideMappingTaskPreemptible,
            diskGB=guideMappingTaskDiskGB,
            memoryGB=guideMappingTaskMemoryGB,
            maxRetries=guideMappingTaskMaxRetries,
            diskType=guideMappingTaskDiskType,
            cpus=guideMappingTaskCpus
    }

    
	# TODO: Since we have a MAP, perhaps someone can input a TSV of the sample sheet, then we can return a final table with the samples attached. Or we can return a table without the sample sheet. Make this into another workflow for preprocessing all outputs

    output {
        #
        # Guide Mapping Outputs
        #
        Map[String, Array[Pair[Pair[AnnotatedSample,Array[String]], File]]] output_screen_countResults_map = mappingWorkflow.output_screen_countResults_map
        Map[String, Array[Pair[Pair[AnnotatedSample,Array[String]], Map[String, Float]]]] output_screen_editingEfficiencies_map = mappingWorkflow.output_screen_editingEfficiencies_map
        Map[String, Array[Pair[Pair[AnnotatedSample,Array[String]], Map[String, File]]]] output_screen_supplementaryFiles_map = mappingWorkflow.output_screen_supplementaryFiles_map
    }
}


