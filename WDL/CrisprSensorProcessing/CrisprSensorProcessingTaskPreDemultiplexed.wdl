version development

import "BBMapDemultiplex.wdl" as demultiplex
import "CrisprSelfEditMapping.wdl" as mapping

workflow CrisprSensorPreprocessing_Workflow {
    input {
        # Pre demultiplex results
        Map[String, Array[Pair[AnnotatedSample, Array[String]]]] output_screenIdToSampleMap

        #
        #   FOR COUNTING
        #
        #NEW contains the demultiplex to screen ID map
        String? input_screenId
        Map[String, String]? input_i5ToScreenidMap
        Map[String, String]? input_barcodeToScreenidMap
        Map[String, Map[String, String]]? input_i5ToBarcodeToScreenidMap

        #NEW contains the demultiplex to sample Information
        Array[String]? input_sampleInfoVarnames
        Array[String]? input_sampleInfoVars
        Map[String, Array[String]]? input_i5ToSampleInfoVarsMap
        Map[String, Array[String]]? input_barcodeToSampleInfoVarsMap
        Map[String, Map[String, Array[String]]]? input_i5ToBarcodeToSampleInfoVarsMap

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
            input_protospacerHammingThresholdStrict=input_protospacerHammingThresholdStrict
    }

    
	# TODO: Since we have a MAP, perhaps someone can input a TSV of the sample sheet, then we can return a final table with the samples attached. Or we can return a table without the sample sheet. Make this into another workflow for preprocessing all outputs

    output {
        #
        # Guide Mapping Outputs
        #
        Map[String, Array[Pair[Pair[AnnotatedSample,Array[String]], File]]] output_screen_countResults_map = mappingWorkflow.output_screen_countResults_map
    }
}


