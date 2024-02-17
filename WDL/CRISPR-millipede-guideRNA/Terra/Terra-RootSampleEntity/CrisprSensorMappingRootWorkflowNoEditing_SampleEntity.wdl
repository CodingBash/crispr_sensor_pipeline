version development

import "https://api.firecloud.org/ga4gh/v1/tools/pinellolab:CrisprSensorGuideCountTaskNoEditing/versions/1/plain-WDL/descriptor" as count

workflow CrisprSelfEditMappingOrchestratorWorkflowSampleEntity {
    input {
        # TASK PARAMS
        File umiToolsRead1
        File? umiToolsRead2

        File input_whitelistGuideReporterTsv

        String sampleName
        String? input_umiToolsHeaderBarcodeRegex
        String? input_umiToolsUmiPatternRegex
        
        Int? input_surrogateHammingThresholdStrict
        Int? input_barcodeHammingThresholdStrict
        Int? input_protospacerHammingThresholdStrict


        # RUNTIME PARAMS
        String dockerImage = "pinellolab/crispr_selfedit_mapping:release-0.0.142"
        Int preemptible = 1
        Int diskGB = 10
        Int memoryGB = 2
        Int maxRetries = 0
        String diskType = "HDD"
        Int cpus = 1
    }

    call count.CrisprSensorGuideCountTask as GuideCount_ScreenId {
        input:
            screenId=sampleName,
            countInputRead1=umiToolsRead1,
            countInputRead2=umiToolsRead2,
            whitelistGuideReporterTsv=input_whitelistGuideReporterTsv,
            umiToolsHeaderBarcodeRegex=input_umiToolsHeaderBarcodeRegex,
            umiToolsUmiPatternRegex=input_umiToolsUmiPatternRegex,
            surrogateHammingThresholdStrict=input_surrogateHammingThresholdStrict,
            barcodeHammingThresholdStrict=input_barcodeHammingThresholdStrict,
            protospacerHammingThresholdStrict=input_protospacerHammingThresholdStrict,
            dockerImage=dockerImage,
            preemptible=preemptible,
            diskGB=diskGB,
            memoryGB=memoryGB,
            maxRetries=maxRetries,
            diskType=diskType,
            cpus=cpus
    }

    Map[String, File] supplementary_files_dict = {
        "whitelist_guide_reporter_df": GuideCount_ScreenId.whitelist_guide_reporter_df,
        "count_series_result": GuideCount_ScreenId.count_series_result
    }

    File count_result = GuideCount_ScreenId.count_result

    output {
        Map[String, File] output_supplementary_files_dict = supplementary_files_dict
        File output_count_result = count_result
    }

}