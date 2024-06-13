version development

import "https://api.firecloud.org/ga4gh/v1/tools/pinellolab:CrisprSensorGuideCountTask/versions/3/plain-WDL/descriptor" as count

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
        String dockerImage = "pinellolab/crispr_selfedit_mapping:release-0.0.140"
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


    Map[String, Float] editing_efficiency_dict = {"protospacer_editing_efficiency":GuideCount_ScreenId.protospacer_editing_efficiency,
    "surrogate_editing_efficiency": GuideCount_ScreenId.surrogate_editing_efficiency,
    "barcode_editing_efficiency": GuideCount_ScreenId.barcode_editing_efficiency}

    Map[String, File] supplementary_files_dict = {"match_set_whitelist_reporter_observed_sequence_counter_series_results": GuideCount_ScreenId.match_set_whitelist_reporter_observed_sequence_counter_series_results,
    "mutations_results": GuideCount_ScreenId.mutations_results,
    "linked_mutation_counters": GuideCount_ScreenId.linked_mutation_counters,
    "protospacer_total_mutation_histogram_pdf": GuideCount_ScreenId.protospacer_total_mutation_histogram_pdf,
    "surrogate_total_mutation_histogram_pdf": GuideCount_ScreenId.surrogate_total_mutation_histogram_pdf,
    "barcode_total_mutation_histogram_pdf": GuideCount_ScreenId.barcode_total_mutation_histogram_pdf,
    "surrogate_trinucleotide_mutational_signature": GuideCount_ScreenId.surrogate_trinucleotide_mutational_signature,
    "surrogate_trinucleotide_positional_signature": GuideCount_ScreenId.surrogate_trinucleotide_positional_signature,
    "whitelist_guide_reporter_df": GuideCount_ScreenId.whitelist_guide_reporter_df,
    "count_series_result": GuideCount_ScreenId.count_series_result,
    "observed_guide_reporter_umi_counts_inferred": GuideCount_ScreenId.observed_guide_reporter_umi_counts_inferred,
    "quality_control_result": GuideCount_ScreenId.quality_control_result,
    "count_input": GuideCount_ScreenId.count_input
    }

    File count_result = GuideCount_ScreenId.count_result

    output {
        Map[String, Float] output_editing_efficiency_dict = editing_efficiency_dict
        Map[String, File] output_supplementary_files_dict = supplementary_files_dict
        File output_count_result = count_result
    }

}