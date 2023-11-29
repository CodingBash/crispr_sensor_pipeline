version development

import "BBMapDemultiplex.wdl" as demultiplex

task GuideCount {
    input {
        File countInputRead1
        File? countInputRead2
        File whitelistGuideReporterTsv
        String? umiToolsHeaderBarcodeRegex
        String? umiToolsUmiPatternRegex
        Int? surrogateHammingThresholdStrict
        Int? barcodeHammingThresholdStrict
        Int? protospacerHammingThresholdStrict
    }

    command <<<
        python <<CODE
            import crispr_ambiguous_mapping
            import pandas as pd
            
            whitelist_guide_reporter_df = pd.read_table("~{whitelistGuideReporterTsv}")

            result = crispr_ambiguous_mapping.mp.get_whitelist_reporter_counts_from_umitools_output(
                whitelist_guide_reporter_df=whitelist_guide_reporter_df, 
                fastq_r1_fn='~{countInputRead1}', 
                fastq_r2_fn=~{if defined(countInputRead2) then "'~{countInputRead2}'" else "None" },
                barcode_pattern_regex=~{if defined(umiToolsHeaderBarcodeRegex) then "'~{umiToolsHeaderBarcodeRegex}'" else "None" },
                umi_pattern_regex=~{if defined(umiToolsUmiPatternRegex) then "'~{umiToolsUmiPatternRegex}'" else "None" },
                surrogate_hamming_threshold_strict=~{if defined(surrogateHammingThresholdStrict) then "~{surrogateHammingThresholdStrict}" else "None" },
                barcode_hamming_threshold_strict =~{if defined(barcodeHammingThresholdStrict) then "~{barcodeHammingThresholdStrict}" else "None" },
                protospacer_hamming_threshold_strict=~{if defined(protospacerHammingThresholdStrict) then "~{protospacerHammingThresholdStrict}" else "None" },
                cores=1)


            crispr_ambiguous_mapping.ut.save_or_load_pickle("./", "result", py_object = result, date_string="")
        CODE
    >>>

    output {
        File count_result = "result_.pickle"
    }

    runtime {
        docker: "pinellolab/crispr_selfedit_mapping:release-0.0.106a"
    }
}

workflow CrisprSelfEditMappingOrchestratorWorkflow {
    input {
        Map[String, Array[AnnotatedSample]] input_screenIdToSampleMap

        File? input_whitelistGuideReporterTsv
        Map[String, File]? input_screenIdToWhitelistGuideReporterTsv
        Map[String, File]? input_screenIdToGuideAnnotationsTsv

        String? input_umiToolsHeaderBarcodeRegex
        String? input_umiToolsUmiPatternRegex
        
        Int? input_surrogateHammingThresholdStrict
        Int? input_barcodeHammingThresholdStrict
        Int? input_protospacerHammingThresholdStrict
    }

    scatter(input_screenIdToSamplePair in as_pairs(input_screenIdToSampleMap)){
        String screenId = input_screenIdToSamplePair.left
        Array[AnnotatedSample] screenAnnotatedSamples = input_screenIdToSamplePair.right

        scatter(annotatedSample in screenAnnotatedSamples){

            # Select the whitelist guide reporter TSV to use for mapping
            if (defined(input_screenIdToWhitelistGuideReporterTsv)){
                Map[String, File] input_screenIdToWhitelistGuideReporterTsv_defined = select_first([input_screenIdToWhitelistGuideReporterTsv])
                File screen_whitelistGuideReporterTsv = input_screenIdToWhitelistGuideReporterTsv_defined[screenId]
            }
            File selected_whitelistGuideReporterTsv = select_first([screen_whitelistGuideReporterTsv, input_whitelistGuideReporterTsv])
            
            #
            #   Perform guide mapping of sample
            #
            call GuideCount as GuideCount_ScreenId {
                input:
                    countInputRead1=annotatedSample.read1,
                    countInputRead2=annotatedSample.read2,
                    whitelistGuideReporterTsv=selected_whitelistGuideReporterTsv,
                    umiToolsHeaderBarcodeRegex=input_umiToolsHeaderBarcodeRegex,
                    umiToolsUmiPatternRegex=input_umiToolsUmiPatternRegex,
                    surrogateHammingThresholdStrict=input_surrogateHammingThresholdStrict,
                    barcodeHammingThresholdStrict=input_barcodeHammingThresholdStrict,
                    protospacerHammingThresholdStrict=input_protospacerHammingThresholdStrict
            }
        }

        Array[File] screen_countResults = GuideCount_ScreenId.count_result
        Pair[String, Array[File]] screen_countResults_pair = (screenId, screen_countResults)
    }

    Map[String, Array[File]] screen_countResults_map = as_map(screen_countResults_pair)

    output {
        Map[String, Array[File]] output_screen_countResults_map = screen_countResults_map
    }

}