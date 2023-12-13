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
            barcode_pattern_regex=~{if defined(umiToolsHeaderBarcodeRegex) then "~{umiToolsHeaderBarcodeRegex}" else "None" },
            umi_pattern_regex=~{if defined(umiToolsUmiPatternRegex) then "~{umiToolsUmiPatternRegex}" else "None" },
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
        docker: "pinellolab/crispr_selfedit_mapping:release-0.0.108a"
         memory: "4G"
    }
}

workflow CrisprSelfEditMappingOrchestratorWorkflow {
    input {
        Map[String, Array[Pair[AnnotatedSample, Array[String]]]] input_screenIdToSampleMap

        File? input_whitelistGuideReporterTsv
        Map[String, File]? input_screenIdToWhitelistGuideReporterTsv
        Map[String, File]? input_screenIdToGuideAnnotationsTsv

        String? input_umiToolsHeaderBarcodeRegex
        String? input_umiToolsUmiPatternRegex
        
        Int? input_surrogateHammingThresholdStrict
        Int? input_barcodeHammingThresholdStrict
        Int? input_protospacerHammingThresholdStrict
    }

    #
    #   IMPORTANT NOTE: If input_whitelistGuideReporterTsv is provided but not input_screenIdToWhitelistGuideReporterTsv (i.e. same guide table for all samples), then workflow will fail. Will end up with duplicate code, but use if statement to do happy and exception path.
    #
    # Iterate through all screenId-sample pairs
    scatter(input_screenIdToSamplePair in as_pairs(input_screenIdToSampleMap)){
        String screenId = input_screenIdToSamplePair.left
        Array[Pair[AnnotatedSample,Array[String]]] screenAnnotatedSamples = input_screenIdToSamplePair.right
        
        if(screenId != "None"){
            # Iterate through each saple of the screen ID
            scatter(annotatedSamplePair in screenAnnotatedSamples){
                AnnotatedSample annotatedSample = annotatedSamplePair.left
                Array[String] sampleInfoVars = annotatedSamplePair.right

                # Select the whitelist guide reporter TSV to use for mapping
                Map[String, File] input_screenIdToWhitelistGuideReporterTsv_defined = select_first([input_screenIdToWhitelistGuideReporterTsv])
                File screen_whitelistGuideReporterTsv = input_screenIdToWhitelistGuideReporterTsv_defined[screenId]
                
                #
                #   Perform guide mapping of sample
                #
                call GuideCount as GuideCount_ScreenId {
                    input:
                        countInputRead1=annotatedSample.read1,
                        countInputRead2=annotatedSample.read2,
                        whitelistGuideReporterTsv=screen_whitelistGuideReporterTsv,
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

        

        # TODO: Perform ADATA/BDATA for each screen here! Will use the sampleInfoVars for the sample , Array[Array[String]] sampleInfoVarsScreenList = sampleInfoVars
    }

    Map[String, Array[File]] screen_countResults_map = as_map(select_all(screen_countResults_pair))

    output {
        Map[String, Array[File]] output_screen_countResults_map = screen_countResults_map
    }

}