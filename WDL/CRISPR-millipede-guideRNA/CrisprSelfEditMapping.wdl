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
            cores=9)

        
        match_set_whitelist_reporter_observed_sequence_counter_series_results = crispr_ambiguous_mapping.postprocessing.get_matchset_alleleseries(result.observed_guide_reporter_umi_counts_inferred, "protospacer_match_surrogate_match_barcode_match", contains_surrogate=True, contains_barcode=True, contains_umi=True) # TODO: HARDCODED contains_*- eventually this should be within the result object
        mutations_results = crispr_ambiguous_mapping.postprocessing.get_mutation_profile(match_set_whitelist_reporter_observed_sequence_counter_series_results, whitelist_reporter_df=whitelist_guide_reporter_df, contains_surrogate=True, contains_barcode=True) # TODO: HARDCODED - eventually this should be within the result object
        linked_mutation_counters = crispr_ambiguous_mapping.postprocessing.tally_linked_mutation_count_per_sequence(mutations_results=mutations_results, contains_surrogate = True, contains_barcode = True)# TODO: HARDCODED - eventually this should be within the result object
        crispr_ambiguous_mapping.visualization.plot_mutation_count_histogram(linked_mutation_counters.protospacer_total_mutation_counter, filename="protospacer_total_mutation_histogram.pdf")
        crispr_ambiguous_mapping.visualization.plot_mutation_count_histogram(linked_mutation_counters.surrogate_total_mutation_counter, filename="surrogate_total_mutation_histogram.pdf")
        crispr_ambiguous_mapping.visualization.plot_mutation_count_histogram(linked_mutation_counters.barcode_total_mutation_counter, filename="barcode_total_mutation_histogram.pdf")
        
        with open("protospacer_editing_efficiency.txt", "w") as text_file:
            print(crispr_ambiguous_mapping.utility.calculate_average_editing_frequency(linked_mutation_counters.protospacer_total_mutation_counter), file=text_file)
        with open("surrogate_editing_efficiency.txt", "w") as text_file:
            print(crispr_ambiguous_mapping.utility.calculate_average_editing_frequency(linked_mutation_counters.surrogate_total_mutation_counter), file=text_file)
        with open("barcode_editing_efficiency.txt", "w") as text_file:
            print(crispr_ambiguous_mapping.utility.calculate_average_editing_frequency(linked_mutation_counters.barcode_total_mutation_counter), file=text_file)
            
        crispr_ambiguous_mapping.visualization.plot_trinucleotide_mutational_signature(unlinked_mutations_df=mutations_results.ambiguous_accepted_umi_noncollapsed_mutations.all_observed_surrogate_unlinked_mutations_df, label=screen_name, filename="surrogate_trinucleotide_mutational_signature.pdf")
        crispr_ambiguous_mapping.visualization.plot_positional_mutational_signature(unlinked_mutations_df=mutations_results.ambiguous_accepted_umi_noncollapsed_mutations.all_observed_surrogate_unlinked_mutations_df, label=screen_name, min_position = 6, max_position=20, filename="surrogate_trinucleotide_positional_signature.pdf")
        
        crispr_ambiguous_mapping.utility.save_or_load_pickle("./", "match_set_whitelist_reporter_observed_sequence_counter_series_results", py_object = match_set_whitelist_reporter_observed_sequence_counter_series_results, date_string="")
        crispr_ambiguous_mapping.utility.save_or_load_pickle("./", "mutations_results", py_object = mutations_results, date_string="")
        crispr_ambiguous_mapping.utility.save_or_load_pickle("./", "linked_mutation_counters", py_object = linked_mutation_counters, date_string="")
        crispr_ambiguous_mapping.utility.save_or_load_pickle("./", "whitelist_guide_reporter_df", py_object = whitelist_guide_reporter_df, date_string="")
        crispr_ambiguous_mapping.utility.save_or_load_pickle("./", "result", py_object = result, date_string="")
        crispr_ambiguous_mapping.utility.save_or_load_pickle("./", "protospacer_match_count_result", py_object = result.all_match_set_whitelist_reporter_counter_series_results.protospacer_match, date_string="")
        crispr_ambiguous_mapping.utility.save_or_load_pickle("./", "protospacer_match_barcode_match_count_result", py_object = result.all_match_set_whitelist_reporter_counter_series_results.protospacer_match_barcode_match, date_string="")
        crispr_ambiguous_mapping.utility.save_or_load_pickle("./", "protospacer_match_surrogate_match_count_result", py_object = result.all_match_set_whitelist_reporter_counter_series_results.protospacer_match_surrogate_match, date_string="")
        crispr_ambiguous_mapping.utility.save_or_load_pickle("./", "protospacer_match_surrogate_match_barcode_match_count_result", py_object = result.all_match_set_whitelist_reporter_counter_series_results.protospacer_match_surrogate_match_barcode_match, date_string="")
        
        CODE
    >>>

    output {
        Float protospacer_editing_efficiency =  read_float("protospacer_editing_efficiency.txt")
        Float surrogate_editing_efficiency = read_float("surrogate_editing_efficiency.txt")
        Float barcode_editing_efficiency = read_float("barcode_editing_efficiency.txt")

        File match_set_whitelist_reporter_observed_sequence_counter_series_results = "match_set_whitelist_reporter_observed_sequence_counter_series_results_.pickle"
        File mutations_results = "mutations_results_.pickle"
        File linked_mutation_counters = "linked_mutation_counters_.pickle"

        File protospacer_total_mutation_histogram_pdf = "protospacer_total_mutation_histogram.pdf"
        File surrogate_total_mutation_histogram_pdf = "surrogate_total_mutation_histogram.pdf"
        File barcode_total_mutation_histogram_pdf = "barcode_total_mutation_histogram.pdf"

        File surrogate_trinucleotide_mutational_signature = "surrogate_trinucleotide_mutational_signature.pdf"
        File surrogate_trinucleotide_positional_signature = "surrogate_trinucleotide_positional_signature.pdf"

        File whitelist_guide_reporter_df = "whitelist_guide_reporter_df_.pickle"
        File count_result = "result_.pickle"
        File protospacer_match_count_result = "protospacer_match_count_result_.pickle"
        File protospacer_match_barcode_match_count_result = "protospacer_match_barcode_match_count_result_.pickle"
        File protospacer_match_surrogate_match_count_result = "protospacer_match_surrogate_match_count_result_.pickle"
        File protospacer_match_surrogate_match_barcode_match_count_result = "protospacer_match_surrogate_match_barcode_match_count_result_.pickle"
    }

    runtime {
        docker: "pinellolab/crispr_selfedit_mapping:release-0.0.130"
        memory: "32G"
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
                "protospacer_match_count_result": GuideCount_ScreenId.protospacer_match_count_result,
                "protospacer_match_barcode_match_count_result": GuideCount_ScreenId.protospacer_match_barcode_match_count_result,
                "protospacer_match_surrogate_match_count_result": GuideCount_ScreenId.protospacer_match_surrogate_match_count_result,
                "protospacer_match_surrogate_match_barcode_match_count_result": GuideCount_ScreenId.protospacer_match_surrogate_match_barcode_match_count_result
                }

                Pair[Pair[AnnotatedSample,Array[String]], File] annotated_count_result = (annotatedSamplePair, GuideCount_ScreenId.count_result)
                Pair[Pair[AnnotatedSample,Array[String]], Map[String, Float]] annotated_editing_efficiencies = (annotatedSamplePair, editing_efficiency_dict)
                Pair[Pair[AnnotatedSample,Array[String]], Map[String, File]] annotated_supplementary_files = (annotatedSamplePair, supplementary_files_dict)

            }
            Array[Pair[Pair[AnnotatedSample,Array[String]], File]] annotated_count_result_list = annotated_count_result
            Pair[String, Array[Pair[Pair[AnnotatedSample,Array[String]], File]]] screen_countResults_pair = (screenId, annotated_count_result_list)

            Array[Pair[Pair[AnnotatedSample,Array[String]], Map[String, Float]]] annotated_editing_efficiencies_list = annotated_editing_efficiencies
            Pair[String, Array[Pair[Pair[AnnotatedSample,Array[String]], Map[String, Float]]]] screen_editingEfficiencies_pair = (screenId, annotated_editing_efficiencies_list)

            Array[Pair[Pair[AnnotatedSample,Array[String]], Map[String, File]]] annotated_supplementary_files_list = annotated_supplementary_files
            Pair[String, Array[Pair[Pair[AnnotatedSample,Array[String]], Map[String, File]]]] screen_supplementaryFiles_pair = (screenId, annotated_supplementary_files_list)
        }

        

        # TODO: Perform ADATA/BDATA for each screen here! Will use the sampleInfoVars for the sample , Array[Array[String]] sampleInfoVarsScreenList = sampleInfoVars
    }

    Map[String, Array[Pair[Pair[AnnotatedSample,Array[String]], File]]] screen_countResults_map = as_map(select_all(screen_countResults_pair))
    Map[String, Array[Pair[Pair[AnnotatedSample,Array[String]], Map[String, Float]]]] screen_editingEfficiencies_map = as_map(select_all(screen_editingEfficiencies_pair))
    Map[String, Array[Pair[Pair[AnnotatedSample,Array[String]], Map[String, File]]]] screen_supplementaryFiles_map = as_map(select_all(screen_supplementaryFiles_pair))

    output {
        Map[String, Array[Pair[Pair[AnnotatedSample,Array[String]], File]]] output_screen_countResults_map = screen_countResults_map
        Map[String, Array[Pair[Pair[AnnotatedSample,Array[String]], Map[String, Float]]]] output_screen_editingEfficiencies_map = screen_editingEfficiencies_map
        Map[String, Array[Pair[Pair[AnnotatedSample,Array[String]], Map[String, File]]]] output_screen_supplementaryFiles_map = screen_supplementaryFiles_map
    }

}