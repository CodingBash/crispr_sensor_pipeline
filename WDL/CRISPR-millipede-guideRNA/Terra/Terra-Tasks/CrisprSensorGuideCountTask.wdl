version development

task CrisprSensorGuideCountTask {
    input {
        String screenId
        File countInputRead1
        File? countInputRead2
        File whitelistGuideReporterTsv
        String? umiToolsHeaderBarcodeRegex
        String? umiToolsUmiPatternRegex
        Int? surrogateHammingThresholdStrict
        Int? barcodeHammingThresholdStrict
        Int? protospacerHammingThresholdStrict

        String dockerImage = "pinellolab/crispr_selfedit_mapping:release-0.0.142"
        Int preemptible = 1
        Int diskGB = 10
        Int memoryGB = 2
        Int maxRetries = 0
        String diskType = "HDD"
        Int cpus = 1
    }

    command <<<
        python <<CODE

        import crispr_ambiguous_mapping
        import pandas as pd
        
        whitelist_guide_reporter_df = pd.read_table("~{whitelistGuideReporterTsv}")

        result = crispr_ambiguous_mapping.mapping.get_whitelist_reporter_counts_from_umitools_output(
            whitelist_guide_reporter_df=whitelist_guide_reporter_df, 
            fastq_r1_fn='~{countInputRead1}', 
            fastq_r2_fn=~{if defined(countInputRead2) then "'~{countInputRead2}'" else "None" },
            barcode_pattern_regex=~{if defined(umiToolsHeaderBarcodeRegex) then "~{umiToolsHeaderBarcodeRegex}" else "None" },
            umi_pattern_regex=~{if defined(umiToolsUmiPatternRegex) then "~{umiToolsUmiPatternRegex}" else "None" },
            surrogate_hamming_threshold_strict=~{if defined(surrogateHammingThresholdStrict) then "~{surrogateHammingThresholdStrict}" else "None" },
            barcode_hamming_threshold_strict =~{if defined(barcodeHammingThresholdStrict) then "~{barcodeHammingThresholdStrict}" else "None" },
            protospacer_hamming_threshold_strict=~{if defined(protospacerHammingThresholdStrict) then "~{protospacerHammingThresholdStrict}" else "None" },
            cores=~{cpus})

        
        match_set_whitelist_reporter_observed_sequence_counter_series_results = crispr_ambiguous_mapping.processing.get_matchset_alleleseries(result.observed_guide_reporter_umi_counts_inferred, "protospacer_match_surrogate_match_barcode_match", contains_surrogate=True, contains_barcode=True, contains_umi=True) # TODO: HARDCODED contains_*- eventually this should be within the result object
        mutations_results = crispr_ambiguous_mapping.processing.get_mutation_profile(match_set_whitelist_reporter_observed_sequence_counter_series_results, whitelist_reporter_df=whitelist_guide_reporter_df, contains_surrogate=True, contains_barcode=True) # TODO: HARDCODED - eventually this should be within the result object
        linked_mutation_counters = crispr_ambiguous_mapping.processing.tally_linked_mutation_count_per_sequence(mutations_results=mutations_results, contains_surrogate = True, contains_barcode = True)# TODO: HARDCODED - eventually this should be within the result object
        crispr_ambiguous_mapping.visualization.plot_mutation_count_histogram(linked_mutation_counters.protospacer_total_mutation_counter, filename="protospacer_total_mutation_histogram.png")
        crispr_ambiguous_mapping.visualization.plot_mutation_count_histogram(linked_mutation_counters.surrogate_total_mutation_counter, filename="surrogate_total_mutation_histogram.png")
        crispr_ambiguous_mapping.visualization.plot_mutation_count_histogram(linked_mutation_counters.barcode_total_mutation_counter, filename="barcode_total_mutation_histogram.png")
        
        with open("protospacer_editing_efficiency.txt", "w") as text_file:
            print(crispr_ambiguous_mapping.utility.calculate_average_editing_frequency(linked_mutation_counters.protospacer_total_mutation_counter), file=text_file)
        with open("surrogate_editing_efficiency.txt", "w") as text_file:
            print(crispr_ambiguous_mapping.utility.calculate_average_editing_frequency(linked_mutation_counters.surrogate_total_mutation_counter), file=text_file)
        with open("barcode_editing_efficiency.txt", "w") as text_file:
            print(crispr_ambiguous_mapping.utility.calculate_average_editing_frequency(linked_mutation_counters.barcode_total_mutation_counter), file=text_file)
            
        crispr_ambiguous_mapping.visualization.plot_trinucleotide_mutational_signature(mutations_results=mutations_results, count_attribute_name="ambiguous_accepted_umi_noncollapsed_mutations", unlinked_mutation_attribute_name = "all_observed_surrogate_unlinked_mutations_df", label='~{screenId}', filename="surrogate_trinucleotide_mutational_signature.png")
        crispr_ambiguous_mapping.visualization.plot_positional_mutational_signature(mutations_results=mutations_results, count_attribute_name="ambiguous_accepted_umi_noncollapsed_mutations", unlinked_mutation_attribute_name = "all_observed_surrogate_unlinked_mutations_df", label='~{screenId}', min_position = 6, max_position=20, filename="surrogate_trinucleotide_positional_signature.png")
        
        crispr_ambiguous_mapping.utility.save_or_load_pickle("./", "match_set_whitelist_reporter_observed_sequence_counter_series_results", py_object = match_set_whitelist_reporter_observed_sequence_counter_series_results, date_string="")
        crispr_ambiguous_mapping.utility.save_or_load_pickle("./", "mutations_results", py_object = mutations_results, date_string="")
        crispr_ambiguous_mapping.utility.save_or_load_pickle("./", "linked_mutation_counters", py_object = linked_mutation_counters, date_string="")
        crispr_ambiguous_mapping.utility.save_or_load_pickle("./", "whitelist_guide_reporter_df", py_object = whitelist_guide_reporter_df, date_string="")
        crispr_ambiguous_mapping.utility.save_or_load_pickle("./", "result", py_object = result, date_string="")
        crispr_ambiguous_mapping.utility.save_or_load_pickle("./", "count_series_result", py_object = result.all_match_set_whitelist_reporter_counter_series_results, date_string="")
        
        CODE
    >>>

    output {
        Float protospacer_editing_efficiency =  read_float("protospacer_editing_efficiency.txt")
        Float surrogate_editing_efficiency = read_float("surrogate_editing_efficiency.txt")
        Float barcode_editing_efficiency = read_float("barcode_editing_efficiency.txt")

        File match_set_whitelist_reporter_observed_sequence_counter_series_results = "match_set_whitelist_reporter_observed_sequence_counter_series_results_.pickle"
        File mutations_results = "mutations_results_.pickle"
        File linked_mutation_counters = "linked_mutation_counters_.pickle"

        File protospacer_total_mutation_histogram_pdf = "protospacer_total_mutation_histogram.png"
        File surrogate_total_mutation_histogram_pdf = "surrogate_total_mutation_histogram.png"
        File barcode_total_mutation_histogram_pdf = "barcode_total_mutation_histogram.png"

        File surrogate_trinucleotide_mutational_signature = "surrogate_trinucleotide_mutational_signature.png"
        File surrogate_trinucleotide_positional_signature = "surrogate_trinucleotide_positional_signature.png"

        File whitelist_guide_reporter_df = "whitelist_guide_reporter_df_.pickle"
        File count_result = "result_.pickle"
        File count_series_result = "count_series_result_.pickle"
    }

    runtime {
        docker: "${dockerImage}"
        preemptible: "${preemptible}"
        maxRetries: "${maxRetries}"
        memory: "${memoryGB} GB"
        disks: "local-disk ${diskGB} ${diskType}"
        cpu: "${cpus}"
    }
}