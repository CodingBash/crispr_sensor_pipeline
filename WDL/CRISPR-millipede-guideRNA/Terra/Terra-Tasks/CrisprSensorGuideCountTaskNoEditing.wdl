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

        
        crispr_ambiguous_mapping.utility.save_or_load_pickle("./", "whitelist_guide_reporter_df", py_object = whitelist_guide_reporter_df, date_string="")
        crispr_ambiguous_mapping.utility.save_or_load_pickle("./", "result", py_object = result, date_string="")
        crispr_ambiguous_mapping.utility.save_or_load_pickle("./", "count_series_result", py_object = result.all_match_set_whitelist_reporter_counter_series_results, date_string="")
        
        CODE
    >>>

    output {
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