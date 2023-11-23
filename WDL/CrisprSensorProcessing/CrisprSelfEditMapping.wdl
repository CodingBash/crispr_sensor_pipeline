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
        Pair[DemultiplexedFiles, UndeterminedFiles]? input_DemultiplexedResult_i5
        Map[String, Pair[IndexPair, Pair[DemultiplexedFiles, UndeterminedFiles]]]? input_readIndexMap_i5_Barcode_Map
        Pair[DemultiplexedFiles, UndeterminedFiles]? input_DemultiplexedResult_Barcode
        File input_nonDemultiplexedRead1
        File? input_nonDemultiplexedRead2

        File? input_whitelistGuideReporterTsv
        Map[String, File]? input_i5Only_whitelistGuideReporterTsv
        Map[String, File]? input_barcodeOnly_whitelistGuideReporterTsv
        Map[String, Map[String, File]]? input_i5Barcode_whitelistGuideReporterTsv

        String? input_umiToolsHeaderBarcodeRegex
        String? input_umiToolsUmiPatternRegex
        
        Int? input_surrogateHammingThresholdStrict
        Int? input_barcodeHammingThresholdStrict
        Int? input_protospacerHammingThresholdStrict
    }

    Boolean i5_demultiplexed = defined(input_DemultiplexedResult_i5) && !defined(input_readIndexMap_i5_Barcode_Map)
    Boolean i5_Barcode_demultiplexed = defined(input_readIndexMap_i5_Barcode_Map)
    Boolean barcode_demultiplexed = defined(input_DemultiplexedResult_Barcode)
    Boolean not_demultiplexed = !i5_demultiplexed && !i5_Barcode_demultiplexed && !barcode_demultiplexed

    #
    #   Only demultiplexed by i5
    #
    if (i5_demultiplexed){
        Pair[DemultiplexedFiles, UndeterminedFiles] output_DemultiplexedResult_i5_defined = select_first([input_DemultiplexedResult_i5])
        
        #
        #   Scatter through the i5 indices
        #
        scatter(demultiplexedFiles_i5 in as_pairs(output_DemultiplexedResult_i5_defined.left.demultiplexedFiles)){ 
            
            String demultiplexedFiles_i5_Index = demultiplexedFiles_i5.left
            IndexPair demultiplexedFiles_i5_IndexPair = demultiplexedFiles_i5.right
            
            #
            #   Choose the whitelist_guide_reporter_tsv
            #
            if(defined(input_i5Only_whitelistGuideReporterTsv)){
                Map[String, File] input_i5Only_whitelistGuideReporterTsv_defined = select_first([input_i5Only_whitelistGuideReporterTsv])
                File input_i5Only_whitelistGuideReporterTsv_value = input_i5Only_whitelistGuideReporterTsv_defined[demultiplexedFiles_i5_Index]
            }
            File whitelistGuideReporterTsv_i5 = select_first([input_i5Only_whitelistGuideReporterTsv_value, input_whitelistGuideReporterTsv])
            
            #
            #   Perform guide mapping
            #
            call GuideCount as GuideCount_i5 {
                input:
                    countInputRead1=demultiplexedFiles_i5_IndexPair.read1,
                    countInputRead2=demultiplexedFiles_i5_IndexPair.read2,
                    whitelistGuideReporterTsv=whitelistGuideReporterTsv_i5,
                    umiToolsHeaderBarcodeRegex=input_umiToolsHeaderBarcodeRegex,
                    umiToolsUmiPatternRegex=input_umiToolsUmiPatternRegex,
                    surrogateHammingThresholdStrict=input_surrogateHammingThresholdStrict,
                    barcodeHammingThresholdStrict=input_barcodeHammingThresholdStrict,
                    protospacerHammingThresholdStrict=input_protospacerHammingThresholdStrict
            }

            Pair[String, File] GuideCount_i5_count_result_pair = (demultiplexedFiles_i5_Index, GuideCount_i5.count_result)
        }

        Array[Pair[String, File]] GuideCount_i5_count_result_pair_list = GuideCount_i5_count_result_pair
        Map[String, File] GuideCount_i5_count_result_map = as_map(select_all(GuideCount_i5_count_result_pair_list))
    }
    

    #
    #   Demultiplexed by i5 and barcode
    #
    if (i5_Barcode_demultiplexed){
        Map[String, Pair[IndexPair, Pair[DemultiplexedFiles, UndeterminedFiles]]] output_readIndexMap_i5_Barcode_Map_defined = select_first([input_readIndexMap_i5_Barcode_Map])
        
        #
        #   Scatter through the i5 indices
        #
        scatter(output_readIndexMap_i5_Barcode_Map_defined_i5Pair in as_pairs(output_readIndexMap_i5_Barcode_Map_defined)) { 
            String output_readIndexMap_i5_Barcode_Map_defined_i5Index = output_readIndexMap_i5_Barcode_Map_defined_i5Pair.left
            Pair[DemultiplexedFiles, UndeterminedFiles] output_readIndexMap_i5_Barcode_Map_defined_i5_barcodeFiles = output_readIndexMap_i5_Barcode_Map_defined_i5Pair.right.right


            #
            #   Scatter through the barcode indices
            #
            scatter(demultiplexedFiles_i5_Barcode in as_pairs(output_readIndexMap_i5_Barcode_Map_defined_i5_barcodeFiles.left.demultiplexedFiles)){ # ITERATE THROUGH i5 INDEXES
            
                String demultiplexedFiles_i5_BarcodeIndex = demultiplexedFiles_i5_Barcode.left
                IndexPair demultiplexedFiles_i5_Barcode_IndexPair = demultiplexedFiles_i5_Barcode.right
                
                #
                #   Choose the whitelist_guide_reporter_tsv
                #
                if(defined(input_i5Barcode_whitelistGuideReporterTsv)){
                    Map[String, Map[String, File]] input_i5Barcode_whitelistGuideReporterTsv_defined = select_first([input_i5Barcode_whitelistGuideReporterTsv])
                    File input_i5Barcode_whitelistGuideReporterTsv_value = input_i5Barcode_whitelistGuideReporterTsv_defined[output_readIndexMap_i5_Barcode_Map_defined_i5Index][demultiplexedFiles_i5_BarcodeIndex]
                }
                if(defined(input_i5Only_whitelistGuideReporterTsv)){
                    Map[String, File] input_i5Only_whitelistGuideReporterTsv_defined = select_first([input_i5Only_whitelistGuideReporterTsv])
                    File input_i5Only_whitelistGuideReporterTsv_value = input_i5Only_whitelistGuideReporterTsv_defined[output_readIndexMap_i5_Barcode_Map_defined_i5Index]
                }
                if(defined(input_barcodeOnly_whitelistGuideReporterTsv)){
                    Map[String, File] input_barcodeOnly_whitelistGuideReporterTsv_defined = select_first([input_barcodeOnly_whitelistGuideReporterTsv])
                    File input_barcodeOnly_whitelistGuideReporterTsv_value = input_barcodeOnly_whitelistGuideReporterTsv_defined[demultiplexedFiles_i5_BarcodeIndex]
                }
                File whitelistGuideReporterTsv_i5_Barcode = select_first([input_i5Barcode_whitelistGuideReporterTsv_value, input_i5Only_whitelistGuideReporterTsv_value, input_barcodeOnly_whitelistGuideReporterTsv_value, input_whitelistGuideReporterTsv])
                
                #
                #   Perform guide mapping of sample
                #
                call GuideCount as GuideCount_i5_Barcode {
                    input:
                        countInputRead1=demultiplexedFiles_i5_Barcode_IndexPair.read1,
                        countInputRead2=demultiplexedFiles_i5_Barcode_IndexPair.read2,
                        whitelistGuideReporterTsv=whitelistGuideReporterTsv_i5_Barcode,
                        umiToolsHeaderBarcodeRegex=input_umiToolsHeaderBarcodeRegex,
                        umiToolsUmiPatternRegex=input_umiToolsUmiPatternRegex,
                        surrogateHammingThresholdStrict=input_surrogateHammingThresholdStrict,
                        barcodeHammingThresholdStrict=input_barcodeHammingThresholdStrict,
                        protospacerHammingThresholdStrict=input_protospacerHammingThresholdStrict
                }
                
                Pair[String, File] GuideCount_i5_Barcode_count_result_pair = (demultiplexedFiles_i5_BarcodeIndex, GuideCount_i5_Barcode.count_result)
            }
            Array[Pair[String, File]] GuideCount_i5_Barcode_count_result_pair_list = GuideCount_i5_Barcode_count_result_pair
            Map[String, File] GuideCount_i5_Barcode_count_result_map = as_map(select_all(GuideCount_i5_Barcode_count_result_pair_list))

            Pair[String, Map[String, File]] GuideCount_i5_Barcode_count_result_map_pair = (output_readIndexMap_i5_Barcode_Map_defined_i5Index, GuideCount_i5_Barcode_count_result_map)
        }
        Array[Pair[String, Map[String, File]]] GuideCount_i5_Barcode_count_result_map_pair_list = GuideCount_i5_Barcode_count_result_map_pair
        Map[String, Map[String, File]] GuideCount_i5_Barcode_count_result_nested_map = as_map(select_all(GuideCount_i5_Barcode_count_result_map_pair_list))

    }

    #
    #   Only demultiplexed by barcode
    #
    if (barcode_demultiplexed){
        Pair[DemultiplexedFiles, UndeterminedFiles] output_DemultiplexedResult_Barcode_defined = select_first([input_DemultiplexedResult_Barcode])
        
        #
        #   Scatter through the i5 indices
        #
        scatter(demultiplexedFiles_barcode in as_pairs(output_DemultiplexedResult_Barcode_defined.left.demultiplexedFiles)){ 
            
            String demultiplexedFiles_barcode_Index = demultiplexedFiles_barcode.left
            IndexPair demultiplexedFiles_barcode_IndexPair = demultiplexedFiles_barcode.right
            
            #
            #   Choose the whitelist_guide_reporter_tsv
            #
            if(defined(input_barcodeOnly_whitelistGuideReporterTsv)){
                Map[String, File] input_barcodeOnly_whitelistGuideReporterTsv_defined = select_first([input_barcodeOnly_whitelistGuideReporterTsv])
                File input_barcodeOnly_whitelistGuideReporterTsv_value = input_barcodeOnly_whitelistGuideReporterTsv_defined[demultiplexedFiles_barcode_Index]
            }
            File whitelistGuideReporterTsv_barcode = select_first([input_barcodeOnly_whitelistGuideReporterTsv_value, input_whitelistGuideReporterTsv])
            
            #
            #   Perform guide mapping
            #
            call GuideCount as GuideCount_Barcode {
                input:
                    countInputRead1=demultiplexedFiles_barcode_IndexPair.read1,
                    countInputRead2=demultiplexedFiles_barcode_IndexPair.read2,
                    whitelistGuideReporterTsv=whitelistGuideReporterTsv_barcode,
                    umiToolsHeaderBarcodeRegex=input_umiToolsHeaderBarcodeRegex,
                    umiToolsUmiPatternRegex=input_umiToolsUmiPatternRegex,
                    surrogateHammingThresholdStrict=input_surrogateHammingThresholdStrict,
                    barcodeHammingThresholdStrict=input_barcodeHammingThresholdStrict,
                    protospacerHammingThresholdStrict=input_protospacerHammingThresholdStrict
            }

            Pair[String, File] GuideCount_Barcode_count_result_pair = (demultiplexedFiles_barcode_Index, GuideCount_Barcode.count_result)
        }

        Array[Pair[String, File]] GuideCount_Barcode_count_result_pair_list = GuideCount_Barcode_count_result_pair
        Map[String, File] GuideCount_Barcode_count_result_map = as_map(select_all(GuideCount_Barcode_count_result_pair_list))

    }

    if (not_demultiplexed){
        #
        #   Choose the whitelist_guide_reporter_tsv
        #
        File whitelistGuideReporterTsv_nonIndexed = select_first([input_whitelistGuideReporterTsv])
        
        #
        #   Perform guide mapping
        #
        call GuideCount as GuideCount_NonIndexed {
            input:
                countInputRead1=input_nonDemultiplexedRead1,
                countInputRead2=input_nonDemultiplexedRead2,
                whitelistGuideReporterTsv=whitelistGuideReporterTsv_nonIndexed,
                umiToolsHeaderBarcodeRegex=input_umiToolsHeaderBarcodeRegex,
                umiToolsUmiPatternRegex=input_umiToolsUmiPatternRegex,
                surrogateHammingThresholdStrict=input_surrogateHammingThresholdStrict,
                barcodeHammingThresholdStrict=input_barcodeHammingThresholdStrict,
                protospacerHammingThresholdStrict=input_protospacerHammingThresholdStrict
        }
    }

    output {
        Map[String, File]? output_GuideCount_i5_count_result_map = GuideCount_i5_count_result_map 
        Map[String, Map[String, File]]? output_GuideCount_i5_Barcode_count_result_nested_map = GuideCount_i5_Barcode_count_result_nested_map
        Map[String, File]? output_GuideCount_Barcode_count_result_map = GuideCount_Barcode_count_result_map
    }

}