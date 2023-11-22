version development

import "BBMapDemultiplex.wdl" as demultiplex

task UmiToolsExtractTask {
    input {
        File inputRead1
        File? inputRead2
        String r1U6BarcodeExtractionRegex
        String? r24ntBarcodeExtractionRegex
        String sampleName
    }

	String logFileFn="${sampleName}.umitools_log.txt"
	String outputRead1Fn="${sampleName}.umitools.R1.fastq"
	String outputRead2Fn="${sampleName}.umitools.R2.fastq"
	String outputFilteredRead1Fn="${sampleName}.umitools.filtered.R1.fastq"
	String outputFilteredRead2Fn="${sampleName}.umitools.filtered.R2.fastq"
	
    Boolean inputRead2Defined = defined(inputRead2)
    Boolean r24ntBarcodeExtractionRegexDefined = defined(r24ntBarcodeExtractionRegex)

	command <<<
        if [ ~{inputRead2Defined} ] && [ ~{r24ntBarcodeExtractionRegexDefined} ]; then
            echo "Running R1+R2 Mode"
            umi_tools extract --extract-method=regex \
            --stdin ~{inputRead1} \
            --read2-in ~{inputRead2} \
            --bc-pattern ~{r1U6BarcodeExtractionRegex} \
            --bc-pattern2 ~{r24ntBarcodeExtractionRegex} \
            -L ~{logFileFn} \
            --stdout ~{outputRead1Fn} \
            --read2-out ~{outputRead2Fn} \
            --filtered-out ~{outputFilteredRead1Fn} \
            --filtered-out2 ~{outputFilteredRead2Fn}
            echo "Completed"
        else
            echo "Running R1-only Mode"
            umi_tools extract --extract-method=regex \
            --stdin ~{inputRead1} \
            --bc-pattern=~{r1U6BarcodeExtractionRegex} \
            -L ~{logFileFn} \
            --stdout ~{outputRead1Fn} \
            --filtered-out ~{outputFilteredRead1Fn}
            echo "Completed"
        fi

        awk '{print $1/4}' <(wc -l < ~{outputRead1Fn}) > output_read_count.txt
        awk '{print $1/4}' <(wc -l < ~{outputFilteredRead1Fn}) > output_filtered_read_count.txt

        echo "Output read count:"
        cat output_read_count.txt
        echo "Output filtered read count:"
        cat output_filtered_read_count.txt
    >>>

    output {
		File outputRead1="${outputRead1Fn}"
		File? outputRead2="${outputRead2Fn}"
		File outputFilteredRead1="${outputFilteredRead1Fn}"
		File? outputFilteredRead2="${outputFilteredRead2Fn}"
		File logFile="${logFileFn}"
		Float outputRead1ReadCount=read_float('output_read_count.txt')
		Float outputFilteredRead1ReadCount=read_float('output_filtered_read_count.txt')
	}

	runtime {
		docker: "quay.io/biocontainers/umi_tools:1.1.4--py38he5da3d1_2"
	}
}

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

workflow CrisprSensorPreprocessing_Workflow {
    input {
        File rawFastqR1
        File? rawFastqR2
        String r1U6BarcodeExtractionRegex # TODO: Rename to distinguish between the other REGEX
        String? r24ntBarcodeExtractionRegex # TODO: Rename to distinguish between the other REGEX
        String sampleName

        File? i5IndexStringTextFile
        Array[String]? i5IndexList
        File? barcodeIndexStringTextFile
        Array[String]? barcodeIndexList
        
        Int i5Hamming 
        Int barcodeHamming


        File? input_whitelistGuideReporterTsv
        Map[String, File]? input_i5Only_whitelistGuideReporterTsv
        Map[String, File]? input_barcodeOnly_whitelistGuideReporterTsv
        Map[String, Map[String, File]]? input_i5Barrcode_whitelistGuideReporterTsv

        String? input_umiToolsHeaderBarcodeRegex
        String? input_umiToolsUmiPatternRegex
        
        Int? input_surrogateHammingThresholdStrict
        Int? input_barcodeHammingThresholdStrict
        Int? input_protospacerHammingThresholdStrict
    }

    # Extract relevant sequences using UMI tools
    call UmiToolsExtractTask {
        input:
            inputRead1=rawFastqR1,
            inputRead2=rawFastqR2,
            r1U6BarcodeExtractionRegex=r1U6BarcodeExtractionRegex,
            r24ntBarcodeExtractionRegex=r24ntBarcodeExtractionRegex,
            sampleName=sampleName
    }
	
    # Demultiplex the UMI tools result
    call demultiplex.BBMapDemultiplexOrchestratorWorkflow as demultiplexWorkflow {
        input:
            inputRead1=UmiToolsExtractTask.outputRead1,
            inputRead2=UmiToolsExtractTask.outputRead2,
            i5IndexStringTextFile=i5IndexStringTextFile,
            i5IndexList=i5IndexList,
            barcodeIndexStringTextFile=barcodeIndexStringTextFile,
            barcodeIndexList=barcodeIndexList,
            i5Hamming=i5Hamming,
            barcodeHamming=barcodeHamming,
            sampleName=sampleName
    }

    Boolean i5_demultiplexed = defined(demultiplexWorkflow.output_DemultiplexedResult_i5) && !defined(demultiplexWorkflow.output_readIndexMap_i5_Barcode_Map)
    Boolean i5_Barcode_demultiplexed = defined(demultiplexWorkflow.output_readIndexMap_i5_Barcode_Map)
    Boolean barcode_demultiplexed = defined(demultiplexWorkflow.output_DemultiplexedResult_Barcode)
    Boolean not_demultiplexed = !i5_demultiplexed && !i5_Barcode_demultiplexed && !barcode_demultiplexed

    if (i5_demultiplexed){
        Pair[DemultiplexedFiles, UndeterminedFiles] output_DemultiplexedResult_i5_defined = select_first([demultiplexWorkflow.output_DemultiplexedResult_i5])
        
        #
        #   Scatter through the i5 indices
        #
        scatter(demultiplexedFiles_i5 in output_DemultiplexedResult_i5_defined.left.demultiplexedFiles){ 
            
            String demultiplexedFiles_i5_Index = demultiplexedFiles_i5.left
            IndexPair demultiplexedFiles_i5_IndexPair = demultiplexedFiles_i5.right
            
            #
            #   Choose the whitelist_guide_reporter_tsv
            #
            if defined(input_i5Only_whitelistGuideReporterTsv){
                File input_i5Only_whitelistGuideReporterTsv_value = input_i5Only_whitelistGuideReporterTsv[demultiplexedFiles_i5_Index]
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
        }

        # TODO: Create Map from the count result in a similar format to the demultiplex result.
        GuideCount_i5.count_result

    }
    
    if (i5_Barcode_demultiplexed){
        Map[String, Pair[IndexPair, Pair[DemultiplexedFiles, UndeterminedFiles]]] output_readIndexMap_i5_Barcode_Map_defined = select_first([demultiplexWorkflow.output_readIndexMap_i5_Barcode_Map])
        
        #
        #   Scatter through the i5 indices
        #
        scatter(output_readIndexMap_i5_Barcode_Map_defined_i5Pair in output_readIndexMap_i5_Barcode_Map_defined) { 
            String output_readIndexMap_i5_Barcode_Map_defined_i5Index = output_readIndexMap_i5_Barcode_Map_defined_i5Pair.left
            Pair[DemultiplexedFiles, UndeterminedFiles] output_readIndexMap_i5_Barcode_Map_defined_i5_barcodeFiles = output_readIndexMap_i5_Barcode_Map_defined_i5Pair.right.right


            #
            #   Scatter through the barcode indices
            #
            scatter(demultiplexedFiles_i5_Barcode in output_readIndexMap_i5_Barcode_Map_defined_i5_barcodeFiles.left.demultiplexedFiles){ # ITERATE THROUGH i5 INDEXES
            
                String demultiplexedFiles_i5_BarcodeIndex = demultiplexedFiles_i5_Barcode.left
                IndexPair demultiplexedFiles_i5_Barcode_IndexPair = demultiplexedFiles_i5_Barcode.right
                
                #
                #   Choose the whitelist_guide_reporter_tsv
                #
                if defined(input_i5Barcode_whitelistGuideReporterTsv){
                    File input_i5Barcode_whitelistGuideReporterTsv_value = input_i5Barcode_whitelistGuideReporterTsv[output_readIndexMap_i5_Barcode_Map_defined_i5Index][demultiplexedFiles_i5_BarcodeIndex]
                }
                if defined(input_i5Only_whitelistGuideReporterTsv){
                    File input_i5Only_whitelistGuideReporterTsv_value = input_i5Only_whitelistGuideReporterTsv[output_readIndexMap_i5_Barcode_Map_defined_i5Index]
                }
                if defined(input_barcodeOnly_whitelistGuideReporterTsv){
                    File input_barcodeOnly_whitelistGuideReporterTsv_value = input_barcodeOnly_whitelistGuideReporterTsv[demultiplexedFiles_i5_BarcodeIndex]
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
            }

            # TODO: Create Map from the count result in a similar format to the demultiplex result.
        }

        # TODO: Create Map from the count result in a similar format to the demultiplex result.
    }

    if (barcode_demultiplexed){
        Pair[DemultiplexedFiles, UndeterminedFiles] output_DemultiplexedResult_Barcode_defined = select_first([demultiplexWorkflow.output_DemultiplexedResult_Barcode])
        
        #
        #   Scatter through the i5 indices
        #
        scatter(demultiplexedFiles_barcode in output_DemultiplexedResult_Barcode_defined.left.demultiplexedFiles){ 
            
            String demultiplexedFiles_barcode_Index = demultiplexedFiles_barcode.left
            IndexPair demultiplexedFiles_barcode_IndexPair = demultiplexedFiles_barcode.right
            
            #
            #   Choose the whitelist_guide_reporter_tsv
            #
            if defined(input_barcodeOnly_whitelistGuideReporterTsv){
                File input_barcodeOnly_whitelistGuideReporterTsv_value = input_barcodeOnly_whitelistGuideReporterTsv[demultiplexedFiles_barcode_Index]
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
        }

        # TODO: Create Map from the count result in a similar format to the demultiplex result.
        GuideCount_Barcode.count_result

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
                countInputRead1=UmiToolsExtractTask.outputRead1,
                countInputRead2=UmiToolsExtractTask.outputRead2,
                whitelistGuideReporterTsv=whitelistGuideReporterTsv_nonIndexed,
                umiToolsHeaderBarcodeRegex=input_umiToolsHeaderBarcodeRegex,
                umiToolsUmiPatternRegex=input_umiToolsUmiPatternRegex,
                surrogateHammingThresholdStrict=input_surrogateHammingThresholdStrict,
                barcodeHammingThresholdStrict=input_barcodeHammingThresholdStrict,
                protospacerHammingThresholdStrict=input_protospacerHammingThresholdStrict
        }
    }
    # LEFTOFF - Organize the guide count results and set as output, then test. After testing, do the TODO before to use the sample sheet

	# TODO: Since we have a MAP, perhaps someone can input a TSV of the sample sheet, then we can return a final table with the samples attached. Or we can return a table without the sample sheet.

    output {
        File umiToolsExtractedOutputRead1 = UmiToolsExtractTask.outputRead1
        File? umiToolsExtractedOutputRead2 = UmiToolsExtractTask.outputRead2
        File umiToolsExtractedOutputFilteredRead1=UmiToolsExtractTask.outputFilteredRead1
        File? umiToolsExtractedOutputFilteredRead2=UmiToolsExtractTask.outputFilteredRead2
        File umiToolsExtractTaskLogFile=UmiToolsExtractTask.logFile
        Float umiToolsExtractOutputRead1ReadCount = UmiToolsExtractTask.outputRead1ReadCount
        Float umiToolsExtractOutputFilteredRead1ReadCount = UmiToolsExtractTask.outputFilteredRead1ReadCount

        Pair[DemultiplexedFiles, UndeterminedFiles]? output_DemultiplexedResult_i5 = demultiplexWorkflow.output_DemultiplexedResult_i5
        Map[String, Pair[IndexPair, Pair[DemultiplexedFiles, UndeterminedFiles]]]? output_readIndexMap_i5_Barcode_Map = demultiplexWorkflow.output_readIndexMap_i5_Barcode_Map
        Pair[DemultiplexedFiles, UndeterminedFiles]? output_DemultiplexedResult_Barcode = demultiplexWorkflow.output_DemultiplexedResult_Barcode
    }
}


