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
        DemultiplexedFiles demultiplexedFiles
        File whitelistGuideReporterTsv
        String? umiToolsHeaderBarcodeRegex
        String? umiToolsUmiPatternRegex
        Int? surrogateHammingThresholdStrict
        Int? barcodeHammingThresholdStrict
        Int? protospacerHammingThresholdStrict
    }

    IndexPair indexPair = demultiplexedFiles.demultiplexedFiles.right
    File read1 = indexPair.read1
    File? read2 = indexPair.read2


    command <<<
        python <<CODE
            import crispr_ambiguous_mapping
            import pandas as pd
            
            whitelist_guide_reporter_df = pd.read_table("~{whitelistGuideReporterTsv}")

            result = crispr_ambiguous_mapping.mp.get_whitelist_reporter_counts_from_umitools_output(
                whitelist_guide_reporter_df=whitelist_guide_reporter_df, 
                fastq_r1_fn='~{read1}', 
                fastq_r2_fn=~{if defined(read2) then "'~{read2}'" else "None" },
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
        String r1U6BarcodeExtractionRegex
        String? r24ntBarcodeExtractionRegex
        String sampleName

        File? i5IndexStringTextFile
        Array[String]? i5IndexList
        File? barcodeIndexStringTextFile
        Array[String]? barcodeIndexList
        
        Int i5Hamming
        Int barcodeHamming
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

    if defined(demultiplexWorkflow.output_DemultiplexedResult_i5){
        Pair[DemultiplexedFiles, UndeterminedFiles] output_DemultiplexedResult_i5_defined = select_first([demultiplexWorkflow.output_DemultiplexedResult_i5])
        call demultiplex.BBMapDemultiplexOrchestratorWorkflow as demultiplexWorkflow {
            input:
                demultiplexedFiles=demultiplexedFiles,
                whitelistGuideReporterTsv=whitelistGuideReporterTsv,
                umiToolsHeaderBarcodeRegex=umiToolsHeaderBarcodeRegex,
                umiToolsUmiPatternRegex=umiToolsUmiPatternRegex,
                surrogateHammingThresholdStrict=surrogateHammingThresholdStrict,
                barcodeHammingThresholdStrict=surrogateHammingThresholdStrict,
                protospacerHammingThresholdStrict=protospacerHammingThresholdStrict
        }
        # LEFTOFF, The inputs above needs to be updated to be workflow updates. Continue on the below defined blocks, and add the outputs. Figure out how to pass guide libraries as a Map.
        output_DemultiplexedResult_i5_defined.left

    }
    
    if defined(demultiplexWorkflow.output_readIndexMap_i5_Barcode_Map){
        Map[String, Pair[IndexPair, Pair[DemultiplexedFiles, UndeterminedFiles]]] output_readIndexMap_i5_Barcode_Map_defined = select_first([demultiplexWorkflow.output_readIndexMap_i5_Barcode_Map])
    }

    if defined(demultiplexWorkflow.output_DemultiplexedResult_Barcode){
        Pair[DemultiplexedFiles, UndeterminedFiles] output_DemultiplexedResult_Barcode_defined = select_first([demultiplexWorkflow.output_DemultiplexedResult_Barcode])
    }


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


