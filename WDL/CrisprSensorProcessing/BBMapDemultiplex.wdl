version development

import "BBMapDemultiplexPostprocessWorkflow.wdl" as postprocess


task BBMapDemultiplexRunnerTask {
        input {
            File processedRead1
            File? processedRead2
            Array[String] indexStringList
            String delimiter
            Int column
            Int hamming
            String demultiplexedOutputFilenamePattern
            String demultiplexedUndeterminedFilenamePattern
        }
        
        command <<<
            demuxbyname.sh in=~{processedRead1} \
            out=~{if defined(processedRead2) then "~{demultiplexedOutputFilenamePattern}" else "~{sub(demultiplexedOutputFilenamePattern, '#', '1')}" } \
                outu=~{if defined(processedRead2) then "~{demultiplexedUndeterminedFilenamePattern}" else "~{sub(demultiplexedUndeterminedFilenamePattern, '#', '1')}" } \
                delimiter='~{delimiter}' \
                column=~{column} \
                names=~{sep(",", indexStringList)} \
                hdist=~{hamming}
        >>>
        

        output {
            Array[File]? demultiplexedOutputR1Files = glob(sub(sub(demultiplexedOutputFilenamePattern,"%", "*"), "#", "1"))
            Array[File]? demultiplexedOutputR2Files = glob(sub(sub(demultiplexedOutputFilenamePattern,"%", "*"),"#", "2"))
            File? undeterminedOutputR1File = "${sub(demultiplexedUndeterminedFilenamePattern, "#", "1")}"
            File? undeterminedOutputR2File = "${sub(demultiplexedUndeterminedFilenamePattern, "#", "2")}"
        }

        runtime {
                docker: "quay.io/biocontainers/bbmap:39.01--h92535d8_1"
        }
}

workflow BBMapDemultiplexOrchestratorWorkflow {
    
    input {
        File inputRead1
        File? inputRead2

        File? i5IndexStringTextFile
        Array[String]? i5IndexList

        File? barcodeIndexStringTextFile
        Array[String]? barcodeIndexList
        
        Int i5Hamming
        Int barcodeHamming

        String sampleName
    }

    # Process the indices
    if(!defined(i5IndexList)){
        if(defined(i5IndexStringTextFile)){
            Array[String] i5IndexStringTextFileList = read_lines(select_first([i5IndexStringTextFile]))
        }
    }
    if(!defined(barcodeIndexList)){
        if(defined(barcodeIndexStringTextFile)){
            Array[String] barcodeIndexStringTextFileList = read_lines(select_first([barcodeIndexStringTextFile]))
        }
    }
	
    # Prepare the demultiplex filename patterns. The "root" filenames are there to be used with "basename" function to extract the indices from the globbed filenames post-demultiplex
    String rootDemultiplexedI5OutputFilenamePattern = ".${sampleName}.Demultiplexed.R#.fastq"
    String demultiplexedI5OutputFilenamePattern = "%" + rootDemultiplexedI5OutputFilenamePattern
    String demultiplexedI5UndeterminedFilenamePattern = "${sampleName}.Undetermined.R#.fastq"

    String rootDemultiplexedi5DemultiplexedBarcodeOutputFilenamePattern = ".'i5Index'.${sampleName}.Demultiplexed.R#.fastq"
    String demultiplexedi5DemultiplexedBarcodeOutputFilenamePattern = "%" + rootDemultiplexedi5DemultiplexedBarcodeOutputFilenamePattern
    String demultiplexedI5DemultiplexedBarcodeUndeterminedFilenamePattern = "'i5Index'.${sampleName}.Undetermined.R#.fastq"
    
    # Boolean indicators if i5 or barcode is available
    Boolean i5IndexAvailable = defined(i5IndexList) || defined(i5IndexStringTextFileList) 
    Boolean barcodeIndexAvailable = defined(barcodeIndexList) || defined(barcodeIndexStringTextFileList)

    # If i5 indices available, demultiplex the i5
    if (i5IndexAvailable) {
        # Perform i5 demultiplex
        call BBMapDemultiplexRunnerTask as BBMapDemultiplexRunnerTask_i5 {
            input:
                    processedRead1=inputRead1,
                    processedRead2=inputRead2,
                    indexStringList=select_first([i5IndexList, i5IndexStringTextFileList]),
                    delimiter="\\\\+",
                    column=2,
                    hamming=i5Hamming,
                    demultiplexedOutputFilenamePattern=demultiplexedI5OutputFilenamePattern,
                    demultiplexedUndeterminedFilenamePattern=demultiplexedI5UndeterminedFilenamePattern
        }
        
        # Organize the i5 results into the DemultiplexResult struct
        call postprocess.BBMapDemultiplexRunnerPostprocessWorkflow as BBMapDemultiplexRunnerPostprocessWorkflow_i5 {
            input:
                demultiplexedOutputR1Files=select_first([BBMapDemultiplexRunnerTask_i5.demultiplexedOutputR1Files]),
                demultiplexedOutputR2Files=BBMapDemultiplexRunnerTask_i5.demultiplexedOutputR2Files,
                undeterminedR1File= select_first([BBMapDemultiplexRunnerTask_i5.undeterminedOutputR1File]),
                undeterminedR2File= BBMapDemultiplexRunnerTask_i5.undeterminedOutputR2File,
                rootDemultiplexedOutputFilenamePattern=rootDemultiplexedI5OutputFilenamePattern
        }
        
        
        Pair[DemultiplexedFiles, UndeterminedFiles] demultiplexedResult_i5 = BBMapDemultiplexRunnerPostprocessWorkflow_i5.output_demultiplexedResults

        # If barcode indices available, demultiplex the barcode
        if(barcodeIndexAvailable){
            
            # For each demultiplexed i5, demultiplex by barcode
            scatter(readIndexPair_i5 in as_pairs(demultiplexedResult_i5.left.demultiplexedFiles)){

                # Perform barcode demultiplex
                call BBMapDemultiplexRunnerTask as BBMapDemultiplexRunnerTask_i5_Barcode {
                    input:
                        processedRead1=readIndexPair_i5.right.read1,
                        processedRead2=readIndexPair_i5.right.read2,
                        indexStringList=select_first([barcodeIndexList, barcodeIndexStringTextFileList]),
                        delimiter="_",
                        column=2,
                        hamming=barcodeHamming,
                        demultiplexedOutputFilenamePattern=sub(demultiplexedi5DemultiplexedBarcodeOutputFilenamePattern, "'i5Index'", readIndexPair_i5.right.index),
                        demultiplexedUndeterminedFilenamePattern=sub(demultiplexedI5DemultiplexedBarcodeUndeterminedFilenamePattern,"'i5Index'", readIndexPair_i5.right.index)
                }

                # Organize the barcode results into the DemultiplexResult struct
                call postprocess.BBMapDemultiplexRunnerPostprocessWorkflow as BBMapDemultiplexRunnerPostprocessWorkflow_i5_Barcode {
                    input:
                        demultiplexedOutputR1Files=select_first([BBMapDemultiplexRunnerTask_i5_Barcode.demultiplexedOutputR1Files]),
                        demultiplexedOutputR2Files=BBMapDemultiplexRunnerTask_i5_Barcode.demultiplexedOutputR2Files,
                        undeterminedR1File= select_first([BBMapDemultiplexRunnerTask_i5_Barcode.undeterminedOutputR1File]),
                        undeterminedR2File= BBMapDemultiplexRunnerTask_i5_Barcode.undeterminedOutputR2File,
                        rootDemultiplexedOutputFilenamePattern=sub(rootDemultiplexedi5DemultiplexedBarcodeOutputFilenamePattern, "'i5Index'", readIndexPair_i5.right.index)
                }
                
                Pair[DemultiplexedFiles, UndeterminedFiles] demultiplexedResult_i5_Barcode = BBMapDemultiplexRunnerPostprocessWorkflow_i5_Barcode.output_demultiplexedResults

                # Prepare the barcode result as a pair for proper coercion into a map outside of scatter
                Pair[String, Pair[IndexPair, Pair[DemultiplexedFiles, UndeterminedFiles]]]? readIndexMap_i5_Barcode = (readIndexPair_i5.left, (readIndexPair_i5.right, demultiplexedResult_i5_Barcode))
            }
            # Coerce the array of barcode results into a list of pairs, then coerce into a map
            Array[Pair[String, Pair[IndexPair, Pair[DemultiplexedFiles, UndeterminedFiles]]]?] readIndexMap_i5_Barcode_List = readIndexMap_i5_Barcode
            Map[String, Pair[IndexPair, Pair[DemultiplexedFiles, UndeterminedFiles]]] readIndexMap_i5_Barcode_Map = as_map(select_all(readIndexMap_i5_Barcode_List))
        }
    } 

    # If there is no i5 indices available
    if (!i5IndexAvailable) {
        # But there is barcode indices available, then only demultiplex by the barcode
        if(barcodeIndexAvailable){

            # Perform barcode demultiplex
            call BBMapDemultiplexRunnerTask as BBMapDemultiplexRunnerTask_Barcode {
                    input:
                        processedRead1=inputRead1,
                        processedRead2=inputRead2,
                        indexStringList=select_first([barcodeIndexList, barcodeIndexStringTextFileList]),
                        delimiter="_",
                        column=2,
                        hamming=barcodeHamming,
                        demultiplexedOutputFilenamePattern=sub(demultiplexedi5DemultiplexedBarcodeOutputFilenamePattern,"'i5Index'", "NA"),
                        demultiplexedUndeterminedFilenamePattern=sub(demultiplexedI5DemultiplexedBarcodeUndeterminedFilenamePattern,"'i5Index'", "NA")
            }

            # Organize the barcode results into the DemultiplexResult struct
            call postprocess.BBMapDemultiplexRunnerPostprocessWorkflow as BBMapDemultiplexRunnerPostprocessWorkflow_Barcode {
                input:
                    demultiplexedOutputR1Files=select_first([BBMapDemultiplexRunnerTask_Barcode.demultiplexedOutputR1Files]),
                    demultiplexedOutputR2Files=BBMapDemultiplexRunnerTask_Barcode.demultiplexedOutputR2Files,
                    undeterminedR1File=select_first([BBMapDemultiplexRunnerTask_Barcode.undeterminedOutputR1File]),
                    undeterminedR2File=BBMapDemultiplexRunnerTask_Barcode.undeterminedOutputR2File,
                    rootDemultiplexedOutputFilenamePattern=sub(rootDemultiplexedi5DemultiplexedBarcodeOutputFilenamePattern,"'i5Index'", "NA")
            }
            
             Pair[DemultiplexedFiles, UndeterminedFiles] demultiplexedResult_Barcode = BBMapDemultiplexRunnerPostprocessWorkflow_Barcode.output_demultiplexedResults
        }
    }

    output {
        Pair[DemultiplexedFiles, UndeterminedFiles]? output_DemultiplexedResult_i5 = demultiplexedResult_i5
        Map[String, Pair[IndexPair, Pair[DemultiplexedFiles, UndeterminedFiles]]]? output_readIndexMap_i5_Barcode_Map = readIndexMap_i5_Barcode_Map
        Pair[DemultiplexedFiles, UndeterminedFiles]? output_DemultiplexedResult_Barcode = demultiplexedResult_Barcode
    }
}