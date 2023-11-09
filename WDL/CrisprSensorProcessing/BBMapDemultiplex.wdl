version development

import "BBMapDemultiplexPostprocessWorkflow.wdl" as postprocess

struct DemultiplexResult {
    Map[String, IndexPair] demultiplexedFiles
    File undeterminedR1File
    File? undeterminedR2File
}

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
        
        command {
            demuxbyname.sh in=${processedRead1} \ 
                ${if defined(processedRead2) then 'in2=${processedRead2}' else ' ' } \
                out=${demultiplexedOutputFilenamePattern} \
                outu=${demultiplexedUndeterminedFilenamePattern} \ 
                delimiter='${delimiter}' \
                column=${column} \
                names=${sep(",", indexStringList)} \
                hdist=${hamming}
        }
        

        output {
            Array[File]? demultiplexedOutputR1Files = glob(sub(sub(demultiplexedOutputFilenamePattern,"%", "*"), "#", "1"))
            Array[File]? demultiplexedOutputR2Files = glob(sub(sub(demultiplexedOutputFilenamePattern,"%", "*"),"#", "2"))
            File? undeterminedOutputR1File = sub(demultiplexedUndeterminedFilenamePattern, "#", "1")
            File? undeterminedOutputR2File = sub(demultiplexedUndeterminedFilenamePattern, "#", "1")
        }

        runtime {
                docker: "quay.io/biocontainers/bbmap:39.01--h92535d8_1"
        }
}

# TODO: Deal with the undetermined files. Will probabily need another struct that captures the Map and undetermined files.
workflow BBMapDemultiplexOrchestratorWorkflow {
    
    input {
        File inputRead1
        File? inputRead2

        Array[String]? i5IndexListFinal
        Array[String]? barcodeIndexListFinal
        
        Int i5Hamming
        Int barcodeHamming

        String sampleName
    }

    String rootDemultiplexedI5OutputFilenamePattern = ".${sampleName}.Demultiplexed.R#.fastq"
    String demultiplexedI5OutputFilenamePattern = "%" + rootDemultiplexedI5OutputFilenamePattern
    String demultiplexedI5UndeterminedFilenamePattern = "undetermined" + rootDemultiplexedI5OutputFilenamePattern

    String rootDemultiplexedi5DemultiplexedBarcodeOutputFilenamePattern = ".{i5Index}.${sampleName}.Demultiplexed.R#.fastq"
    String demultiplexedi5DemultiplexedBarcodeOutputFilenamePattern = "%" + rootDemultiplexedi5DemultiplexedBarcodeOutputFilenamePattern
    String demultiplexedI5DemultiplexedBarcodeUndeterminedFilenamePattern = "undetermined" + rootDemultiplexedi5DemultiplexedBarcodeOutputFilenamePattern
    
    if (defined(i5IndexListFinal)) {
        if(defined(inputRead2)){
            # Perform i5 demultiplex on R1&R2
            call BBMapDemultiplexRunnerTask as BBMapDemultiplexRunnerTask_i5 {
                input:
                        processedRead1=inputRead1,
                        processedRead2=inputRead2,
                        indexStringList=select_first([i5IndexListFinal]),
                        delimiter="\\+",
                        column=2,
                        hamming=i5Hamming,
                        demultiplexedOutputFilenamePattern=demultiplexedI5OutputFilenamePattern,
                        demultiplexedUndeterminedFilenamePattern=demultiplexedI5UndeterminedFilenamePattern
            }
            
            if (defined(BBMapDemultiplexRunnerTask_i5.demultiplexedOutputR1Files)){
                call postprocess.BBMapDemultiplexRunnerPostprocessWorkflow as BBMapDemultiplexRunnerPostprocessWorkflow_i5 {
                    input:
                        demultiplexedOutputR1Files=select_first([BBMapDemultiplexRunnerTask_i5.demultiplexedOutputR1Files]),
                        demultiplexedOutputR2Files=BBMapDemultiplexRunnerTask_i5.demultiplexedOutputR2Files,
                        rootDemultiplexedOutputFilenamePattern=rootDemultiplexedI5OutputFilenamePattern
                }
                Map[String, IndexPair] readIndexMap_i5_Map = BBMapDemultiplexRunnerPostprocessWorkflow_i5.output_readIndexMap
                DemultiplexResult demultiplexResult_i5 = {
                    "demultiplexedFiles":readIndexMap_i5_Map,
                    "undeterminedR1File":select_first([BBMapDemultiplexRunnerTask_i5.undeterminedOutputR1File]),
                    "undeterminedR2File":BBMapDemultiplexRunnerTask_i5.undeterminedOutputR2File
                }

                if(defined(barcodeIndexListFinal)){
                    # For each i5 demultiplex, demultiplex by barcode
                    scatter(readIndexPair_i5 in as_pairs(readIndexMap_i5_Map)){
                        call BBMapDemultiplexRunnerTask as BBMapDemultiplexRunnerTask_i5_Barcode {
                            input:
                                processedRead1=readIndexPair_i5.right.read1,
                                processedRead2=readIndexPair_i5.right.read2,
                                indexStringList=select_first([barcodeIndexListFinal]),
                                delimiter="_",
                                column=2,
                                hamming=barcodeHamming,
                                demultiplexedOutputFilenamePattern=sub(demultiplexedi5DemultiplexedBarcodeOutputFilenamePattern, "{i5Index}", readIndexPair_i5.right.index),
                                demultiplexedUndeterminedFilenamePattern=sub(demultiplexedI5DemultiplexedBarcodeUndeterminedFilenamePattern,"{i5Index}", readIndexPair_i5.right.index)
                        }

                        #BBMapDemultiplexRunnerTask_i5_Barcode.undeterminedOutputR1File
                        #BBMapDemultiplexRunnerTask_i5_Barcode.undeterminedOutputR2File
                        if (defined(BBMapDemultiplexRunnerTask_i5_Barcode.demultiplexedOutputR1Files)){
                            call postprocess.BBMapDemultiplexRunnerPostprocessWorkflow as BBMapDemultiplexRunnerPostprocessWorkflow_i5_Barcode {
                                input:
                                    demultiplexedOutputR1Files=select_first([BBMapDemultiplexRunnerTask_i5_Barcode.demultiplexedOutputR1Files]),
                                    demultiplexedOutputR2Files=BBMapDemultiplexRunnerTask_i5_Barcode.demultiplexedOutputR2Files,
                                    rootDemultiplexedOutputFilenamePattern=sub(rootDemultiplexedi5DemultiplexedBarcodeOutputFilenamePattern, "{i5Index}", readIndexPair_i5.right.index)
                            }
                            Map[String, IndexPair] readIndexMap_i5_Barcode_Map = BBMapDemultiplexRunnerPostprocessWorkflow_i5_Barcode.output_readIndexMap
                            DemultiplexResult demultiplexResult_i5_Barcode = {
                                "demultiplexedFiles":readIndexMap_i5_Barcode_Map,
                                "undeterminedR1File":select_first([BBMapDemultiplexRunnerTask_i5_Barcode.undeterminedOutputR1File]),
                                "undeterminedR2File":BBMapDemultiplexRunnerTask_i5_Barcode.undeterminedOutputR2File
                            }
                            Pair[String, Pair[IndexPair, DemultiplexResult]]? readIndexMap_i5_Barcode = (readIndexPair_i5.left, (readIndexPair_i5.right, demultiplexResult_i5_Barcode))


                        }
                    }
                    Array[Pair[String, Pair[IndexPair, DemultiplexResult]]?] readIndexMap_i5_Barcode_List = readIndexMap_i5_Barcode
                    Map[String, Pair[IndexPair, DemultiplexResult]] readIndexMap_i5_Barcode_Map = as_map(select_all(readIndexMap_i5_Barcode_List))
                }
            }
        }
    } 
    if (!defined(i5IndexListFinal)) {
        if(defined(barcodeIndexListFinal)){
            call BBMapDemultiplexRunnerTask as BBMapDemultiplexRunnerTask_Barcode {
                    input:
                        processedRead1=inputRead1,
                        processedRead2=inputRead2,
                        indexStringList=select_first([barcodeIndexListFinal]),
                        delimiter="_",
                        column=2,
                        hamming=barcodeHamming,
                        demultiplexedOutputFilenamePattern=sub(demultiplexedi5DemultiplexedBarcodeOutputFilenamePattern,"{i5Index}", "NA"),
                        demultiplexedUndeterminedFilenamePattern=sub(demultiplexedI5DemultiplexedBarcodeUndeterminedFilenamePattern,"{i5Index}", "NA")
            }
            #BBMapDemultiplexRunnerTask_Barcode.undeterminedOutputR1File
            #BBMapDemultiplexRunnerTask_Barcode.undeterminedOutputR2File
            if (defined(BBMapDemultiplexRunnerTask_Barcode.demultiplexedOutputR1Files)){
                call postprocess.BBMapDemultiplexRunnerPostprocessWorkflow as BBMapDemultiplexRunnerPostprocessWorkflow_Barcode {
                    input:
                        demultiplexedOutputR1Files=select_first([BBMapDemultiplexRunnerTask_Barcode.demultiplexedOutputR1Files]),
                        demultiplexedOutputR2Files=BBMapDemultiplexRunnerTask_Barcode.demultiplexedOutputR2Files,
                        rootDemultiplexedOutputFilenamePattern=sub(rootDemultiplexedi5DemultiplexedBarcodeOutputFilenamePattern,"{i5Index}", "NA")
                }
                Map[String, IndexPair] readIndexMap_Barcode_Map = BBMapDemultiplexRunnerPostprocessWorkflow_Barcode.output_readIndexMap
                DemultiplexResult demultiplexResult_Barcode = {
                    "demultiplexedFiles":readIndexMap_Barcode_Map,
                    "undeterminedR1File":select_first([BBMapDemultiplexRunnerTask_Barcode.undeterminedOutputR1File]),
                    "undeterminedR2File":BBMapDemultiplexRunnerTask_Barcode.undeterminedOutputR2File
                }
            }


        }
    }
    output {
        DemultiplexResult? output_DemultiplexResult_i5 = demultiplexResult_i5
        Map[String, Pair[IndexPair, DemultiplexResult]]? output_readIndexMap_i5_Barcode_Map = readIndexMap_i5_Barcode_Map
        DemultiplexResult? output_DemultiplexResult_Barcode = demultiplexResult_Barcode
    }
    
}