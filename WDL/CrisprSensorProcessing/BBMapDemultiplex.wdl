version development

import "BBMapDemultiplexPostprocessWorkflow.wdl" as postprocess

struct AnnotatedSample {
    String screenId
    String? i5Index
    String? barcodeIndex
    File read1
    File? read2
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

        # Mapping from indices to screen name
        String? screenId
        Map[String, String]? i5ToScreenidMap
        Map[String, String]? barcodeToScreenidMap
        Map[String, Map[String, String]]? i5ToBarcodeToScreenidMap

        # Mapping from indices to sample annotations
        Array[String]? sampleInfoVars
        Map[String, Array[String]]? i5ToSampleInfoVarsMap
        Map[String, Array[String]]? barcodeToSampleInfoVarsMap
        Map[String, Map[String, Array[String]]]? i5ToBarcodeToSampleInfoVarsMap
    }



    #
    # Process the indices
    #
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
	
    #
    # Prepare the demultiplex filename patterns. The "root" filenames are there to be used with "basename" function to extract the indices from the globbed filenames post-demultiplex
    #
    String rootDemultiplexedI5OutputFilenamePattern = ".${sampleName}.Demultiplexed.R#.fastq"
    String demultiplexedI5OutputFilenamePattern = "%" + rootDemultiplexedI5OutputFilenamePattern
    String demultiplexedI5UndeterminedFilenamePattern = "${sampleName}.Undetermined.R#.fastq"

    String rootDemultiplexedi5DemultiplexedBarcodeOutputFilenamePattern = ".'i5Index'.${sampleName}.Demultiplexed.R#.fastq"
    String demultiplexedi5DemultiplexedBarcodeOutputFilenamePattern = "%" + rootDemultiplexedi5DemultiplexedBarcodeOutputFilenamePattern
    String demultiplexedI5DemultiplexedBarcodeUndeterminedFilenamePattern = "'i5Index'.${sampleName}.Undetermined.R#.fastq"
    
    # Boolean indicators if i5 or barcode is available
    Boolean i5IndexAvailable = defined(i5IndexList) || defined(i5IndexStringTextFileList) 
    Boolean barcodeIndexAvailable = defined(barcodeIndexList) || defined(barcodeIndexStringTextFileList)

    #
    # If i5 indices available, demultiplex the i5
    #
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
        
        #
        # If barcode indices available, demultiplex the barcode
        #
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


            #
            # Convert into Map of screenID to sample.
            #
            scatter(readIndexMap_i5_Barcode_Map_as_pairs in as_pairs(readIndexMap_i5_Barcode_Map)){
                String readIndexMap_i5_Barcode_Map_as_pairs_i5Index = readIndexMap_i5_Barcode_Map_as_pairs.left
                Pair[DemultiplexedFiles, UndeterminedFiles] readIndexMap_i5_Barcode_Map_as_pairs_barcodeFiles = readIndexMap_i5_Barcode_Map_as_pairs.right.right

                scatter(demultiplexedFiles_i5_Barcode in as_pairs(readIndexMap_i5_Barcode_Map_as_pairs_barcodeFiles.left.demultiplexedFiles)){ # ITERATE THROUGH i5 INDEXES
            
                    String demultiplexedFiles_i5_BarcodeIndex = demultiplexedFiles_i5_Barcode.left
                    IndexPair demultiplexedFiles_i5_Barcode_IndexPair = demultiplexedFiles_i5_Barcode.right
                    

                    #
                    #   Get the sample screen ID
                    #
                    if(defined(i5ToBarcodeToScreenidMap)){
                        Map[String, Map[String, String]] i5ToBarcodeToScreenidMap_defined = select_first([i5ToBarcodeToScreenidMap])
                        String i5ToBarcodeToScreenidMap_screenId = i5ToBarcodeToScreenidMap_defined[readIndexMap_i5_Barcode_Map_as_pairs_i5Index][demultiplexedFiles_i5_BarcodeIndex]
                    }
                    if(defined(i5ToScreenidMap)){
                        Map[String, String] i5ToScreenidMap_defined = select_first([i5ToScreenidMap])
                        String i5ToScreenidMap_screenId = i5ToScreenidMap_defined[readIndexMap_i5_Barcode_Map_as_pairs_i5Index]
                    }
                    if(defined(barcodeToScreenidMap)){
                        Map[String, String] barcodeToScreenidMap_defined = select_first([barcodeToScreenidMap])
                        String barcodeToScreenidMap_screenId = barcodeToScreenidMap_defined[demultiplexedFiles_i5_BarcodeIndex]
                    }
                    String i5Barcode_screenId_selected = select_first([i5ToBarcodeToScreenidMap_screenId, i5ToScreenidMap_screenId, barcodeToScreenidMap_screenId, screenId, "nullScreenId"])
                    

                    #
                    #   Get the sample annotation
                    #
                    if(defined(i5ToBarcodeToSampleInfoVarsMap)){
                        Map[String, Map[String, Array[String]]] i5ToBarcodeToSampleInfoVarsMap_defined = select_first([i5ToBarcodeToSampleInfoVarsMap])
                        Array[String] i5ToBarcodeToSampleInfoVarsMap_sampleAnnotation = i5ToBarcodeToSampleInfoVarsMap_defined[readIndexMap_i5_Barcode_Map_as_pairs_i5Index][demultiplexedFiles_i5_BarcodeIndex]
                    }
                    if(defined(i5ToSampleInfoVarsMap)){
                        Map[String, Array[String]] i5ToSampleInfoVarsMap_defined = select_first([i5ToSampleInfoVarsMap])
                        Array[String] i5ToSampleInfoVarsMap_sampleAnnotation = i5ToSampleInfoVarsMap_defined[readIndexMap_i5_Barcode_Map_as_pairs_i5Index]
                    }
                    if(defined(barcodeToSampleInfoVarsMap)){
                        Map[String, Array[String]] barcodeToSampleInfoVarsMap_defined = select_first([barcodeToSampleInfoVarsMap])
                        Array[String] barcodeToSampleInfoVarsMap_sampleAnnotation = barcodeToSampleInfoVarsMap_defined[demultiplexedFiles_i5_BarcodeIndex]
                    }
                    Array[String] i5Barcode_screenAnnotation_selected = select_first([i5ToBarcodeToSampleInfoVarsMap_sampleAnnotation, i5ToSampleInfoVarsMap_sampleAnnotation, barcodeToSampleInfoVarsMap_sampleAnnotation, sampleInfoVars])


                    AnnotatedSample i5Barcode_annotatedSample = {
                        "screenId":i5Barcode_screenId_selected,
                        "i5Index":readIndexMap_i5_Barcode_Map_as_pairs_i5Index,
                        "barcodeIndex": demultiplexedFiles_i5_BarcodeIndex,
                        "read1":demultiplexedFiles_i5_Barcode_IndexPair.read1,
                        "read2":demultiplexedFiles_i5_Barcode_IndexPair.read2
                    }
                    Pair[String, Pair[AnnotatedSample, Array[String]]] i5Barcode_screenIdAnnotatedSamplePair = (i5Barcode_screenId_selected, (i5Barcode_annotatedSample, i5Barcode_screenAnnotation_selected))
                }
            }

            Array[Pair[String, Pair[AnnotatedSample, Array[String]]]] i5Barcode_screenIdToSampleArray = flatten(i5Barcode_screenIdAnnotatedSamplePair)
            Map[String, Array[Pair[AnnotatedSample, Array[String]]]] i5Barcode_screenIdToSampleMap = collect_by_key(i5Barcode_screenIdToSampleArray)
        }



        #
        #   If barcode is note available, still convert into Map of screen ID to sample
        #
        if(!barcodeIndexAvailable){
            #
            # Convert into Map of screenID to sample.
            #
            scatter(demultiplexedFiles_i5_pair in as_pairs(demultiplexedResult_i5.left.demultiplexedFiles)){
                String demultiplexedFiles_i5_Index = demultiplexedFiles_i5_pair.left
                IndexPair demultiplexedFiles_i5_IndexPair = demultiplexedFiles_i5_pair.right

                #
                #   Get the sample screen ID
                #
                if(defined(i5ToScreenidMap)){
                    Map[String, String] i5ToScreenidMap_defined = select_first([i5ToScreenidMap])
                    String i5ToScreenidMap_screenId = i5ToScreenidMap_defined[demultiplexedFiles_i5_Index]
                }
                String i5_screenId_selected = select_first([i5ToScreenidMap_screenId, screenId, "nullScreenId"])
                
                #
                #   Get the sample annotation
                #
                if(defined(i5ToSampleInfoVarsMap)){
                    Map[String, Array[String]] i5ToSampleInfoVarsMap_defined = select_first([i5ToSampleInfoVarsMap])
                    Array[String] i5ToSampleInfoVarsMap_sampleAnnotation = i5ToSampleInfoVarsMap_defined[demultiplexedFiles_i5_Index]
                }
                Array[String] i5_screenAnnotation_selected = select_first([i5ToSampleInfoVarsMap_sampleAnnotation, sampleInfoVars])


                AnnotatedSample i5_annotatedSample = {
                    "screenId":i5_screenId_selected,
                    "i5Index":demultiplexedFiles_i5_Index,
                    "read1":demultiplexedFiles_i5_IndexPair.read1,
                    "read2":demultiplexedFiles_i5_IndexPair.read2
                }
                Pair[String, Pair[AnnotatedSample, Array[String]]] i5_screenIdAnnotatedSamplePair = (i5_screenId_selected, (i5_annotatedSample,i5_screenAnnotation_selected))
            }

            Array[Pair[String, Pair[AnnotatedSample, Array[String]]]] i5_screenIdToSampleArray = i5_screenIdAnnotatedSamplePair
            Map[String, Array[Pair[AnnotatedSample, Array[String]]]] i5_screenIdToSampleMap = collect_by_key(i5_screenIdToSampleArray)
        }
    } 

    #
    # If there is no i5 indices available
    #
    if (!i5IndexAvailable) {
        #
        # But there is barcode indices available, then only demultiplex by the barcode
        #
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


            #
            # Convert into Map of screenID to sample.
            #
            scatter(demultiplexedFiles_barcode in as_pairs(demultiplexedResult_Barcode.left.demultiplexedFiles)){ 
                String demultiplexedFiles_barcode_Index = demultiplexedFiles_barcode.left
                IndexPair demultiplexedFiles_barcode_IndexPair = demultiplexedFiles_barcode.right

                #
                #   Get the sample screen ID
                #
                if(defined(barcodeToScreenidMap)){
                    Map[String, String] barcodeToScreenidMap_defined = select_first([barcodeToScreenidMap])
                    String barcodeToScreenidMap_screenId = barcodeToScreenidMap_defined[demultiplexedFiles_barcode_Index]
                }
                String barcode_screenId_selected = select_first([barcodeToScreenidMap_screenId, screenId, "nullScreenId"])
                
                #
                #   Get the sample annotation
                #
                if(defined(barcodeToSampleInfoVarsMap)){
                    Map[String, Array[String]] barcodeToSampleInfoVarsMap_defined = select_first([barcodeToSampleInfoVarsMap])
                    Array[String] barcodeToSampleInfoVarsMap_sampleAnnotation = barcodeToSampleInfoVarsMap_defined[demultiplexedFiles_barcode_Index]
                }
                Array[String] barcode_screenAnnotation_selected = select_first([barcodeToSampleInfoVarsMap_sampleAnnotation, sampleInfoVars])


                AnnotatedSample barcode_annotatedSample = {
                    "screenId":barcode_screenId_selected,
                    "barcodeIndex": demultiplexedFiles_barcode_Index,
                    "read1":demultiplexedFiles_barcode_IndexPair.read1,
                    "read2":demultiplexedFiles_barcode_IndexPair.read2
                }
                Pair[String, Pair[AnnotatedSample, Array[String]]] barcode_screenIdAnnotatedSamplePair = (barcode_screenId_selected, (barcode_annotatedSample,barcode_screenAnnotation_selected))
            }

            Array[Pair[String, Pair[AnnotatedSample, Array[String]]]] barcode_screenIdToSampleArray = barcode_screenIdAnnotatedSamplePair
            Map[String, Array[Pair[AnnotatedSample, Array[String]]]] barcode_screenIdToSampleMap = collect_by_key(barcode_screenIdToSampleArray)
        }
    }

    #
    # Prepare sample annotation for nonDemultiplexed sample (requires that the single sample info vars are provided)
    #
    if (defined(sampleInfoVars)){
        Array[String] sampleInfoVars_defined = select_first([sampleInfoVars])
        
        String nonDemultiplexed_screenId = select_first([screenId, "nullScreenId"])
            AnnotatedSample nonDemultiplexed_annotatedSample = {
                    "screenId":nonDemultiplexed_screenId,
                    "read1":inputRead1,
                    "read2":inputRead2
                }

        Map[String, Array[Pair[AnnotatedSample, Array[String]]]] nonDemultiplexed_screenIdToSampleMap = {nonDemultiplexed_screenId:[(nonDemultiplexed_annotatedSample,sampleInfoVars_defined)]}
    }



    # Select the final screenIdToSampleMap
    # Map[String, Array[Pair[AnnotatedSample, Array[String]]]] final_screenIdToSampleMap = select_first([i5Barcode_screenIdToSampleMap, i5_screenIdToSampleMap, barcode_screenIdToSampleMap, nonDemultiplexed_screenIdToSampleMap])
    Map[String, Array[Pair[AnnotatedSample, Array[String]]]] final_screenIdToSampleMap = select_first([nonDemultiplexed_screenIdToSampleMap])
    output {
        Map[String, Array[Pair[AnnotatedSample, Array[String]]]] output_screenIdToSampleMap = final_screenIdToSampleMap
    }
}