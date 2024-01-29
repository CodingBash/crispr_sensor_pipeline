version development

import "https://api.firecloud.org/ga4gh/v1/tools/pinellolab:BBMapDemultiplexRunnerPostprocessWorkflow/versions/1/plain-WDL/descriptor" as postprocess
import "https://api.firecloud.org/ga4gh/v1/tools/pinellolab:BBMapDemultiplexRunnerTask/versions/3/plain-WDL/descriptor" as demultiplexRunner

struct AnnotatedSample {
    String screenId
    String? index1
    String? index2
    File read1
    File? read2
}

workflow BBMapDemultiplexOrchestratorWorkflow {
    
    input {
        # SPECIFIC READS TO DEMULTIPLEX
        File inputRead1
        File? inputRead2

        # SPECIFY INDICES 
        File? index1StringTextFile
        Array[String]? index1List

        File? index2StringTextFile
        Array[String]? index2List
        
        # SPECIFY HAMMING ACCEPTANCE
        Int index1Hamming
        Int index2Hamming

        # SPECIFY SAMPLE NAME
        String sampleName

        # SPECIFY HOW TO PARSE INDEX1 FROM HEADER
        String index1BBMapHeaderDelimiter = "\\\\+"
        Int index1BBMapHeaderColumn = 2

        # SPECIFY HOW TO PARSE INDEX2 FROM HEADER
        String index2BBMapHeaderDelimiter = "_"
        Int index2BBMapHeaderColumn = 2
        
        #
        # FOR POSTPROCESS ORGANIZATION
        #
        # Mapping from indices to screen name
        String? screenId
        Map[String, String]? index1ToScreenidMap
        Map[String, String]? index2ToScreenidMap
        Map[String, Map[String, String]]? index1ToIndex2ToScreenidMap

        # Mapping from indices to sample annotations
        Array[String]? sampleInfoVars
        Map[String, Array[String]]? index1ToSampleInfoVarsMap
        Map[String, Array[String]]? index2ToSampleInfoVarsMap
        Map[String, Map[String, Array[String]]]? index1ToIndex2ToSampleInfoVarsMap

        # RUNTIME PARAMS
        Int preemptible = 1
        Int diskGB = 10
        Int memoryGB = 2
        Int maxRetries = 0
        String diskType = "SSD"
    }


    #
    # Process the indices
    #
    if(!defined(index1List)){
        if(defined(index1StringTextFile)){
            Array[String] index1StringTextFileList = read_lines(select_first([index1StringTextFile]))
        }
    }
    if(!defined(index2List)){
        if(defined(index2StringTextFile)){
            Array[String] index2StringTextFileList = read_lines(select_first([index2StringTextFile]))
        }
    }
	
    #
    # Prepare the demultiplex filename patterns. The "root" filenames are there to be used with "basename" function to extract the indices from the globbed filenames post-demultiplex
    #
    String rootDemultiplexedIndex1OutputFilenamePattern = ".${sampleName}.Demultiplexed.R#.fastq"
    String demultiplexedIndex1OutputFilenamePattern = "%" + rootDemultiplexedIndex1OutputFilenamePattern
    String demultiplexedIndex1UndeterminedFilenamePattern = "${sampleName}.Undetermined.R#.fastq"

    String rootDemultiplexedIndex1DemultiplexedIndex2OutputFilenamePattern = ".'index1'.${sampleName}.Demultiplexed.R#.fastq"
    String demultiplexedIndex1DemultiplexedIndex2OutputFilenamePattern = "%" + rootDemultiplexedIndex1DemultiplexedIndex2OutputFilenamePattern
    String demultiplexedIndex1DemultiplexedIndex2UndeterminedFilenamePattern = "'index1'.${sampleName}.Undetermined.R#.fastq"
    
    # Boolean indicators if index1 or index2 is available
    Boolean index1Available = defined(index1List) || defined(index1StringTextFileList) 
    Boolean index2Available = defined(index2List) || defined(index2StringTextFileList)

    #
    # If index 1 indices available, demultiplex the index1
    #
    if (index1Available) {
        # Perform index1 demultiplex
        call demultiplexRunner.BBMapDemultiplexRunnerTask as BBMapDemultiplexRunnerTask_index1 {
            input:
                    processedRead1=inputRead1,
                    processedRead2=inputRead2,
                    indexStringList=select_first([index1List, index1StringTextFileList]),
                    delimiter=index1BBMapHeaderDelimiter,
                    column=index1BBMapHeaderColumn,
                    hamming=index1Hamming,
                    demultiplexedOutputFilenamePattern=demultiplexedIndex1OutputFilenamePattern,
                    demultiplexedUndeterminedFilenamePattern=demultiplexedIndex1UndeterminedFilenamePattern,
                    preemptible=preemptible,
                    diskGB=diskGB,
                    memoryGB=memoryGB,
                    maxRetries=maxRetries,
                    diskType=diskType
        }
        
        # Organize the index1 results into the DemultiplexResult struct
        call postprocess.BBMapDemultiplexRunnerPostprocessWorkflow as BBMapDemultiplexRunnerPostprocessWorkflow_index1 {
            input:
                demultiplexedOutputR1Files=select_first([BBMapDemultiplexRunnerTask_index1.demultiplexedOutputR1Files]),
                demultiplexedOutputR2Files=BBMapDemultiplexRunnerTask_index1.demultiplexedOutputR2Files,
                undeterminedR1File= select_first([BBMapDemultiplexRunnerTask_index1.undeterminedOutputR1File]),
                undeterminedR2File= BBMapDemultiplexRunnerTask_index1.undeterminedOutputR2File,
                rootDemultiplexedOutputFilenamePattern=rootDemultiplexedIndex1OutputFilenamePattern
        }
        
        
        Pair[DemultiplexedFiles, UndeterminedFiles] demultiplexedResult_index1 = BBMapDemultiplexRunnerPostprocessWorkflow_index1.output_demultiplexedResults
        
        #
        # If index2 indices available, demultiplex the index2
        #
        if(index2Available){
            
            # For each demultiplexed index1, demultiplex by index2
            scatter(readIndexPair_index1 in as_pairs(demultiplexedResult_index1.left.demultiplexedFiles)){

                # Perform index2 demultiplex
                call demultiplexRunner.BBMapDemultiplexRunnerTask as BBMapDemultiplexRunnerTask_index1_index2 {
                    input:
                        processedRead1=readIndexPair_index1.right.read1,
                        processedRead2=readIndexPair_index1.right.read2,
                        indexStringList=select_first([index2List, index2StringTextFileList]),
                        delimiter=index2BBMapHeaderDelimiter,
                        column=index2BBMapHeaderColumn,
                        hamming=index2Hamming,
                        demultiplexedOutputFilenamePattern=sub(demultiplexedIndex1DemultiplexedIndex2OutputFilenamePattern, "'index1'", readIndexPair_index1.right.index),
                        demultiplexedUndeterminedFilenamePattern=sub(demultiplexedIndex1DemultiplexedIndex2UndeterminedFilenamePattern,"'index1'", readIndexPair_index1.right.index),
                        preemptible=preemptible,
                        diskGB=diskGB,
                        memoryGB=memoryGB,
                        maxRetries=maxRetries,
                        diskType=diskType
                }

                # Organize the index2 results into the DemultiplexResult struct
                call postprocess.BBMapDemultiplexRunnerPostprocessWorkflow as BBMapDemultiplexRunnerPostprocessWorkflow_index1_index2 {
                    input:
                        demultiplexedOutputR1Files=select_first([BBMapDemultiplexRunnerTask_index1_index2.demultiplexedOutputR1Files]),
                        demultiplexedOutputR2Files=BBMapDemultiplexRunnerTask_index1_index2.demultiplexedOutputR2Files,
                        undeterminedR1File= select_first([BBMapDemultiplexRunnerTask_index1_index2.undeterminedOutputR1File]),
                        undeterminedR2File= BBMapDemultiplexRunnerTask_index1_index2.undeterminedOutputR2File,
                        rootDemultiplexedOutputFilenamePattern=sub(rootDemultiplexedIndex1DemultiplexedIndex2OutputFilenamePattern, "'index1'", readIndexPair_index1.right.index)
                }
                
                Pair[DemultiplexedFiles, UndeterminedFiles] demultiplexedResult_index1_index2 = BBMapDemultiplexRunnerPostprocessWorkflow_index1_index2.output_demultiplexedResults

                # Prepare the index2 result as a pair for proper coercion into a map outside of scatter
                Pair[String, Pair[IndexPair, Pair[DemultiplexedFiles, UndeterminedFiles]]]? readIndexMap_index1_index2 = (readIndexPair_index1.left, (readIndexPair_index1.right, demultiplexedResult_index1_index2))
            }
            # Coerce the array of index2 results into a list of pairs, then coerce into a map
            Array[Pair[String, Pair[IndexPair, Pair[DemultiplexedFiles, UndeterminedFiles]]]?] readIndexMap_index1_index2_List = readIndexMap_index1_index2
            Map[String, Pair[IndexPair, Pair[DemultiplexedFiles, UndeterminedFiles]]] readIndexMap_index1_index2_Map = as_map(select_all(readIndexMap_index1_index2_List))

            # VERY IMPORTANT NOTE (maybe causing bug in later larger runs): In test run, these two scatters use ~100G of JVM heap, unsure why. Got some insight here (since the WriteMetadataActor is the culprit) https://github.com/bcbio/bcbio-nextgen/issues/2697. An idea for debugging is to use a hard disk SQL (instead of in-memory) and observe what is being outputted. Another simpler solution is to change localization https://cromwell.readthedocs.io/en/stable/Configuring/#local-filesystem-options
            #
            # Convert into Map of screenID to sample.
            #
            scatter(readIndexMap_index1_index2_Map_as_pairs in as_pairs(readIndexMap_index1_index2_Map)){
                String readIndexMap_index1_index2_Map_as_pairs_INDEX1_STRING = readIndexMap_index1_index2_Map_as_pairs.left
                Pair[DemultiplexedFiles, UndeterminedFiles] readIndexMap_index1_index2_Map_as_pairs_index2Files = readIndexMap_index1_index2_Map_as_pairs.right.right

                scatter(demultiplexedFiles_index1_index2 in as_pairs(readIndexMap_index1_index2_Map_as_pairs_index2Files.left.demultiplexedFiles)){ # ITERATE THROUGH index1 INDEXES
            
                    String demultiplexedFiles_index1_index2_INDEX2_STRING = demultiplexedFiles_index1_index2.left
                    IndexPair demultiplexedFiles_index1_index2_IndexPair = demultiplexedFiles_index1_index2.right
                    
                    #
                    #   Get the sample screen ID
                    #
                    Map[String, Map[String, String]] index1ToIndex2ToScreenidMap_defined = select_first([index1ToIndex2ToScreenidMap])
                    String index1ToIndex2ToScreenidMap_screenId = index1ToIndex2ToScreenidMap_defined[readIndexMap_index1_index2_Map_as_pairs_INDEX1_STRING][demultiplexedFiles_index1_index2_INDEX2_STRING]
                    

                    #
                    #   Get the sample annotation
                    #
                    Map[String, Map[String, Array[String]]] index1ToIndex2ToSampleInfoVarsMap_defined = select_first([index1ToIndex2ToSampleInfoVarsMap])
                    Array[String] index1ToIndex2ToSampleInfoVarsMap_sampleAnnotation = [""]


                    AnnotatedSample index1_index2_annotatedSample = {
                        "screenId":index1ToIndex2ToScreenidMap_screenId,
                        "index1":readIndexMap_index1_index2_Map_as_pairs_INDEX1_STRING,
                        "index2": demultiplexedFiles_index1_index2_INDEX2_STRING,
                        "read1":demultiplexedFiles_index1_index2_IndexPair.read1,
                        "read2":demultiplexedFiles_index1_index2_IndexPair.read2
                    }
                    Pair[String, Pair[AnnotatedSample, Array[String]]] index1_index2_screenIdAnnotatedSamplePair = (index1ToIndex2ToScreenidMap_screenId, (index1_index2_annotatedSample, index1ToIndex2ToSampleInfoVarsMap_sampleAnnotation))
                }
            }

            Array[Pair[String, Pair[AnnotatedSample, Array[String]]]] index1_index2_screenIdToSampleArray = flatten(index1_index2_screenIdAnnotatedSamplePair)
            Map[String, Array[Pair[AnnotatedSample, Array[String]]]] index1_index2_screenIdToSampleMap = collect_by_key(index1_index2_screenIdToSampleArray)
        }



        #
        #   If index 2 is not available, still convert into Map of screen ID to sample
        #
        if(!index2Available){
            #
            # Convert into Map of screenID to sample.
            #
            scatter(demultiplexedFiles_index1_pair in as_pairs(demultiplexedResult_index1.left.demultiplexedFiles)){
                String demultiplexedFiles_INDEX1_STRING = demultiplexedFiles_index1_pair.left
                IndexPair demultiplexedFiles_index1_IndexPair = demultiplexedFiles_index1_pair.right

                #
                #   Get the sample screen ID
                #
                Map[String, String] index1ToScreenidMap_defined = select_first([index1ToScreenidMap])
                String index1ToScreenidMap_screenId = index1ToScreenidMap_defined[demultiplexedFiles_INDEX1_STRING]
                
                #
                #   Get the sample annotation
                #
                Map[String, Array[String]] index1ToSampleInfoVarsMap_defined = select_first([index1ToSampleInfoVarsMap])
                Array[String] index1ToSampleInfoVarsMap_sampleAnnotation = index1ToSampleInfoVarsMap_defined[demultiplexedFiles_INDEX1_STRING]

                AnnotatedSample index1_annotatedSample = {
                    "screenId":index1ToScreenidMap_screenId,
                    "index1":demultiplexedFiles_INDEX1_STRING,
                    "read1":demultiplexedFiles_index1_IndexPair.read1,
                    "read2":demultiplexedFiles_index1_IndexPair.read2
                }
                Pair[String, Pair[AnnotatedSample, Array[String]]] index1_screenIdAnnotatedSamplePair = (index1ToScreenidMap_screenId, (index1_annotatedSample,index1ToSampleInfoVarsMap_sampleAnnotation))
            }

            Array[Pair[String, Pair[AnnotatedSample, Array[String]]]] index1_screenIdToSampleArray = index1_screenIdAnnotatedSamplePair
            Map[String, Array[Pair[AnnotatedSample, Array[String]]]] index1_screenIdToSampleMap = collect_by_key(index1_screenIdToSampleArray)
        }
    } 

    #
    # If there is no index1 indices available
    #
    if (!index1Available) {
        #
        # But there is index2 indices available, then only demultiplex by the index2
        #
        if(index2Available){

            # Perform index2 demultiplex
            call demultiplexRunner.BBMapDemultiplexRunnerTask as BBMapDemultiplexRunnerTask_index2 {
                    input:
                        processedRead1=inputRead1,
                        processedRead2=inputRead2,
                        indexStringList=select_first([index2List, index2StringTextFileList]),
                        delimiter=index2BBMapHeaderDelimiter,
                        column=index2BBMapHeaderColumn,
                        hamming=index2Hamming,
                        demultiplexedOutputFilenamePattern=sub(demultiplexedIndex1DemultiplexedIndex2OutputFilenamePattern,"'index1'", "NA"),
                        demultiplexedUndeterminedFilenamePattern=sub(demultiplexedIndex1DemultiplexedIndex2UndeterminedFilenamePattern,"'index1'", "NA"),
                        preemptible=preemptible,
                        diskGB=diskGB,
                        memoryGB=memoryGB,
                        maxRetries=maxRetries,
                        diskType=diskType
            }

            # Organize the index2 results into the DemultiplexResult struct
            call postprocess.BBMapDemultiplexRunnerPostprocessWorkflow as BBMapDemultiplexRunnerPostprocessWorkflow_index2 {
                input:
                    demultiplexedOutputR1Files=select_first([BBMapDemultiplexRunnerTask_index2.demultiplexedOutputR1Files]),
                    demultiplexedOutputR2Files=BBMapDemultiplexRunnerTask_index2.demultiplexedOutputR2Files,
                    undeterminedR1File=select_first([BBMapDemultiplexRunnerTask_index2.undeterminedOutputR1File]),
                    undeterminedR2File=BBMapDemultiplexRunnerTask_index2.undeterminedOutputR2File,
                    rootDemultiplexedOutputFilenamePattern=sub(rootDemultiplexedIndex1DemultiplexedIndex2OutputFilenamePattern,"'index1'", "NA")
            }
            
             Pair[DemultiplexedFiles, UndeterminedFiles] demultiplexedResult_index2 = BBMapDemultiplexRunnerPostprocessWorkflow_index2.output_demultiplexedResults

            #
            #Convert into Map of screenID to sample.
           # 
            scatter(demultiplexedFiles_index2 in as_pairs(demultiplexedResult_index2.left.demultiplexedFiles)){ 
                String demultiplexedFiles_INDEX2_STRING = demultiplexedFiles_index2.left
                IndexPair demultiplexedFiles_index2_IndexPair = demultiplexedFiles_index2.right

                #
                #   Get the sample screen ID
                #
                Map[String, String] index2ToScreenidMap_defined = select_first([index2ToScreenidMap])
                String index2ToScreenidMap_screenId = index2ToScreenidMap_defined[demultiplexedFiles_INDEX2_STRING]
                
                #
                #   Get the sample annotation
                #
                Map[String, Array[String]] index2ToSampleInfoVarsMap_defined = select_first([index2ToSampleInfoVarsMap])
                Array[String] index2ToSampleInfoVarsMap_sampleAnnotation = index2ToSampleInfoVarsMap_defined[demultiplexedFiles_INDEX2_STRING]


                AnnotatedSample index2_annotatedSample = {
                    "screenId":index2ToScreenidMap_screenId,
                    "index2": demultiplexedFiles_INDEX2_STRING,
                    "read1":demultiplexedFiles_index2_IndexPair.read1,
                    "read2":demultiplexedFiles_index2_IndexPair.read2
                }
                Pair[String, Pair[AnnotatedSample, Array[String]]] index2_screenIdAnnotatedSamplePair = (index2ToScreenidMap_screenId, (index2_annotatedSample,index2ToSampleInfoVarsMap_sampleAnnotation))
            }

            Array[Pair[String, Pair[AnnotatedSample, Array[String]]]] index2_screenIdToSampleArray = index2_screenIdAnnotatedSamplePair
            Map[String, Array[Pair[AnnotatedSample, Array[String]]]] index2_screenIdToSampleMap = collect_by_key(index2_screenIdToSampleArray)
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

    Map[String, Array[Pair[AnnotatedSample, Array[String]]]] final_screenIdToSampleMap = select_first([index1_index2_screenIdToSampleMap, index1_screenIdToSampleMap, index2_screenIdToSampleMap, nonDemultiplexed_screenIdToSampleMap])
    
    # Select the final screenIdToSampleMap
    output {
        Map[String, Array[Pair[AnnotatedSample, Array[String]]]] output_screenIdToSampleMap = final_screenIdToSampleMap
    }
}