version development

import "https://api.firecloud.org/ga4gh/v1/tools/pinellolab:UmiToolsExtractTask/versions/3/plain-WDL/descriptor" as umi_tools_extract
import "https://api.firecloud.org/ga4gh/v1/tools/pinellolab:BBMapDemultiplexOrchestratorWorkflow/versions/6/plain-WDL/descriptor" as demultiplex

workflow CrisprMillipedeGuideDemultiplexParticipantEntity {
    input {
        # UMI TOOLS INPUT
        File? rawFastqR1 
        File? rawFastqR2
        String? r1UmiToolsParseRegex # TODO: Rename to distinguish between the other REGEX
        String? r2UmiToolsParseRegex # TODO: Rename to distinguish between the other REGEX
        String sampleName

        Int umiToolsExtractTaskPreemptible = 1
	    Int umiToolsExtractTaskDiskGB = 10
        Int umiToolsExtractTaskMemoryGB = 2
        Int umiToolsExtractTaskMaxRetries = 0
        String umiToolsExtractTaskDiskType = "SSD"
        
        File? umiToolsFastqR1
		File? umiToolsFastqR2


        #
        #   FOR DEMULTIPLEX
        #
        File? index1StringTextFile
        Array[String]? index1List
        File? index2StringTextFile
        Array[String]? index2List
        
        Int index1Hamming 
        Int index2Hamming

        String index1BBMapHeaderDelimiter = "\\\\+"
        Int index1BBMapHeaderColumn = 2

        String index2BBMapHeaderDelimiter = "_"
        Int index2BBMapHeaderColumn = 2

        #NEW contains the demultiplex to screen ID map
        String? input_screenId
        Map[String, String]? input_index1ToScreenidMap
        Map[String, String]? input_index2ToScreenidMap
        Map[String, Map[String, String]]? input_index1ToIndex2ToScreenidMap

        #NEW contains the demultiplex to sample Information
        Array[String]? input_sampleInfoVarnames
        Array[String]? input_sampleInfoVars
        Map[String, Array[String]]? input_index1ToSampleInfoVarsMap
        Map[String, Array[String]]? input_index2ToSampleInfoVarsMap
        Map[String, Map[String, Array[String]]]? input_index1ToIndex2ToSampleInfoVarsMap

        Int demultiplexTaskPreemptible = 1
        Int demultiplexTaskDiskGB = 10
        Int demultiplexTaskMemoryGB = 2
        Int demultiplexTaskMaxRetries = 0
        String demultiplexTaskDiskType = "SSD"
    }

    # Extract relevant sequences using UMI tools
    # TODO: Once we have the seqspec tool, perhaps that plus the extract task can be moddularized into a new workflow.
    if(!defined(umiToolsFastqR1)){
        call umi_tools_extract.UmiToolsExtractTask {
            input:
                inputRead1=select_first([rawFastqR1]),
                inputRead2=rawFastqR2,
                r1UmiToolsParseRegex=select_first([r1UmiToolsParseRegex]),
                r2UmiToolsParseRegex=r2UmiToolsParseRegex,
                sampleName=sampleName,
                preemptible=umiToolsExtractTaskPreemptible,
                diskGB=umiToolsExtractTaskDiskGB,
                memoryGB=umiToolsExtractTaskMemoryGB,
                maxRetries=umiToolsExtractTaskMaxRetries,
                diskType=umiToolsExtractTaskDiskType
        }
    
        # Demultiplex the UMI tools result
        call demultiplex.BBMapDemultiplexOrchestratorWorkflow as demultiplexWorkflowUmiToolsPerformed {
            input:
                inputRead1=UmiToolsExtractTask.outputRead1,
                inputRead2=UmiToolsExtractTask.outputRead2,
                index1StringTextFile=index1StringTextFile,
                index1List=index1List,
                index2StringTextFile=index2StringTextFile,
                index2List=index2List,
                index1Hamming=index1Hamming,
                index2Hamming=index2Hamming,
                sampleName=sampleName,
                index1BBMapHeaderDelimiter=index1BBMapHeaderDelimiter,
                index1BBMapHeaderColumn=index1BBMapHeaderColumn,
                index2BBMapHeaderDelimiter=index2BBMapHeaderDelimiter,
                index2BBMapHeaderColumn=index2BBMapHeaderColumn,

                # Mapping from indices to screen name
                screenId=input_screenId,
                index2ToScreenidMap=input_index1ToScreenidMap,
                index2ToScreenidMap=input_index2ToScreenidMap,
                index1ToIndex2ToScreenidMap=input_index1ToIndex2ToScreenidMap,

                # Mapping from indices to sample annotations
                sampleInfoVars=input_sampleInfoVars,
                index1ToSampleInfoVarsMap=input_index1ToSampleInfoVarsMap,
                index2ToSampleInfoVarsMap=input_index2ToSampleInfoVarsMap,
                index1ToIndex2ToSampleInfoVarsMap=input_index1ToIndex2ToSampleInfoVarsMap,

                # Runtime param
                preemptible=demultiplexTaskPreemptible,
                diskGB=demultiplexTaskDiskGB,
                memoryGB=demultiplexTaskMemoryGB,
                maxRetries=demultiplexTaskMaxRetries,
                diskType=demultiplexTaskDiskType
        }
    } 
    if(defined(umiToolsFastqR1)){
        call demultiplex.BBMapDemultiplexOrchestratorWorkflow as demultiplexWorkflowUmiToolsSkipped {
            input:
                inputRead1=select_first([umiToolsFastqR1]),
                inputRead2=umiToolsFastqR2,
                index1StringTextFile=index1StringTextFile,
                index1List=index1List,
                index2StringTextFile=index2StringTextFile,
                index2List=index2List,
                index1Hamming=index1Hamming,
                index2Hamming=index2Hamming,
                sampleName=sampleName,
                index1BBMapHeaderDelimiter=index1BBMapHeaderDelimiter,
                index1BBMapHeaderColumn=index1BBMapHeaderColumn,
                index2BBMapHeaderDelimiter=index2BBMapHeaderDelimiter,
                index2BBMapHeaderColumn=index2BBMapHeaderColumn,

                # Mapping from indices to screen name
                screenId=input_screenId,
                index1ToScreenidMap=input_index1ToScreenidMap,
                index2ToScreenidMap=input_index2ToScreenidMap,
                index1ToIndex2ToScreenidMap=input_index1ToIndex2ToScreenidMap,

                # Mapping from indices to sample annotations
                sampleInfoVars=input_sampleInfoVars,
                index1ToSampleInfoVarsMap=input_index1ToSampleInfoVarsMap,
                index2ToSampleInfoVarsMap=input_index2ToSampleInfoVarsMap,
                index1ToIndex2ToSampleInfoVarsMap=input_index1ToIndex2ToSampleInfoVarsMap,

                # Runtime params
                preemptible=demultiplexTaskPreemptible,
                diskGB=demultiplexTaskDiskGB,
                memoryGB=demultiplexTaskMemoryGB,
                maxRetries=demultiplexTaskMaxRetries,
                diskType=demultiplexTaskDiskType
        }
    }

    Map[String, Array[Pair[AnnotatedSample, Array[String]]]] output_screenIdToSampleMap_final = select_first([demultiplexWorkflowUmiToolsSkipped.output_screenIdToSampleMap, demultiplexWorkflowUmiToolsPerformed.output_screenIdToSampleMap])

    output {
        #
        #   UMI Tools Outputs
        #
        File? umiToolsExtractedOutputRead1 = UmiToolsExtractTask.outputRead1
        File? umiToolsExtractedOutputRead2 = UmiToolsExtractTask.outputRead2
        File? umiToolsExtractedOutputFilteredRead1=UmiToolsExtractTask.outputFilteredRead1
        File? umiToolsExtractedOutputFilteredRead2=UmiToolsExtractTask.outputFilteredRead2
        File? umiToolsExtractTaskLogFile=UmiToolsExtractTask.logFile
        Float? umiToolsExtractOutputRead1ReadCount = UmiToolsExtractTask.outputRead1ReadCount
        Float? umiToolsExtractOutputFilteredRead1ReadCount = UmiToolsExtractTask.outputFilteredRead1ReadCount

        #
        #   Demultiplexing Outputs
        #
        Map[String, Array[Pair[AnnotatedSample, Array[String]]]] output_screenIdToSampleMap = output_screenIdToSampleMap_final
    }
}


