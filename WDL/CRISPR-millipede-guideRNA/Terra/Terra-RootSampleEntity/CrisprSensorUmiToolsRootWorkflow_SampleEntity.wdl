version development

import "https://api.firecloud.org/ga4gh/v1/tools/pinellolab:UmiToolsExtractTask/versions/3/plain-WDL/descriptor" as umi_tools_extract

workflow CrisprMillipedeGuideUmiToolsOrchestratorWorkflowSampleEntity {
    input {
        # UMI TOOLS INPUT
        File rawFastqR1 
        File? rawFastqR2
        String r1UmiToolsParseRegex # TODO: Rename to distinguish between the other REGEX
        String? r2UmiToolsParseRegex # TODO: Rename to distinguish between the other REGEX
        String sampleName

        Int umiToolsExtractTaskPreemptible = 1
	    Int umiToolsExtractTaskDiskGB = 10
        Int umiToolsExtractTaskMemoryGB = 2
        Int umiToolsExtractTaskMaxRetries = 0
        String umiToolsExtractTaskDiskType = "SSD"
    }

    # Extract relevant sequences using UMI tools
    # TODO: Once we have the seqspec tool, perhaps that plus the extract task can be moddularized into a new workflow.
    call umi_tools_extract.UmiToolsExtractTask {
        input:
            inputRead1=rawFastqR1,
            inputRead2=rawFastqR2,
            r1UmiToolsParseRegex=r1UmiToolsParseRegex,
            r2UmiToolsParseRegex=r2UmiToolsParseRegex,
            sampleName=sampleName,
            preemptible=umiToolsExtractTaskPreemptible,
            diskGB=umiToolsExtractTaskDiskGB,
            memoryGB=umiToolsExtractTaskMemoryGB,
            maxRetries=umiToolsExtractTaskMaxRetries,
            diskType=umiToolsExtractTaskDiskType
    }
    
    output {
        #
        #   UMI Tools Outputs
        #
        File? umiToolsExtractedOutputRead1 = UmiToolsExtractTask.outputRead1
        File? umiToolsExtractedOutputRead2 = UmiToolsExtractTask.outputRead2
        File? umiToolsExtractedOutputFilteredRead1 = UmiToolsExtractTask.outputFilteredRead1
        File? umiToolsExtractedOutputFilteredRead2 = UmiToolsExtractTask.outputFilteredRead2
        File? umiToolsExtractTaskLogFile = UmiToolsExtractTask.logFile
        Float? umiToolsExtractOutputRead1ReadCount = UmiToolsExtractTask.outputRead1ReadCount
        Float? umiToolsExtractOutputFilteredRead1ReadCount = UmiToolsExtractTask.outputFilteredRead1ReadCount
    }
}


