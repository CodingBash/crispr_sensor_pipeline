version development


struct IndexPair {
    String index
    File read1
    File? read2
}

struct UndeterminedFiles {
    File undeterminedR1File
    File? undeterminedR2File
}

struct DemultiplexedFiles {
    Map[String, IndexPair] demultiplexedFiles
}

workflow BBMapDemultiplexRunnerPostprocessWorkflow {
    input {
        Array[File] demultiplexedOutputR1Files
        Array[File]? demultiplexedOutputR2Files
        File undeterminedR1File
        File? undeterminedR2File
        String rootDemultiplexedOutputFilenamePattern
    }


    # Create a Map of the index to R1 files
    scatter(demultiplexedOutputR1File in demultiplexedOutputR1Files){
        String demultiplexedOutputR1File_index = basename(demultiplexedOutputR1File,sub(rootDemultiplexedOutputFilenamePattern, "#", "1"))
        Pair[String, File] demultiplexedOutputR1FileIndexPair = (demultiplexedOutputR1File_index, demultiplexedOutputR1File)
    }

    Map[String, File] demultiplexedOutputR1FileMap = as_map(demultiplexedOutputR1FileIndexPair)



    # Create a Map of the index to R2 files (if exists)
    if (defined(demultiplexedOutputR2Files)){
        scatter(demultiplexedOutputR2File in select_first([demultiplexedOutputR2Files])){
            String demultiplexedOutputR2File_index = basename(demultiplexedOutputR2File, sub(rootDemultiplexedOutputFilenamePattern, "#", "2"))
            Pair[String, File] demultiplexedOutputR2FileIndexPair = (demultiplexedOutputR2File_index, demultiplexedOutputR2File)
        }
        Map[String, File] demultiplexedOutputR2FileMap = as_map(demultiplexedOutputR2FileIndexPair)
    }


    # Create a pairing of the index to the R1 and R2 files as an IndexPair Struct
    if (defined(demultiplexedOutputR2FileMap)){
        Map[String,File] demultiplexedOutputR2FileMap_defined = select_first([demultiplexedOutputR2FileMap])

        scatter(demultiplexedOutputR1FilePair_PE in as_pairs(demultiplexedOutputR1FileMap)){
            String demultiplexedOutputR1FilePair_PE_index = demultiplexedOutputR1FilePair_PE.left
            File demultiplexedOutputR1FilePair_PE_r1 = demultiplexedOutputR1FileMap[demultiplexedOutputR1FilePair_PE_index]
            File demultiplexedOutputR1FilePair_PE_r2 = demultiplexedOutputR2FileMap_defined[demultiplexedOutputR1FilePair_PE_index]
            IndexPair demultiplexedOutputR1FilePair_PE_readIndex = {
                "index":demultiplexedOutputR1FilePair_PE_index,
                "read1":demultiplexedOutputR1FilePair_PE_r1,
                "read2":demultiplexedOutputR1FilePair_PE_r2
            }
            Pair[String, IndexPair] readIndexPair_PE = (demultiplexedOutputR1FilePair_PE_index, demultiplexedOutputR1FilePair_PE_readIndex)
        }
    }

    # Create a pairing of the index to the R1 files as an IndexPair Struct (if R2 not exists)
    if (!defined(demultiplexedOutputR2FileMap)){
            scatter(demultiplexedOutputR1FilePair_SE in as_pairs(demultiplexedOutputR1FileMap)){
                String demultiplexedOutputR1FilePair_SE_index = demultiplexedOutputR1FilePair_SE.left
                File demultiplexedOutputR1FilePair_SE_r1 = demultiplexedOutputR1FileMap[demultiplexedOutputR1FilePair_SE_index]
                IndexPair demultiplexedOutputR1FilePair_SE_readIndex = {
                    "index":demultiplexedOutputR1FilePair_SE_index,
                    "read1":demultiplexedOutputR1FilePair_SE_r1
                }
                Pair[String, IndexPair] readIndexPair_SE = (demultiplexedOutputR1FilePair_SE_index, demultiplexedOutputR1FilePair_SE_readIndex)
            }
    }

    # Coerce Pairs to Maps
    Array[Pair[String, IndexPair]] readIndexPairList = select_first([readIndexPair_PE, readIndexPair_SE])
    Map[String, IndexPair] readIndexMap = as_map(readIndexPairList)
    
    # Create final object
    UndeterminedFiles undeterminedFiles = {
            "undeterminedR1File": undeterminedR1File,
            "undeterminedR2File": undeterminedR2File
    }
    DemultiplexedFiles demultiplexedFiles = {
            "demultiplexedFiles": readIndexMap
    }
    
    Pair[DemultiplexedFiles, UndeterminedFiles] demultiplexedResults = (demultiplexedFiles, undeterminedFiles)
   
    output {
         Pair[DemultiplexedFiles, UndeterminedFiles] output_demultiplexedResults = demultiplexedResults
    }
   
 }