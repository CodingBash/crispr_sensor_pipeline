version development

task BBMapDemultiplexRunnerTask {
        input {
            # TASK PARAMS
            File processedRead1
            File? processedRead2
            Array[String] indexStringList
            String delimiter
            Int column
            Int hamming
            String demultiplexedOutputFilenamePattern
            String demultiplexedUndeterminedFilenamePattern

            # RUNTIME PARAMS
            Int preemptible = 1
            Int diskGB = 10
            Int memoryGB = 2
            Int maxRetries = 0
            String diskType = "SSD"
        }

        # BOOLEAN TO DETERMINE WHETHER TO RUN IN R1-only or R1+R2 MODE
        Boolean inputRead2Defined = defined(processedRead2)
        
        # EXECUTE BBMAP
        command <<<
            if [ ~{inputRead2Defined} ]; then
                echo "Running R1+R2 Mode"

                demuxbyname.sh in=~{processedRead1} \
                    in2=~{processedRead2} \
                    out=~{demultiplexedOutputFilenamePattern} \
                    outu=~{demultiplexedUndeterminedFilenamePattern} \
                    delimiter='~{delimiter}' \
                    column=~{column} \
                    names=~{sep(",", indexStringList)} \
                    hdist=~{hamming}

                echo "Completed"
            else
                echo "Running R1-only Mode"

                demuxbyname.sh in=~{processedRead1} \
                    out=~{sub(demultiplexedOutputFilenamePattern, '#', '1')} \
                    outu=~{sub(demultiplexedUndeterminedFilenamePattern, '#', '1')} \
                    delimiter='~{delimiter}' \
                    column=~{column} \
                    names=~{sep(",", indexStringList)} \
                    hdist=~{hamming}

                echo "Completed"
            fi
        >>>
        

        output {
            Array[File]? demultiplexedOutputR1Files = glob(sub(sub(demultiplexedOutputFilenamePattern,"%", "*"), "#", "1"))
            Array[File]? demultiplexedOutputR2Files = glob(sub(sub(demultiplexedOutputFilenamePattern,"%", "*"),"#", "2"))
            File? undeterminedOutputR1File = "${sub(demultiplexedUndeterminedFilenamePattern, "#", "1")}"
            File? undeterminedOutputR2File = "${sub(demultiplexedUndeterminedFilenamePattern, "#", "2")}"
        }

        runtime {
                docker: "quay.io/biocontainers/bbmap:39.01--h92535d8_1"
                preemptible: "${preemptible}"
                maxRetries: "${maxRetries}"
                memory: "${memoryGB} GB"
                disks: "local-disk ${diskGB} ${diskType}"
                cpu: 1
        }
}