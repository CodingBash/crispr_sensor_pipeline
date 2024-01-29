version development

task UmiToolsExtractTask {
    input {
        # TASK PARAMS
        File inputRead1
        File? inputRead2
        String r1UmiToolsParseRegex
        String? r2UmiToolsParseRegex
        String sampleName
        
        #RUNTIME INPUT PARAMS
	    Int preemptible = 1
	    Int diskGB = 10
        Int memoryGB = 2
        Int maxRetries = 0
        String diskType = "SSD"
    }

    # SET OUTPUT FILENAMES
	String logFileFn="${sampleName}.umitools_log.txt"
	String outputRead1Fn="${sampleName}.umitools.R1.fastq"
	String outputRead2Fn="${sampleName}.umitools.R2.fastq"
	String outputFilteredRead1Fn="${sampleName}.umitools.filtered.R1.fastq"
	String outputFilteredRead2Fn="${sampleName}.umitools.filtered.R2.fastq"
	
    # SET BOOLEANS TO DETERMINE R1-only vs. R1+R2 RUN
    Boolean inputRead2Defined = defined(inputRead2)
    Boolean r2UmiToolsParseRegexDefined = defined(r2UmiToolsParseRegex)

    # EXECUTE UMI-TOOLS PARRSING
	command <<<
        if [ ~{inputRead2Defined} ] && [ ~{r2UmiToolsParseRegexDefined} ]; then
            echo "Running R1+R2 Mode"
            umi_tools extract --extract-method=regex \
            --stdin ~{inputRead1} \
            --read2-in ~{inputRead2} \
            --bc-pattern ~{r1UmiToolsParseRegex} \
            --bc-pattern2 ~{r2UmiToolsParseRegex} \
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
            --bc-pattern=~{r1UmiToolsParseRegex} \
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
		preemptible: "${preemptible}"
        maxRetries: "${maxRetries}"
		memory: "${memoryGB} GB"
		disks: "local-disk ${diskGB} ${diskType}"
        cpu: 1
        
	}
}