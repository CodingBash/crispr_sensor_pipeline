version development

struct CrispressoSampleSettings {
    String sampleId

    # Metadata for downstream modelling purposes
    Array[String] sampleMetadata

    # Sample mode to run in
    String? guideID # if empty, then will run in pooled-mode, else, will run in single-mode
    
    # Sequence files
    File? premerged_FASTQ
    File? R1_FASTQ
    File? R2_FASTQ

    Boolean isUnedited
}

struct CrispressoGuideInformation {
    String id
    String protospacerTarget
    Int coordinate
}


task Crispresso2BatchPooledTask {
    input {
        String sampleName
    }

    Boolean r24ntBarcodeExtractionRegexDefined = defined(r24ntBarcodeExtractionRegex)

	command <<<
    >>>

    output {
		Float outputFilteredRead1ReadCount=read_float('output_filtered_read_count.txt')
	}

	runtime {
		docker: "pinellolab/crispresso2:v2.2.14"
	}
}

task Crispresso2BatchSingleGuideTask {
    input {
        File settingsFile
        String inputAmpliconSequence
        String ampliconName
        String outputName
        Int qualityThreshold
        Int excludeBpFromLeft
        Int excludeBpFromRight
        Float minFrequencyAllelesAroundCutToPlot
        Int maxRowsAllelesAroundCutToPlot
        
        Int singleGuideWindowSize
        Int quantificationWindow

        Int coresPerRun

    }

    Boolean r24ntBarcodeExtractionRegexDefined = defined(r24ntBarcodeExtractionRegex)

	command <<<
        # TODO: May be too hard to parse the outputs without some custom post-processing code to convert output files to JSON. Might need to run in single sample mode, or retrieve per-sample names as list and parse, likely single sample mode since unsure how to do the latter.
        # Actually I may be able to parse outputs using some glob and sub action (similar to DemultiplexPostprocessWorkflow)
        !CRISPRessoBatch -bs ~{settingsFile} -a ~{inputAmpliconSequence} \
        -an ~{ampliconName} -q ~{qualityThreshold} \
        --exclude_bp_from_left ~{excludeBpFromLeft} \
        --exclude_bp_from_right ~{excludeBpFromRight} --no_rerun -n ~{outputName} \
        --min_frequency_alleles_around_cut_to_plot ~{minFrequencyAllelesAroundCutToPlot} --max_rows_alleles_around_cut_to_plot ~{maxRowsAllelesAroundCutToPlot} -p {coresPerRun}  \
        --plot_window_size ~{singleGuideWindowSize} --base_editor_output -w {quantificationWindow} -bo "output/"

        # ZIP output
    >>>

    output {
        File completeOutputZip = "output.zip"
		File = "output/CRISPRessoBatch_on_${outputName}/....."
	}

	runtime {
		docker: "pinellolab/crispresso2:v2.2.14"
	}
}


workflow CrisprMillipedeDirect_Workflow {
    input {
        # Inputs for Crispresso2Batch Settings File
        Array[String] sampleMetadataNames
        Array[CrispressoSampleSettings] crispressoBatchSampleSettingsList
        Array[CrispressoGuideInformation] guideList

        # Crispresso2Batch Command
        String fullAmpliconSequence
        String? read1AmpliconSequence # Must provide if run_singleEndFASTQ and read1 FASTQs are provided
        String? read2AmpliconSequence # Must provide if run_singleEndFASTQ and read2 FASTQs are provided
        String? psuedoMiddleGuideProtospacer # Must provide if any pooled samples are provided

        Int? qualityThreshold = 30
        Int? excludeBpFromLeft = 3
        Int? excludeBpFromRight = 3
        Float? minFrequencyAllelesAroundCutToPlot = 0.001
        Int? maxRowsAllelesAroundCutToPlot = 500
        Int? coresPerRun = 1

        Int? singleGuideWindowSize = 4
        Int? quantificationWindow = 0

        Int read1Length
        Int read2Length

        Boolean? run_premergedFASTQ = true
        Boolean? run_singleEndFASTQ = true
    }

    #
    # RUN POOLED MODE ON POOLED SAMPLES
    #

    # Get the set of pooled samples 
    scatter(crispressoBatchSampleSettings in crispressoBatchSampleSettingsList){
        if(!defined(crispressoBatchSampleSettings.guideID)){
            Pair[String, CrispressoSampleSettings] pooledCrispressoBatchSampleSettings = (crispressoBatchSampleSettings.sampleId, crispressoBatchSampleSettings)    
        }
    }
    Map[String, CrispressoSampleSettings] pooledCrispressoBatchSampleSettingsList = as_map(select_all(pooledCrispressoBatchSampleSettings))

    # Generate settings file

    #
    # RUN SINGLE MODE ON SINGLE GUIDE SAMPLES (USING PSUEDO-MIDDLE PROTOSPACER)
    #
    scatter(crispressoBatchSampleSettings in crispressoBatchSampleSettingsList){
        if(defined(crispressoBatchSampleSettings.guideID)){
            Pair[String, CrispressoSampleSettings] sgCrispressoBatchSampleSettings = (crispressoBatchSampleSettings.sampleId, crispressoBatchSampleSettings)    
        }
    }
    Map[String, CrispressoSampleSettings] sgCrispressoBatchSampleSettingsList = as_map(select_all(sgCrispressoBatchSampleSettings))

    # Generate settings file


    output {
    }
}


