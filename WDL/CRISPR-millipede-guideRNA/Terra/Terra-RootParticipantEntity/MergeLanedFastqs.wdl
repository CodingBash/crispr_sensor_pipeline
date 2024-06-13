version development
workflow MergeLanedFastqs {
    input {
        Array[File] individualRead1Fastqs
        Array[File]? individualRead2Fastqs

        # Allowing ability to concat using zcat instead of cat
        Boolean uncompressConcatenation = false

         # Constants for calculating resources
        Int estimatedeMemoryBuffer = 2
        Int estimatedDiskGbBuffer = 10
        Int estimatedDockerGbSize = 2

        # Avoid exorbitant resources by setting max resource limit
        Int maxDiskSpaceGB = 500
        Int maxMemoryGB = 10
        Int maxCores = 5

        # Optionally provide explicit resource amounts
        Int? diskGB
        Int? memoryGB
        Int? cpus

        # Provide other resource specifications
        Int preemptible = 1
        Int maxRetries = 0
        String diskType = "HDD"
    }

    Int estimatedR1DiskGB = ceil(estimatedDockerGbSize + (size(individualRead1Fastqs, "GB")*2) + estimatedDiskGbBuffer)
    Int specifiedR1DiskGB = min(select_first([diskGB, estimatedR1DiskGB]), maxDiskSpaceGB)

    Int estimatedMemoryGB = ceil(2 + estimatedeMemoryBuffer)
    Int specifiedMemoryGB = min(select_first([memoryGB, estimatedMemoryGB]), maxMemoryGB)
    Int estimatedCoresNeeded = 1
    Int specifiedCores = min(select_first([cpus, estimatedCoresNeeded]), maxCores)

    call MergeLanedFastqsTask as MergeLanedFastqsTask_R1 {
        input:
            individualReadFastqs=individualRead1Fastqs,
            uncompressConcatenation=uncompressConcatenation,
            preemptible=preemptible,
            maxRetries=maxRetries,
            specifiedMemoryGB=specifiedMemoryGB,
            specifiedDiskGB=specifiedR1DiskGB,
            diskType=diskType,
            specifiedCores=specifiedCores
    }

    if(defined(individualRead2Fastqs)){
        Array[File] individualRead2Fastqs_defined = select_first([individualRead2Fastqs])
        
        Int estimatedR2DiskGB = ceil(estimatedDockerGbSize + (size(individualRead2Fastqs_defined, "GB")*2) + estimatedDiskGbBuffer)
        Int specifiedR2DiskGB = min(select_first([diskGB, estimatedR2DiskGB]), maxDiskSpaceGB)

        call MergeLanedFastqsTask as MergeLanedFastqsTask_R2 {
            input:
                individualReadFastqs=individualRead2Fastqs_defined,
                uncompressConcatenation=uncompressConcatenation,
                preemptible=preemptible,
                maxRetries=maxRetries,
                specifiedMemoryGB=specifiedMemoryGB,
                specifiedDiskGB=specifiedR2DiskGB,
                diskType=diskType,
                specifiedCores=specifiedCores
        }
    }
    output {
        File concatenatedRead1Fastq = MergeLanedFastqsTask_R1.concatenatedReadFastq
        File? concatenatedRead2Fastq = MergeLanedFastqsTask_R2.concatenatedReadFastq
    }
}

task MergeLanedFastqsTask {
    input {
        Array[File] individualReadFastqs

        # Allowing ability to concat using zcat instead of cat
        Boolean uncompressConcatenation = false

        #Specify resources
        Int preemptible
        Int maxRetries
        Int specifiedMemoryGB
        Int specifiedDiskGB
        String diskType
        Int specifiedCores
    }

    command <<<
        if [ ~{uncompressConcatenation} ]; then
            zcat ~{sep(" ", individualReadFastqs)} | gzip > concat.fastq.gz
        else
            cat ~{sep(" ", individualReadFastqs)} | gzip > concat.fastq.gz
        fi
    >>>

    output {
        File concatenatedReadFastq = "concat.fastq.gz"
    }

    runtime {
    	docker: "ubuntu:latest"
        preemptible: "${preemptible}"
        maxRetries: "${maxRetries}"
        memory: "${specifiedMemoryGB} GB"
        disks: "local-disk ${specifiedDiskGB} ${diskType}"
        cpu: "${specifiedCores}"
    }
}


