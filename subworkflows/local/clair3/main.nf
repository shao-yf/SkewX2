process CLAIR3_CALL {
    
    tag "${meta.id}"
    container "hkubal/clair3:latest"
    label "process_high"
    publishDir "${params.outdir}", pattern: "*.html", mode: "copy" // only save the report
    cpus   = "${params.threads ? params.threads : 1}"
    memory "${params.threads} GB" //actually threads/4 ?
      
    input:
    each clair3_args
    tuple   val(meta), 
            path(bam), 
            path(bam_bai), 
            path(ref), 
            path(ref_fai)
    
    output:
    tuple   val(meta), 
            path(bam), 
            path(bam_bai), 
            path(ref), 
            path(ref_fai),  
            path("${meta.id}.vcf.gz"), 
            path("${meta.id}.vcf.gz.tbi")
    
    script:
    // declare optional flags if present in clair3_args
    clair3_opt_ctg_name = clair3_args.ctg_name ? "--ctg_name=${clair3_args.ctg_name}" : ""
    
    //--output=${OUTPUT_DIR}/skew_test_results not define the output dir & no html_report generated
    """
    # Run Clair3
    /opt/bin/run_clair3.sh \\
        --bam_fn=${bam} \\
        --ref_fn=${ref} \\
        --threads=${clair3_args.threads} \\
        --platform=${clair3_args.platform} \\
        --model_path=${clair3_args.model_path} \\
        --output=clair3_output \\
        ${clair3_opt_ctg_name}
    
    # Check for output VCF and move to expected location
    if [ -f "clair3_output/merge_output.vcf.gz" ]; then
        mv clair3_output/merge_output.vcf.gz ${meta.id}.vcf.gz
        mv clair3_output/merge_output.vcf.gz.tbi ${meta.id}.vcf.gz.tbi
    elif [ -f "clair3_output/pileup.vcf.gz" ]; then
        mv clair3_output/pileup.vcf.gz ${meta.id}.vcf.gz
        mv clair3_output/pileup.vcf.gz.tbi ${meta.id}.vcf.gz.tbi
    else
        echo "No VCF output found!" >> clair3_report.txt
        touch ${meta.id}.vcf.gz
        touch ${meta.id}.vcf.gz.tbi
    fi
    
    # Create a simple report
    echo "Clair3 variant calling completed for ${meta.id}" > clair3_report.txt
    echo "Platform: ${clair3_args.platform}" >> clair3_report.txt
    echo "Model: ${clair3_args.model_path}" >> clair3_report.txt
    echo "Threads: ${clair3_args.threads}" >> clair3_report.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(/opt/bin/run_clair3.sh --version 2>&1 | head -1 | sed 's/Clair3 //' || echo "1.0.10")
    END_VERSIONS
    """
}

workflow separated_clair3 {
    
    take:
        clair3_args    // Channel containing Clair3 parameters
        reads        // Channel containing sample BAMs with metadata
        ref      // Channel containing reference genome
    
    main:
        // Combine channels and feed into Clair3
        reads
            .combine(ref.map{tuple(it[1], it[2])})
            .set {ch_clair3_input}

        results = CLAIR3_CALL(clair3_args, ch_clair3_input)
    emit:
        results
}

// set params to a falsey default value
params.model_type = ''
params.regions = ''
params.num_shards = ''
params.make_examples_extra_args = ''
params.call_variants_extra_args = ''
params.postprocess_variants_extra_args = ''
params.bams_dir = ''
params.ref = ''
params.ctg_name = ''
params.platform = ''
params.threads = ''
params.model_path = ''

//haven't made a change yet
def helpMessage() {
    log.info"""
    Wrapper pipeline around Clair3 variant calling.

    Usage:
        nextflow run main.nf --bams_dir <path> --ref <path> [options]

    Required Arguments:
        --bams_dir: Aligned, sorted, indexed BAM file containing the reads we
          want to call. Should be aligned to a reference genome compatible with --ref.
        --ref: Genome reference to use. Must have an associated FAI index as
          well. Supports text or gzipped references. Should match the reference used
          to align the BAM file provided to --reads.

    Optional Arguments:
        --ctg_name: Contig name for Clair3 variant calling (default: chrX)
        --platform: Platform for Clair3 variant calling. Possible options: ont, hifi, ilmn (default: ont)
        --threads: Number of threads for Clair3 variant calling (default: 4)
        --model_path: Model path for Clair3 variant calling (default: /opt/models/r1041_e82_400bps_sup_v500)
    """.stripIndent()
}

workflow {

    if (params.help) {
        helpMessage()
        exit 0
    }
    ch_clair3_args = channel.from([
        ctg_name: params.ctg_name,
        platform: params.platform,
        threads: params.threads,
        model_path: params.model_path
    ])
    ch_samples = channel.fromPath("${params.bams_dir}/*.bam")
        .map{tuple([id: it.baseName], it, "${it}.bai")}
    ch_ref = channel.fromPath("${params.ref}", checkIfExists: true)
        .map{tuple([id: it.baseName], it, "${it}.fai")}
    separated_clair3(ch_clair3_args, ch_samples, ch_ref)
}
