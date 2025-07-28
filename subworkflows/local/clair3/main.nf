process CLAIR3_DIRECT {

    tag "${meta.id}"
    container "hkubal/clair3:latest"
    label "process_single"

    input:
    val clair3_args
    val meta
    path bam
    path bam_bai
    path ref
    path ref_fai

    output:
    tuple val(meta), path("${meta.id}.vcf.gz"), emit: vcf
    tuple val(meta), path("clair3_report.txt"), emit: report
    path "versions.yml", emit: versions
    
    script:
    clair3_opt_ctg_name = clair3_args.ctg_name ? "--ctg_name=${clair3_args.ctg_name}" : ""
    clair3_opt_platform = clair3_args.platform ? "--platform=${clair3_args.platform}" : ""
    clair3_opt_threads = clair3_args.threads ? "--threads=${clair3_args.threads}" : ""
    clair3_opt_model_path = clair3_args.model_path ? "--model_path=${clair3_args.model_path}" : ""

    """
    /opt/bin/run_clair3.sh \
        --bam_fn=${bam} \
        --ref_fn=${ref} \
        --threads=${clair3_args.threads} \
        --platform=${clair3_args.platform} \
        --model_path=${clair3_args.model_path} \
        --output=clair3_output \
        ${clair3_opt_ctg_name} ${clair3_opt_platform} ${clair3_opt_threads} \
        ${clair3_opt_model_path}

    if [ -f "clair3_output/merge_output.vcf.gz" ]; then
        mv clair3_output/merge_output.vcf.gz ${meta.id}.vcf.gz
    elif [ -f "clair3_output/pileup.vcf.gz" ]; then
        mv clair3_output/pileup.vcf.gz ${meta.id}.vcf.gz
    else
        echo "No VCF output found!" >> clair3_report.txt
        touch ${meta.id}.vcf.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: $(run_clair3.sh --version 2>&1 | head -1 | sed 's/Clair3 //' || echo "1.0.10")
    END_VERSIONS
    """
}
    
   



workflow separated_clair3 {

    take:
        clair3_args
        samples
        reference

    main:
        samples
            .combine(reference)
            .map { meta_bam, ref_tuple ->
                def (meta, bam, bam_bai) = meta_bam
                def (ref, ref_fai) = ref_tuple
                tuple(clair3_args, meta, bam, bam_bai, ref, ref_fai)
            }
            | CLAIR3_DIRECT
            | set {clair3_results}

    emit:
        clair3_results
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

def helpMessage() {
    log.info"""
    Wrapper pipeline around Google GPU DeepVariant run_deepvariant script.

    Usage:
        nextflow run main.nf --bams_dir <path> --ref <path> [options]

    Required Arguments:
        --bams_dir: Aligned, sorted, indexed BAM file containing the reads we
          want to call. Should be aligned to a reference genome compatible with --ref.
        --ref: Genome reference to use. Must have an associated FAI index as
          well. Supports text or gzipped references. Should match the reference used
          to align the BAM file provided to --reads.

    Optional Arguments:
        --call_variants_extra_args: A comma-separated list of flag_name=flag_value.
          "flag_name" has to be valid flags for call_variants.py. If the flag_value is
          boolean, it has to be flag_name=true or flag_name=false.
        --make_examples_extra_args: A comma-separated list of flag_name=flag_value.
          "flag_name" has to be valid flags for make_examples.py. If the flag_value is
          boolean, it has to be flag_name=true or flag_name=false.
        --num_shards: Number of shards for make_examples step.
          (default: '1')
          (an integer)
        --postprocess_variants_extra_args: A comma-separated list of
          flag_name=flag_value. "flag_name" has to be valid flags for
          postprocess_variants.py. If the flag_value is boolean, it has to be
          flag_name=true or flag_name=false.
        --regions: Space-separated list of regions we want to process.
          Elements can be region literals (e.g., chr20:10-20) or paths to BED/BEDPE
          files.""".stripIndent()

}

workflow {

    if (params.help) {
        helpMessage()
        exit 0
    }
    ch_dv_args = channel.from(
        [
            model_type: params.model_type,
            regions: params.regions,
            num_shards: params.num_shards,
            make_examples_extra_args: params.make_examples_extra_args,
            call_variants_extra_args: params.call_variants_extra_args,
            postprocess_variants_extra_args: params.postprocess_variants_extra_args
        ]
    )
    ch_samples = channel.fromPath("${params.bams_dir}/*.bam")
        .map{tuple([id: it.baseName], it, "${it}.bai")}
    ch_ref = channel.fromPath("${params.ref}", checkIfExists: true)
        .map{tuple([id: it.baseName], it, "${it}.fai")}
    separated_deepvariant(ch_dv_args, ch_samples, ch_ref)
}
