nextflow.enable.dsl=2

params.crams = false
params.reference = false

if(!params.crams) {
    exit 1, "--crams argument like '/path/to/*.cram' is required"
}
if(!params.reference) {
    exit 1, "--reference argument is required"
}

process strling_extract {
    input:
    tuple val(sample), path(cram), path(crai)
    path(reference)
    path(fai)
    val(proportion_repeat)
    val(min_mapq)

    output:
    tuple val("${sample}"), path(cram), path(crai), path("${sample}.bin"), emit: bin
    path("${sample}.bin"), emit: bin_only

    script:
    """
    strling extract -f $reference -g ${params.str_index} -p $proportion_repeat -q $min_mapq $cram ${sample}.bin
    """
}

process strling_merge {
    input:
    path(bin)
    path(reference)
    path(fai)
    val(window)
    val(min_support)
    val(min_clip)
    val(min_clip_total)
    val(min_mapq)

    output:
    path("joint-bounds.txt"), emit: bounds

    script:
    """
    strling merge -f $reference -w $window -m $min_support -c $min_clip \
        -t $min_clip_total -q $min_mapq -o joint $bin
    """
}

process strling_call {
    publishDir "${params.strling_call_outdir}/"${sample}_outliers", mode: 'link'

    input:
    tuple val(sample), path(cram), path(crai), path(bin)
    path(reference)
    path(fai)
    path(bounds)
    val(min_mapq)
    val(min_support)
    val(min_clip)
    val(min_clip_total)

    output:
    path("${sample}-bounds.txt"), emit: bounds
    path("${sample}-genotype.txt"), emit: genotypes
    path("${sample}-unplaced.txt"), emit: unplaced

    script:
    b = bounds ? "-b $bounds" : ""
    """
    strling call -o $sample $b -m $min_support -c $min_clip -t $min_clip_total \
        -q $min_mapq -v -f $reference $cram $bin
    """
}

process strling_outliers {
    publishDir "${params.outliers_outdir}/outliers", mode: 'link'

    input:
    path(genotypes)
    path(unplaced)
    path(control)
    val(slop)
    val(min_clips)
    val(min_size)

    output:
    path("*STRs.tsv")
    path("control.tsv")
    path("depths.tsv")
    path("unplaced.tsv")

    script:
    c = control ? "--control $control" : ""
    """
    strling-outliers.py --genotypes $genotypes --unplaced $unplaced \
        --emit control.tsv --slop $slop --min_clips $min_clips \
        --min_size $min_size $c
    """
}

workflow {
    crams = Channel.fromPath(params.crams, checkIfExists: true)
        .map { file -> tuple(file.simpleName, file, file + ("${file}".endsWith('.cram') ? '.crai' : '.bai')) }
    fai = "${params.reference}.fai"

    strling_extract(
        crams,
        params.reference,
        fai,
        params.proportion_repeat,
        params.min_mapq
    )
    strling_merge(
        strling_extract.out.bin_only.collect(),
        params.reference,
        fai,
        params.window,
        params.min_support,
        params.min_clip,
        params.min_clip_total,
        params.min_mapq
    )
    strling_call(
        strling_extract.out.bin,
        params.reference,
        fai,
        strling_merge.out.bounds,
        params.min_mapq,
        params.min_support,
        params.min_clip,
        params.min_clip_total
    )
    strling_outliers(
        strling_call.out.genotypes.collect(),
        strling_call.out.unplaced.collect(),
        params.control,
        params.slop,
        params.min_clips,
        params.min_size
    )
}
