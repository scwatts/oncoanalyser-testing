process CUSTOM_EXTRACTCONTIG {
    tag "${contig_name}"
    label 'process_single'

    conda "bwa-mem2=2.2.1 samtools=1.19.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-4dde50190ae599f2bb2027cb2c8763ea00fb5084:544519c4a0ff7e9616a3b44afde1f143c52f10c3-0' :
        'biocontainers/mulled-v2-4dde50190ae599f2bb2027cb2c8763ea00fb5084:544519c4a0ff7e9616a3b44afde1f143c52f10c3-0' }"

    input:
    val contig_name
    path genome_fasta
    path genome_fai
    val run

    output:
    path "*extracted.fa"  , emit: contig
    path "*extracted.fa.*", emit: bwa_index
    path 'versions.yml'   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    samtools faidx \\
        -o ${contig_name}_extracted.fa \\
        ${genome_fasta} \\
        ${contig_name}
    bwa-mem2 index ${contig_name}_extracted.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>/dev/null)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${contig_name}_extracted.fa ${contig_name}_extracted.fa.amb
    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
