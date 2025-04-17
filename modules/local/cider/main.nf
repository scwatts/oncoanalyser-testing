process CIDER {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/scwatts/cider:1.0.3--0'

    input:
    tuple val(meta), path(bam), path(bai)
    val genome_ver
    file human_blastdb

    output:
    tuple val(meta), path('cider/*'), emit: cider_dir
    path 'versions.yml'             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    java \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        -cp /opt/cider/cider.jar com.hartwig.hmftools.cider.CiderApplication \\
            ${args} \\
            -sample ${meta.sample_id} \\
            -bam ${bam} \\
            -blast \$(which blastn | sed 's#/bin/blastn##') \\
            -blast_db ${human_blastdb} \\
            -ref_genome_version ${genome_ver} \\
            -threads ${task.cpus} \\
            -write_cider_bam \\
            -output_dir cider/

    # NOTE(SW): hard coded since there is no reliable way to obtain version information.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cider: 1.0.3
    END_VERSIONS
    """

    stub:
    """
    mkdir -p cider/

    touch cider/${meta.sample_id}.cider.bam
    touch cider/${meta.sample_id}.cider.blastn_match.tsv.gz
    touch cider/${meta.sample_id}.cider.layout.gz
    touch cider/${meta.sample_id}.cider.locus_stats.tsv
    touch cider/${meta.sample_id}.cider.vdj.tsv.gz

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
