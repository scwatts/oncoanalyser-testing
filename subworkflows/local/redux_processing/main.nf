//
// Apply post-alignment processing
//


include { REDUX } from '../../../modules/local/redux/main'

workflow REDUX_PROCESSING {
    take:
    // Sample data
    ch_inputs           // channel: [mandatory] [ meta ]
    ch_dna_tumor        // channel: [mandatory] [ meta, [bam, ...], [bai, ...] ]
    ch_dna_normal       // channel: [mandatory] [ meta, [bam, ...], [bai, ...] ]
    ch_dna_donor        // channel: [mandatory] [ meta, [bam, ...], [bai, ...] ]

    // Reference data
    genome_fasta        // channel: [mandatory] /path/to/genome_fasta
    genome_ver          // channel: [mandatory] genome version
    genome_fai          // channel: [mandatory] /path/to/genome_fai
    genome_dict         // channel: [mandatory] /path/to/genome_dict
    unmap_regions       // channel: [mandatory] /path/to/unmap_regions
    msi_jitter_sites    // channel: [mandatory] /path/to/msi_jitter_sites

    // Params
    sequencing_platform // string:  [mandatory] sequencing platform
    umi_enable          // boolean: [mandatory] enable UMI processing
    umi_duplex_delim    // string:  [optional] UMI duplex delimiter
    targeted_mode       // boolean: [mandatory] Set targeted mode

    main:
    // Select and sort input sources, separating bytumor and normal
    // channel: runnable: [ meta, [bam, ...], [bai, ...] ]
    // channel: skip: [ meta ]
    ch_inputs_tumor_sorted = ch_dna_tumor
        .map { meta, bams, bais ->
            return [
                meta,
                Utils.hasExistingInput(meta, Constants.INPUT.BAM_DNA_TUMOR) ? [Utils.getInput(meta, Constants.INPUT.BAM_DNA_TUMOR)] : bams,
                Utils.hasExistingInput(meta, Constants.INPUT.BAI_DNA_TUMOR) ? [Utils.getInput(meta, Constants.INPUT.BAI_DNA_TUMOR)] : bais,
            ]
        }
        .branch { meta, bams, bais ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.BAM_REDUX_DNA_TUMOR)
            runnable: bams && !has_existing
            skip: true
                return meta
        }

    ch_inputs_normal_sorted = ch_dna_normal
        .map { meta, bams, bais ->
            return [
                meta,
                Utils.hasExistingInput(meta, Constants.INPUT.BAM_DNA_NORMAL) ? [Utils.getInput(meta, Constants.INPUT.BAM_DNA_NORMAL)] : bams,
                Utils.hasExistingInput(meta, Constants.INPUT.BAI_DNA_NORMAL) ? [Utils.getInput(meta, Constants.INPUT.BAI_DNA_NORMAL)] : bais,
            ]
        }
        .branch { meta, bams, bais ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.BAM_REDUX_DNA_NORMAL)
            runnable: bams && !has_existing
            skip: true
                return meta
        }

    ch_inputs_donor_sorted = ch_dna_donor
        .map { meta, bams, bais ->
            return [
                meta,
                Utils.hasExistingInput(meta, Constants.INPUT.BAM_DNA_DONOR) ? [Utils.getInput(meta, Constants.INPUT.BAM_DNA_DONOR)] : bams,
                Utils.hasExistingInput(meta, Constants.INPUT.BAI_DNA_DONOR) ? [Utils.getInput(meta, Constants.INPUT.BAI_DNA_DONOR)] : bais,
            ]
        }
        .branch { meta, bams, bais ->
            def has_existing = Utils.hasExistingInput(meta, Constants.INPUT.BAM_REDUX_DNA_DONOR)
            runnable: bams && !has_existing
            skip: true
            return meta
        }

    // Create process input channel
    // channel: [ meta_redux, [bam, ...], [bai, ...] ]
    ch_redux_inputs = channel.empty()
        .mix(
            ch_inputs_tumor_sorted.runnable.map { meta, bams, bais -> [meta, Utils.getTumorDnaSample(meta), 'tumor', bams, bais] },
            ch_inputs_normal_sorted.runnable.map { meta, bams, bais -> [meta, Utils.getNormalDnaSample(meta), 'normal', bams, bais] },
            ch_inputs_donor_sorted.runnable.map { meta, bams, bais -> [meta, Utils.getDonorDnaSample(meta), 'donor', bams, bais] },
        )
        .map { meta, meta_sample, sample_type, bams, bais ->

            def sample_id = meta_sample.getOrDefault('longitudinal_sample_id', meta_sample['sample_id'])

            def meta_redux = [
                key: meta.group_id,
                id: "${meta.group_id}_${sample_id}",
                sample_id: sample_id,
                sample_type: sample_type,
            ]

            return [meta_redux, bams, bais]
        }

    // Run process
    REDUX(
        ch_redux_inputs,
        genome_fasta,
        genome_ver,
        genome_fai,
        genome_dict,
        unmap_regions,
        msi_jitter_sites,
        sequencing_platform,
        umi_enable,
        umi_duplex_delim,
        targeted_mode,
    )

    // Combine TSV outputs into single channel for processing
    // channel: [ meta, bam, bai, bqr_tsv, jitter_tsv, ms_tsv, bqr_plot ]
    ch_redux_out = WorkflowOncoanalyser.groupByMeta(
        channel.topic('redux_bam'),
        channel.topic('redux_bqr_tsv'),
        channel.topic('redux_dup_freq_tsv'),
        channel.topic('redux_jitter_tsv'),
        channel.topic('redux_ms_tsv'),
        channel.topic('redux_bqr_plot'),
    )

    // Sort into a tumor and normal channel
    // channel: [ meta, bam, bai, bqr_tsv, dup_freq_tsv, jitter_tsv, ms_tsv, bqr_plot ]
    ch_redux_out_sorted = ch_redux_out
        .branch { meta, bam, bai, bqr_tsv, dup_freq_tsv, jitter_tsv, ms_tsv, bqr_plot ->
            assert ['tumor', 'normal', 'donor'].contains(meta.sample_type)
            tumor: meta.sample_type == 'tumor'
            normal: meta.sample_type == 'normal'
            donor: meta.sample_type == 'donor'
            placeholder: true
        }

    // Set outputs, restoring original meta, split into BAMs and TSVs
    // channel: [ meta, bam, bai, bqr_tsv, dup_freq_tsv, jitter_tsv, ms_tsv, bqr_plot ]
    ch_redux_tumor_out = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_redux_out_sorted.tumor, ch_inputs),
            ch_inputs_tumor_sorted.skip.map { meta -> [meta, [], [], [], [], [], [], []] },
        )
        .multiMap { meta, bam, bai, bqr_tsv, dup_freq_tsv, jitter_tsv, ms_tsv, bqr_plot ->
            bam: [meta, bam, bai]
            tsv: [meta, bqr_tsv, dup_freq_tsv, jitter_tsv, ms_tsv]
            plot: [meta, bqr_plot]
        }

    ch_redux_normal_out = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_redux_out_sorted.normal, ch_inputs),
            ch_inputs_normal_sorted.skip.map { meta -> [meta, [], [], [], [], [], [], []] },
        )
        .multiMap { meta, bam, bai, bqr_tsv, dup_freq_tsv, jitter_tsv, ms_tsv, bqr_plot ->
            bam: [meta, bam, bai]
            tsv: [meta, bqr_tsv, dup_freq_tsv, jitter_tsv, ms_tsv]
            plot: [meta, bqr_plot]
        }

    ch_redux_donor_out = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_redux_out_sorted.donor, ch_inputs),
            ch_inputs_donor_sorted.skip.map { meta -> [meta, [], [], [], [], [], [], []] },
        )
        .multiMap { meta, bam, bai, bqr_tsv, dup_freq_tsv, jitter_tsv, ms_tsv, bqr_plot ->
            bam: [meta, bam, bai]
            tsv: [meta, bqr_tsv, dup_freq_tsv, jitter_tsv, ms_tsv]
            plot: [meta, bqr_plot]
        }

    emit:
    dna_tumor_bam  = ch_redux_tumor_out.bam  // channel: [ meta, bam, bai ]
    dna_normal_bam = ch_redux_normal_out.bam // channel: [ meta, bam, bai ]
    dna_donor_bam  = ch_redux_donor_out.bam  // channel: [ meta, bam, bai ]

    dna_tumor_tsv  = ch_redux_tumor_out.tsv  // channel: [ meta, redux_tsv, ... ]
    dna_normal_tsv = ch_redux_normal_out.tsv // channel: [ meta, redux_tsv, ... ]
    dna_donor_tsv  = ch_redux_donor_out.tsv  // channel: [ meta, redux_tsv, ... ]

    dna_tumor_plot  = ch_redux_tumor_out.plot  // channel: [ meta, bqr_plot ]
    dna_normal_plot = ch_redux_normal_out.plot // channel: [ meta, bqr_plot ]
    dna_donor_plot  = ch_redux_donor_out.plot  // channel: [ meta, bqr_plot ]
}
