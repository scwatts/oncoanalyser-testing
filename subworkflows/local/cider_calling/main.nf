//
// CIDER identifies and annotates CDR3 sequences of IG and TCR loci
//


include { CIDER } from '../../../modules/local/cider/main'

workflow CIDER_CALLING {
    take:
    // Sample data
    ch_inputs        // channel: [mandatory] [ meta ]
    ch_tumor_dna_bam // channel: [mandatory] [ meta, bam, bai ]
    ch_tumor_rna_bam // channel: [mandatory] [ meta, bam, bai ]

    // Reference data
    genome_fasta     // channel: [mandatory] /path/to/genome_fasta
    genome_version   // channel: [mandatory] genome version
    genome_dict      // channel: [mandatory] /path/to/genome_dict
    genome_img       // channel: [optional]  /path/to/genome_img

    main:
    // Sort inputs, separate by DNA and RNA
    // channel: runnable: [ meta, bam, bai ]
    // channel: skip: [ meta ]
    ch_inputs_tumor_dna_sorted = ch_tumor_dna_bam
        .map { meta, bam, bai ->
            return [
                meta,
                Utils.selectCurrentOrExisting(bam, meta, Constants.INPUT.BAM_REDUX_DNA_TUMOR),
                bai ?: Utils.getInput(meta, Constants.INPUT.BAI_DNA_TUMOR),
            ]
        }
        .branch { meta, bam, bai ->
            runnable: bam
            skip: true
                return meta
        }

    // channel: runnable: [ meta, bam, bai ]
    // channel: skip: [ meta ]
    ch_inputs_tumor_rna_sorted = ch_tumor_rna_bam
        .map { meta, bam, bai ->
            return [
                meta,
                Utils.selectCurrentOrExisting(bam, meta, Constants.INPUT.BAM_RNA_TUMOR),
                bai ?: Utils.getInput(meta, Constants.INPUT.BAI_RNA_TUMOR),
            ]
        }
        .branch { meta, bam, bai ->
            runnable: bam
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_cider, bam, bai ]
    ch_cider_inputs = channel.empty()
        .mix(
            ch_inputs_tumor_dna_sorted.runnable.map { meta, bam, bai -> [meta, Utils.getTumorDnaSample(meta), bam, bai] },
            ch_inputs_tumor_rna_sorted.runnable.map { meta, bam, bai -> [meta, Utils.getTumorRnaSample(meta), bam, bai] },
        )
        .map { meta, meta_sample, bam, bai ->

            def meta_cider = [
                key: meta.group_id,
                id: "${meta.group_id}_${meta_sample.sample_id}",
                sample_id: meta_sample.sample_id,
            ]

            return [meta_cider, bam, bai]
        }

    // Run process
    CIDER(
        ch_cider_inputs,
        genome_fasta,
        genome_version,
        genome_dict,
        genome_img,
    )
}
