//
// WISP estimates tumor purity in longitudinal samples using WGS data of the primary
//


include { WISP } from '../../../modules/local/wisp/main'

workflow WISP_ANALYSIS {
    take:
    // Sample data
    ch_inputs                  // channel: [mandatory] [ meta ]
    ch_amber_out               // channel: [mandatory] [ meta, amber_dir ]
    ch_cobalt_out              // channel: [mandatory] [ meta, cobalt_dir ]
    ch_sage_somatic_append_out // channel: [mandatory] [ meta, sage_append_dir ]

    // Reference data
    genome_fasta               // channel: [mandatory] /path/to/genome_fasta
    genome_fai                 // channel: [mandatory] /path/to/genome_fai

    // Params
    purity_estimate_run_mode   //  string: [mandatory] Purity estimate run mode

    main:
    // Select input sources and sort
    // channel: runnable: [ meta, ... ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_amber_out,
        ch_cobalt_out,
        ch_sage_somatic_append_out,
    )
        .branch { meta, amber_dir, cobalt_dir, sage_append_dir ->

            def primary_purple_dir = Utils.getInput(meta, Constants.INPUT.PURPLE_DIR)
            def primary_amber_dir = Utils.getInput(meta, Constants.INPUT.AMBER_DIR)

            def runnable
            def purity_estimate_mode = purity_estimate_run_mode

            if (purity_estimate_run_mode == Constants.RunMode.WGTS) {
                runnable = primary_purple_dir && primary_amber_dir && sage_append_dir && amber_dir && cobalt_dir
            } else {
                runnable = primary_purple_dir && sage_append_dir
            }

            runnable: runnable
                return [meta, primary_purple_dir, primary_amber_dir, amber_dir, cobalt_dir, sage_append_dir]
            skip: true
                return meta
        }

    // Create process input channel
    // channel: [ meta_wisp, ... ]
    ch_wisp_inputs = ch_inputs_sorted.runnable

        .map { meta, primary_purple_dir, primary_amber_dir, amber_dir, cobalt_dir, sage_append_dir ->

            def meta_wisp = [
                key: meta.group_id,
                id: meta.group_id,
                subject_id: meta.subject_id,
                primary_id: Utils.getTumorDnaSampleName(meta, primary: true),
                longitudinal_id: Utils.getTumorDnaSampleName(meta, primary: false),
            ]

            return [meta_wisp, primary_purple_dir, primary_amber_dir, amber_dir, cobalt_dir, sage_append_dir]
        }


    // Run process
    WISP(
        ch_wisp_inputs,
        genome_fasta,
        genome_fai,
        purity_estimate_run_mode == Constants.RunMode.TARGETED,
    )
}
