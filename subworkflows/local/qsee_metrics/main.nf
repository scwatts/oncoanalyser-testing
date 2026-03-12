//
// Qsee calculates and visualises QC metrics
//

include { QSEE } from '../../../modules/local/qsee/main'

workflow QSEE_METRICS {
    take:
    // Sample data
    ch_inputs                // channel: [mandatory] [ meta ]
    ch_redux_tsvs_tumor      // channel: [mandatory] [ meta, redux_tsv, ... ]
    ch_redux_tsvs_normal     // channel: [mandatory] [ meta, redux_tsv, ... ]
    ch_bamtools_tumor        // channel: [mandatory] [ meta, metrics_dir ]
    ch_bamtools_normal       // channel: [optional]  [ meta, metrics_dir ]
    ch_cobalt                // channel: [optional]  [ meta, cobalt_dir ]
    ch_esvee                 // channel: [optional]  [ meta, esvee_dir ]
    ch_purple                // channel: [mandatory] [ meta, purple_dir ]

    // Reference data
    driver_gene_panel        // channel: [mandatory] /path/to/driver_gene_panel
    qsee_cohort_percentiles  // channel: [mandatory] /path/to/cohort_percentiles

    // Params
    targeted_mode            // boolean: [mandatory] Set targeted mode

    main:
    // Select and route inputs
    // channel: { meta, redux_tsvs_tumor, redux_tsvs_normal, bamtools_tumor_dir, bamtools_normal_dir, cobalt_dir, esvee_dir, purple_dir }
    ch_inputs_sorted = WorkflowOncoanalyser.groupByMeta(
        ch_redux_tsvs_tumor, ch_redux_tsvs_normal,
        ch_bamtools_tumor, ch_bamtools_normal,
        ch_cobalt,
        ch_esvee,
        ch_purple,
    )
        .map { meta,
            tumor_bqr_tsv , tumor_dup_freq_tsv , tumor_jitter_tsv , tumor_ms_tsv,
            normal_bqr_tsv, normal_dup_freq_tsv, normal_jitter_tsv, normal_ms_tsv,
            bamtools_tumor_dir, bamtools_normal_dir,
            cobalt_dir,
            esvee_dir,
            purple_dir ->

            def tumor_redux_tsvs = [
                tumor_bqr_tsv ?: Utils.getInput(meta, Constants.INPUT.REDUX_BQR_TSV_TUMOR),
                tumor_jitter_tsv ?: Utils.getInput(meta, Constants.INPUT.REDUX_JITTER_TSV_TUMOR),
                tumor_ms_tsv ?: Utils.getInput(meta, Constants.INPUT.REDUX_MS_TSV_TUMOR),
            ]

            def normal_redux_tsvs = [
                normal_bqr_tsv ?: Utils.getInput(meta, Constants.INPUT.REDUX_BQR_TSV_NORMAL),
                normal_jitter_tsv ?: Utils.getInput(meta, Constants.INPUT.REDUX_JITTER_TSV_NORMAL),
                normal_ms_tsv ?: Utils.getInput(meta, Constants.INPUT.REDUX_MS_TSV_NORMAL),
            ]

            tumor_redux_tsvs = tumor_redux_tsvs.findAll { it -> it != [] }
            normal_redux_tsvs = normal_redux_tsvs.findAll { it -> it != [] }

            return [
                meta,
                tumor_redux_tsvs,
                normal_redux_tsvs,
                Utils.selectCurrentOrExisting(bamtools_tumor_dir, meta, Constants.INPUT.BAMTOOLS_DIR_TUMOR),
                Utils.selectCurrentOrExisting(bamtools_normal_dir, meta, Constants.INPUT.BAMTOOLS_DIR_NORMAL),
                Utils.selectCurrentOrExisting(cobalt_dir, meta, Constants.INPUT.COBALT_DIR),
                Utils.selectCurrentOrExisting(esvee_dir, meta, Constants.INPUT.ESVEE_DIR),
                Utils.selectCurrentOrExisting(purple_dir, meta, Constants.INPUT.PURPLE_DIR),
            ]
        }
        .branch { meta, tumor_redux_tsvs, normal_redux_tsvs, bamtools_tumor_dir, bamtools_normal_dir, cobalt_dir, esvee_dir, purple_dir ->
            runnable: bamtools_tumor_dir && purple_dir
            skip: true
                return meta
        }

    // Create process input channel; form metadata
    // channel: [ qsee_meta, redux_tsvs_tumor, redux_tsvs_normal, bamtools_tumor_dir, bamtools_normal_dir, cobalt_dir, esvee_dir, purple_dir ]
    ch_qsee_inputs = ch_inputs_sorted.runnable
        .map { meta, tumor_redux_tsvs, normal_redux_tsvs, bamtools_tumor_dir, bamtools_normal_dir, cobalt_dir, esvee_dir, purple_dir ->

            def meta_qsee = [
                key: meta.group_id,
                id: meta.group_id,
                tumor_id: Utils.getTumorDnaSampleName(meta),
            ]

            if (normal_redux_tsvs || bamtools_normal_dir) {
                meta_qsee.normal_id = Utils.getNormalDnaSampleName(meta)
            }

            return [meta_qsee, tumor_redux_tsvs, normal_redux_tsvs, bamtools_tumor_dir, bamtools_normal_dir, cobalt_dir, esvee_dir, purple_dir]
        }

    // Run process
    QSEE(
        ch_qsee_inputs,
        driver_gene_panel,
        qsee_cohort_percentiles,
        targeted_mode,
    )

    // Set outputs, restoring original meta
    // channel: [ meta, qsee_dir ]
    ch_outputs = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(channel.topic('qsee_dir'), ch_inputs),
            ch_inputs_sorted.skip.map { meta -> [meta, []] },
        )

    emit:
    qsee_dir = ch_outputs  // channel: [ meta, qsee_dir ]
}
