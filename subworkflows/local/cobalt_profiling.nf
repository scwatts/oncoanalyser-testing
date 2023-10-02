//
// COBALT calculates read ratios between tumor and normal samples
//

import Utils

include { COBALT } from '../../modules/local/cobalt/main'

workflow COBALT_PROFILING {
    take:
        // Sample data
        ch_inputs                   // channel: [mandatory] [ meta ]

        // Reference data
        gc_profile                  // channel: [mandatory] /path/to/gc_profile
        diploid_bed                 // channel: [optional]  /path/to/diploid_bed
        target_region_normalisation // channel: [optional]  /path/to/target_region_normalisation

        // Params
        run_config                  // channel: [mandatory] run configuration

    main:
        // Channel for version.yml files
        // channel: [ versions.yml ]
        ch_versions = Channel.empty()

        // Select input sources
        // channel: [ meta_cobalt, tumor_bam, normal_bam, tumor_bai, normal_bai ]
        ch_cobalt_inputs = ch_inputs
            .map { meta ->
                def meta_cobalt = [
                    key: meta.id,
                    id: meta.id,
                    tumor_id: Utils.getTumorDnaSampleName(meta),
                ]
                def tumor_bam = Utils.getTumorDnaBam(meta)

                def normal_bam = []
                def normal_bai = []

                if (run_config.type == Constants.RunType.TUMOR_NORMAL) {

                    assert [Constants.RunMode.DNA, Constants.RunMode.DNA_RNA].contains(run_config.mode)

                    meta_cobalt.normal_id = Utils.getNormalDnaSampleName(meta)
                    normal_bam = Utils.getNormalDnaBam(meta)
                    normal_bai = "${normal_bam}.bai"

                }

                return [meta_cobalt, tumor_bam, normal_bam, "${tumor_bam}.bai", normal_bai]
            }

        // Set reference files on run type
        if (run_config.type != Constants.RunType.TUMOR_ONLY) {
            diploid_bed = []
        }

        // Run process
        COBALT(
            ch_cobalt_inputs,
            gc_profile,
            diploid_bed,
            target_region_normalisation,
        )

        // Set outputs, restoring original meta
        ch_outputs = WorkflowOncoanalyser.restoreMeta(COBALT.out.cobalt_dir, ch_inputs)
        ch_versions = ch_versions.mix(COBALT.out.versions)

    emit:
        cobalt_dir = ch_outputs  // channel: [ meta, cobalt_dir ]

        versions   = ch_versions // channel: [ versions.yml ]
}
