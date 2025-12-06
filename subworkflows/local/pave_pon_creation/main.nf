//
// PAVE PON creation prepares the panel-specific small variant artefact resource
//


include { PAVE_PON_PANEL_CREATION } from '../../../modules/local/pave/pon_creation/main'


workflow PAVE_PON_CREATION {
    take:
    // Sample data
    ch_sage_somatic_vcf // channel: [mandatory] [ meta, sage_somatic_vcf, sage_somatic_tbi ]

    // Reference data
    genome_version      // channel: [mandatory] genome version

    main:
    // Create process input channel
    // channel: [ [sage_vcf, ...], [sage_tbi, ...] ]
    ch_pave_inputs = ch_sage_somatic_vcf
        .map { meta, sage_vcf, sage_tbi ->
            return [
                Utils.selectCurrentOrExisting(sage_vcf, meta, Constants.INPUT.SAGE_VCF_TUMOR),
                Utils.selectCurrentOrExisting(sage_tbi, meta, Constants.INPUT.SAGE_VCF_TBI_TUMOR),
            ]
        }
        .collect(flat: false)
        .map { d -> d.transpose() }

    // Run process
    PAVE_PON_PANEL_CREATION(
        ch_pave_inputs,
        genome_version,
    )
}
