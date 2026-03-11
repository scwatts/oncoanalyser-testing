
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { AMBER_PROFILING          } from '../subworkflows/local/amber_profiling'
include { BAMTOOLS_METRICS         } from '../subworkflows/local/bamtools_metrics'
include { CIDER_CALLING            } from '../subworkflows/local/cider_calling'
include { COBALT_PROFILING         } from '../subworkflows/local/cobalt_profiling'
include { ESVEE_CALLING            } from '../subworkflows/local/esvee_calling'
include { ISOFOX_QUANTIFICATION    } from '../subworkflows/local/isofox_quantification'
include { LILAC_CALLING            } from '../subworkflows/local/lilac_calling'
include { LINX_ANNOTATION          } from '../subworkflows/local/linx_annotation'
include { LINX_PLOTTING            } from '../subworkflows/local/linx_plotting'
include { ORANGE_REPORTING         } from '../subworkflows/local/orange_reporting'
include { PAVE_ANNOTATION          } from '../subworkflows/local/pave_annotation'
include { PEACH_CALLING            } from '../subworkflows/local/peach_calling'
include { PREPARE_REFERENCE        } from '../subworkflows/local/prepare_reference'
include { PREPARE_OUTPUTS_TARGETED } from '../subworkflows/local/prepare_outputs'
include { PURPLE_CALLING           } from '../subworkflows/local/purple_calling'
include { READ_ALIGNMENT_DNA       } from '../subworkflows/local/read_alignment_dna'
include { READ_ALIGNMENT_RNA       } from '../subworkflows/local/read_alignment_rna'
include { REDUX_PROCESSING         } from '../subworkflows/local/redux_processing'
include { SAGE_APPEND              } from '../subworkflows/local/sage_append'
include { SAGE_CALLING             } from '../subworkflows/local/sage_calling'
include { SAGE_PLOTTING            } from '../subworkflows/local/sage_plotting'

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TARGETED {
    take:
    inputs
    run_config
    params

    main:
    // Check input path parameters to see if they exist
    def checkPathParamList = [
        params.isofox_counts,
        params.isofox_gc_ratios,
        params.isofox_gene_ids,
        params.isofox_tpm_norm,
    ]

    checkPathParamList.each { param -> if (param) { file(param, checkIfExists: true) } }

    // Create input channel from parsed CSV
    // channel: [ meta ]
    ch_inputs = channel.fromList(inputs)

    // Set up reference data, assign more human readable variables
    prep_config = WorkflowMain.getPrepConfigFromSamplesheet(run_config)
    PREPARE_REFERENCE(
        prep_config,
        run_config,
        params,
    )
    ref_data = PREPARE_REFERENCE.out
    hmf_data = PREPARE_REFERENCE.out.hmf_data
    panel_data = PREPARE_REFERENCE.out.panel_data

    //
    // SUBWORKFLOW: Run read alignment to generate BAMs
    //
    // channel: [ meta, [bam, ...], [bai, ...] ]
    ch_align_dna_tumor_out = channel.empty()
    ch_align_dna_normal_out = channel.empty()
    ch_align_dna_donor_out = channel.empty()
    ch_align_rna_tumor_out = channel.empty()
    if (run_config.stages.alignment) {

        READ_ALIGNMENT_DNA(
            ch_inputs,
            ref_data.genome_fasta,
            ref_data.genome_bwamem2_index,
            params.max_fastq_records.toInteger(),
            params.fastp_umi_enabled,
            params.fastp_umi_location,
            params.fastp_umi_length,
            params.fastp_umi_skip,
        )

        READ_ALIGNMENT_RNA(
            ch_inputs,
            ref_data.genome_star_index,
        )

        ch_align_dna_tumor_out = ch_align_dna_tumor_out.mix(READ_ALIGNMENT_DNA.out.dna_tumor)
        ch_align_dna_normal_out = ch_align_dna_normal_out.mix(READ_ALIGNMENT_DNA.out.dna_normal)
        ch_align_dna_donor_out = ch_align_dna_donor_out.mix(READ_ALIGNMENT_DNA.out.dna_donor)
        ch_align_rna_tumor_out = ch_align_rna_tumor_out.mix(READ_ALIGNMENT_RNA.out.rna_tumor)

    } else {

        ch_align_dna_tumor_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_align_dna_normal_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_align_dna_donor_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_align_rna_tumor_out = ch_inputs.map { meta -> [meta, [], []] }

    }

    //
    // SUBWORKFLOW: Run REDUX for DNA BAMs
    //
    // channel: [ meta, bam, bai ]
    ch_redux_dna_tumor_bam_out = channel.empty()
    ch_redux_dna_normal_bam_out = channel.empty()
    ch_redux_dna_donor_bam_out = channel.empty()

    // channel: [ meta, bqr_tsv, jitter_tsv, ms_tsv ]
    ch_redux_dna_tumor_tsv_out = channel.empty()
    ch_redux_dna_normal_tsv_out = channel.empty()
    ch_redux_dna_donor_tsv_out = channel.empty()

    // channel: [ meta, bqr_plot ]
    ch_redux_dna_tumor_plot_out = channel.empty()
    ch_redux_dna_normal_plot_out = channel.empty()

    if (run_config.stages.redux) {

        REDUX_PROCESSING(
            ch_inputs,
            ch_align_dna_tumor_out,
            ch_align_dna_normal_out,
            ch_align_dna_donor_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            hmf_data.unmap_regions,
            hmf_data.msi_jitter_sites,
            params.sequencing_platform,
            params.redux_umi_enabled,
            params.redux_umi_duplex_delim,
            true,  // targeted_mode
        )

        ch_redux_dna_tumor_bam_out = ch_redux_dna_tumor_bam_out.mix(REDUX_PROCESSING.out.dna_tumor_bam)
        ch_redux_dna_normal_bam_out = ch_redux_dna_normal_bam_out.mix(REDUX_PROCESSING.out.dna_normal_bam)
        ch_redux_dna_donor_bam_out = ch_redux_dna_donor_bam_out.mix(REDUX_PROCESSING.out.dna_donor_bam)

        ch_redux_dna_tumor_tsv_out = ch_redux_dna_tumor_tsv_out.mix(REDUX_PROCESSING.out.dna_tumor_tsv)
        ch_redux_dna_normal_tsv_out = ch_redux_dna_normal_tsv_out.mix(REDUX_PROCESSING.out.dna_normal_tsv)
        ch_redux_dna_donor_tsv_out = ch_redux_dna_donor_tsv_out.mix(REDUX_PROCESSING.out.dna_donor_tsv)

        ch_redux_dna_tumor_plot_out = ch_redux_dna_tumor_plot_out.mix(REDUX_PROCESSING.out.dna_tumor_plot)
        ch_redux_dna_normal_plot_out = ch_redux_dna_normal_plot_out.mix(REDUX_PROCESSING.out.dna_normal_plot)

    } else {

        ch_redux_dna_tumor_bam_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_redux_dna_normal_bam_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_redux_dna_donor_bam_out = ch_inputs.map { meta -> [meta, [], []] }

        ch_redux_dna_tumor_tsv_out = ch_inputs.map { meta -> [meta, [], [], []] }
        ch_redux_dna_normal_tsv_out = ch_inputs.map { meta -> [meta, [], [], []] }
        ch_redux_dna_donor_tsv_out = ch_inputs.map { meta -> [meta, [], [], []] }

        ch_redux_dna_tumor_plot_out = ch_inputs.map { meta -> [meta, []] }
        ch_redux_dna_normal_plot_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // MODULE: Run Isofox to analyse RNA data
    //

    isofox_counts = params.isofox_counts ? file(params.isofox_counts) : panel_data.isofox_counts
    isofox_gc_ratios = params.isofox_gc_ratios ? file(params.isofox_gc_ratios) : panel_data.isofox_gc_ratios

    isofox_gene_ids = params.isofox_gene_ids ? file(params.isofox_gene_ids) : panel_data.isofox_gene_ids
    isofox_tpm_norm = params.isofox_tpm_norm ? file(params.isofox_tpm_norm) : panel_data.isofox_tpm_norm

    isofox_read_length = params.isofox_read_length != null ? params.isofox_read_length : Constants.DEFAULT_ISOFOX_READ_LENGTH_TARGETED

    // channel: [ meta, isofox_dir ]
    ch_isofox_out = channel.empty()
    if (run_config.stages.isofox) {

        ISOFOX_QUANTIFICATION(
            ch_inputs,
            ch_align_rna_tumor_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            hmf_data.ensembl_data_resources,
            hmf_data.known_fusion_data,
            hmf_data.isofox_gene_distribution,
            hmf_data.isofox_alt_sj_distribution,
            isofox_counts,
            isofox_gc_ratios,
            isofox_gene_ids,
            isofox_tpm_norm,
            params.isofox_functions,
            isofox_read_length,
        )

        ch_isofox_out = ch_isofox_out.mix(ISOFOX_QUANTIFICATION.out.isofox_dir)

    } else {

        ch_isofox_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Run AMBER to obtain b-allele frequencies
    //
    // channel: [ meta, amber_dir ]
    ch_amber_out = channel.empty()
    if (run_config.stages.amber) {

        AMBER_PROFILING(
            ch_inputs,
            ch_redux_dna_tumor_bam_out,
            ch_redux_dna_normal_bam_out,
            ch_redux_dna_donor_bam_out,
            ref_data.genome_version,
            hmf_data.heterozygous_sites,
            panel_data.target_region_bed,
            [],  // tumor_min_depth
        )

        ch_amber_out = ch_amber_out.mix(AMBER_PROFILING.out.amber_dir)

    } else {

        ch_amber_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Run COBALT to obtain read ratios
    //
    // channel: [ meta, cobalt_dir ]
    ch_cobalt_out = channel.empty()
    if (run_config.stages.cobalt) {

        COBALT_PROFILING(
            ch_inputs,
            ch_redux_dna_tumor_bam_out,
            ch_redux_dna_normal_bam_out,
            ref_data.genome_version,
            hmf_data.gc_profile,
            hmf_data.diploid_bed,
            panel_data.target_region_normalisation,
            true,  // targeted_mode
        )

        ch_cobalt_out = ch_cobalt_out.mix(COBALT_PROFILING.out.cobalt_dir)

    } else {

        ch_cobalt_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Call structural variants with ESVEE
    //
    // channel: [ meta, esvee_vcf ]
    ch_esvee_germline_out = channel.empty()
    ch_esvee_somatic_out = channel.empty()
    if (run_config.stages.esvee) {

        ESVEE_CALLING(
            ch_inputs,
            ch_redux_dna_tumor_bam_out,
            ch_redux_dna_normal_bam_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            ref_data.genome_img,
            hmf_data.known_fusions,
            hmf_data.esvee_pon_breakends,
            hmf_data.esvee_pon_breakpoints,
            hmf_data.decoy_sequences_image,
            hmf_data.repeatmasker_annotations,
            hmf_data.unmap_regions,
            panel_data.target_region_bed,
            params.sequencing_platform,
        )

        ch_esvee_germline_out = ch_esvee_germline_out.mix(ESVEE_CALLING.out.germline_vcf)
        ch_esvee_somatic_out = ch_esvee_somatic_out.mix(ESVEE_CALLING.out.somatic_vcf)

    } else {

        ch_esvee_germline_out = ch_inputs.map { meta -> [meta, []] }
        ch_esvee_somatic_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: call SNV, MNV, and small INDELS with SAGE
    //
    // channel: [ meta, sage_vcf, sage_tbi ]
    ch_sage_germline_vcf_out = channel.empty()
    ch_sage_somatic_vcf_out = channel.empty()
    // channel: [ meta, sage_dir ]
    ch_sage_germline_dir_out = channel.empty()
    ch_sage_somatic_dir_out = channel.empty()
    if (run_config.stages.sage) {

        SAGE_CALLING(
            ch_inputs,
            ch_redux_dna_tumor_bam_out,
            ch_redux_dna_normal_bam_out,
            ch_redux_dna_donor_bam_out,
            ch_redux_dna_tumor_tsv_out,
            ch_redux_dna_normal_tsv_out,
            ch_redux_dna_donor_tsv_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            hmf_data.sage_pon,
            hmf_data.sage_known_hotspots_somatic,
            hmf_data.sage_known_hotspots_germline,
            hmf_data.sage_highconf_regions,
            panel_data.driver_gene_panel,
            hmf_data.ensembl_data_resources,
            hmf_data.gnomad_resource,
            params.sequencing_platform,
            true,  // enable_germline
            true,  // targeted_mode
        )

        ch_sage_germline_vcf_out = ch_sage_germline_vcf_out.mix(SAGE_CALLING.out.germline_vcf)
        ch_sage_somatic_vcf_out = ch_sage_somatic_vcf_out.mix(SAGE_CALLING.out.somatic_vcf)
        ch_sage_germline_dir_out = ch_sage_germline_dir_out.mix(SAGE_CALLING.out.germline_dir)
        ch_sage_somatic_dir_out = ch_sage_somatic_dir_out.mix(SAGE_CALLING.out.somatic_dir)

    } else {

        ch_sage_germline_vcf_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_sage_somatic_vcf_out = ch_inputs.map { meta -> [meta, [], []] }
        ch_sage_germline_dir_out = ch_inputs.map { meta -> [meta, []] }
        ch_sage_somatic_dir_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Annotate variants with PAVE
    //
    // channel: [ meta, pave_vcf ]
    ch_pave_germline_out = channel.empty()
    ch_pave_somatic_out = channel.empty()
    if (run_config.stages.pave) {

        PAVE_ANNOTATION(
            ch_inputs,
            ch_sage_germline_vcf_out,
            ch_sage_somatic_vcf_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            panel_data.pon_artefacts,
            hmf_data.sage_pon,
            hmf_data.sage_blocklist_regions,
            hmf_data.sage_blocklist_sites,
            hmf_data.clinvar_annotations,
            hmf_data.segment_mappability,
            panel_data.driver_gene_panel,
            hmf_data.ensembl_data_resources,
            hmf_data.gnomad_resource,
            params.sequencing_platform,
        )

        ch_pave_germline_out = ch_pave_germline_out.mix(PAVE_ANNOTATION.out.germline)
        ch_pave_somatic_out = ch_pave_somatic_out.mix(PAVE_ANNOTATION.out.somatic)

    } else {

        ch_pave_germline_out = ch_inputs.map { meta -> [meta, []] }
        ch_pave_somatic_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Call CNVs, infer purity and ploidy, and recover low quality SVs with PURPLE
    //
    // channel: [ meta, purple_dir ]
    ch_purple_out = channel.empty()
    if (run_config.stages.purple) {

        PURPLE_CALLING(
            ch_inputs,
            ch_amber_out,
            ch_cobalt_out,
            ch_pave_somatic_out,
            ch_pave_germline_out,
            ch_esvee_somatic_out,
            ch_esvee_germline_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            hmf_data.gc_profile,
            hmf_data.sage_known_hotspots_somatic,
            hmf_data.sage_known_hotspots_germline,
            panel_data.driver_gene_panel,
            hmf_data.ensembl_data_resources,
            hmf_data.germline_amp_del_freq,
            panel_data.target_region_bed,
            panel_data.target_region_ratios,
            panel_data.target_region_msi_indels,
        )

        ch_purple_out = ch_purple_out.mix(PURPLE_CALLING.out.purple_dir)

    } else {

        ch_purple_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Append read data to SAGE VCF
    //
    // channel: [ meta, sage_append_vcf ]
    ch_sage_somatic_append_out = channel.empty()
    ch_sage_germline_append_out = channel.empty()
    if (run_config.stages.orange) {

        SAGE_APPEND(
            ch_inputs,
            ch_purple_out,
            ch_inputs.map { meta -> [meta, [], []] },      // ch_tumor_redux_bam
            ch_inputs.map { meta -> [meta, [], [], []] },  // ch_tumor_redux_tsv
            ch_align_rna_tumor_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            params.sequencing_platform,
            true,  // enable_germline
            true,  // targeted_mode
            false,  // purity_estimate_mode
        )

        ch_sage_somatic_append_out = ch_sage_somatic_append_out.mix(SAGE_APPEND.out.somatic_dir)
        ch_sage_germline_append_out = ch_sage_germline_append_out.mix(SAGE_APPEND.out.germline_dir)

    } else {

        ch_sage_somatic_append_out = ch_inputs.map { meta -> [meta, []] }
        ch_sage_germline_append_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Visualise SAGE variants
    //
    if (run_config.stages.sage_vis) {

        SAGE_PLOTTING(
            ch_inputs,
            ch_redux_dna_tumor_bam_out,
            ch_redux_dna_normal_bam_out,
            ch_redux_dna_donor_bam_out,
            ch_redux_dna_tumor_tsv_out,
            ch_redux_dna_normal_tsv_out,
            ch_redux_dna_donor_tsv_out,
            ch_purple_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            ref_data.genome_dict,
            hmf_data.sage_pon,
            hmf_data.sage_known_hotspots_somatic,
            hmf_data.sage_highconf_regions,
            hmf_data.ensembl_data_resources,
        )

    }

    //
    // SUBWORKFLOW: Group structural variants into higher order events with LINX
    //
    // channel: [ meta, linx_annotation_dir ]
    ch_linx_somatic_out = channel.empty()
    ch_linx_germline_out = channel.empty()
    if (run_config.stages.linx) {

        LINX_ANNOTATION(
            ch_inputs,
            ch_purple_out,
            ref_data.genome_version,
            hmf_data.ensembl_data_resources,
            hmf_data.known_fusion_data,
            panel_data.driver_gene_panel,
        )

        ch_linx_somatic_out = ch_linx_somatic_out.mix(LINX_ANNOTATION.out.somatic)
        ch_linx_germline_out = ch_linx_germline_out.mix(LINX_ANNOTATION.out.germline)

    } else {

        ch_linx_somatic_out = ch_inputs.map { meta -> [meta, []] }
        ch_linx_germline_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Visualise LINX annotations
    //
    // channel: [ meta, linx_visualiser_dir ]
    ch_linx_somatic_visualiser_dir_out = channel.empty()
    if (run_config.stages.linx) {

        LINX_PLOTTING(
            ch_inputs,
            ch_linx_somatic_out,
            ch_amber_out,
            ch_cobalt_out,
            ch_purple_out,
            ref_data.genome_version,
            hmf_data.ensembl_data_resources,
        )

        ch_linx_somatic_visualiser_dir_out = ch_linx_somatic_visualiser_dir_out.mix(LINX_PLOTTING.out.visualiser_dir)

    } else {

        ch_linx_somatic_visualiser_dir_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Run Bam Tools to generate stats required for downstream processes
    //
    // channel: [ meta, metrics ]
    ch_bamtools_somatic_out = channel.empty()
    ch_bamtools_germline_out = channel.empty()
    if (run_config.stages.bamtools) {

        BAMTOOLS_METRICS(
            ch_inputs,
            ch_redux_dna_tumor_bam_out,
            ch_redux_dna_normal_bam_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            panel_data.driver_gene_panel,
            hmf_data.ensembl_data_resources,
            panel_data.target_region_bed,
        )

        ch_bamtools_somatic_out = ch_bamtools_somatic_out.mix(BAMTOOLS_METRICS.out.somatic)
        ch_bamtools_germline_out = ch_bamtools_germline_out.mix(BAMTOOLS_METRICS.out.germline)

    } else {

        ch_bamtools_somatic_out = ch_inputs.map { meta -> [meta, []] }
        ch_bamtools_germline_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Run CIDER to identify and annotate CDR3 sequences of IG and TCR loci
    //
    if (run_config.stages.cider) {

        CIDER_CALLING(
            ch_inputs,
            ch_redux_dna_tumor_bam_out,
            ch_align_rna_tumor_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_dict,
            ref_data.genome_img,
        )

    }

    //
    // SUBWORKFLOW: Run LILAC for HLA typing and somatic CNV and SNV calling
    //
    // channel: [ meta, lilac_dir ]
    ch_lilac_out = channel.empty()
    if (run_config.stages.lilac) {

        LILAC_CALLING(
            ch_inputs,
            ch_redux_dna_tumor_bam_out,
            ch_redux_dna_normal_bam_out,
            ch_align_rna_tumor_out,
            ch_purple_out,
            ref_data.genome_fasta,
            ref_data.genome_version,
            ref_data.genome_fai,
            hmf_data.lilac_resources,
            true,  // targeted_mode,
            params.sequencing_platform,
        )

        ch_lilac_out = ch_lilac_out.mix(LILAC_CALLING.out.lilac_dir)

    } else {

        ch_lilac_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Run PEACH to call germline haplotypes and report pharmacogenomics
    //
    // channel: [ meta, peach_dir ]
    ch_peach_out = channel.empty()
    if (run_config.stages.peach) {

        PEACH_CALLING(
            ch_inputs,
            ch_purple_out,
            hmf_data.peach_haplotypes,
            hmf_data.peach_haplotype_functions,
            hmf_data.peach_drug_info,
        )

        ch_peach_out = ch_peach_out.mix(PEACH_CALLING.out.peach_dir)

    } else {

        ch_peach_out = ch_inputs.map { meta -> [meta, []] }

    }

    //
    // SUBWORKFLOW: Run ORANGE to generate static PDF report
    //
    if (run_config.stages.orange) {

        // Create placeholder channels for empty remaining channels
        ch_chord_out = ch_inputs.map { meta -> [meta, []] }
        ch_cuppa_out = ch_inputs.map { meta -> [meta, []] }
        ch_sigs_out = ch_inputs.map { meta -> [meta, []] }
        ch_virusinterpreter_out = ch_inputs.map { meta -> [meta, []] }

        ORANGE_REPORTING(
            ch_inputs,
            ch_redux_dna_tumor_plot_out,
            ch_redux_dna_normal_plot_out,
            ch_bamtools_somatic_out,
            ch_bamtools_germline_out,
            ch_sage_somatic_dir_out,
            ch_sage_germline_dir_out,
            ch_sage_somatic_append_out,
            ch_sage_germline_append_out,
            ch_purple_out,
            ch_linx_somatic_out,
            ch_linx_somatic_visualiser_dir_out,
            ch_linx_germline_out,
            ch_virusinterpreter_out,
            ch_chord_out,
            ch_sigs_out,
            ch_lilac_out,
            ch_cuppa_out,
            ch_peach_out,
            ch_isofox_out,
            ref_data.genome_version,
            hmf_data.disease_ontology,
            hmf_data.cohort_mapping,
            panel_data.driver_gene_panel,
            hmf_data.sigs_etiology,
            true,  // targeted_mode
        )

    }

    //
    // SUBWORKFLOW: Prepare results for publishing
    //
    PREPARE_OUTPUTS_TARGETED()

    //
    // TASK: Aggregate software versions
    //
    def topic_versions = channel.topic('versions')
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(topic_versions.versions_file)
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'software_versions.yml',
            sort: true,
            newLine: true,
        )

    emit:
    results = PREPARE_OUTPUTS_TARGETED.out.results
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
