/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_dbcan_pipeline'

// new modules added by Xinpeng

include { KRAKEN2_DB_PREPARATION          } from '../modules/local/kraken2_db/'

include { KRAKEN2_KRAKEN2                 } from '../modules/nf-core/kraken2/kraken2/main'
include { KRAKENTOOLS_EXTRACTKRAKENREADS  } from '../modules/nf-core/krakentools/extractkrakenreads/main'
include { MEGAHIT                         } from '../modules/nf-core/megahit/main'
include { PRODIGAL                        } from '../modules/nf-core/prodigal/main'


// new subworkflows added by Xinpeng

include { FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS } from '../subworkflows/nf-core/fastq_extract_kraken_krakentools'
include { FASTQC_TRIMGALORE                } from '../subworkflows/local/fastqc_trimgalore'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Prepare the project parameters and databases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if(params.kraken_db_archive){
    ch_kraken2_db_file = file(params.kraken_db_archive, checkIfExists: true)
} else {
    ch_kraken2_db_file = []
}

if(params.kraken_tax){
    ch_kraken2_tax = params.kraken_tax
} else {
    ch_kraken2_tax = []
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DBCAN {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run FastQC
    //

    FASTQC_TRIMGALORE (
        ch_samplesheet,
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.fastqc_zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    //
    // MODULE: Kraken2 Build Database
    //
if (!ch_kraken2_db_file.isEmpty()) {
    if (ch_kraken2_db_file.extension in ['gz', 'tgz']) {
        // Expects to be tar.gz!
        ch_db_for_kraken2 = KRAKEN2_DB_PREPARATION(ch_kraken2_db_file).db
    } else if (ch_kraken2_db_file.isDirectory()) {
        // Directly used as database path
        ch_db_for_kraken2 = Channel.fromPath(ch_kraken2_db_file)
    } else {
        ch_db_for_kraken2 = Channel.empty()
    }
} else {
    ch_db_for_kraken2 = Channel.empty()
}

    //
    // Subworkflow: Extract classified Kraken2 reads by taxonomic id
    //
    FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS (
            FASTQC_TRIMGALORE.out.reads,
            ch_db,
            params.tax_id)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS.out.report)
    ch_versions = ch_versions.mix ( FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS.out.versions )

    // MODULE: Megahit to assemble metagenomics
    MEGAHIT ( FASTQ_EXTRACT_KRAKEN_KRAKENTOOLS.out.extracted_kraken2_reads )
    ch_versions = ch_versions.mix ( MEGAHIT.out.versions )

    //
    // MODULE: Prodigal to find genes in bacteria and archaea
    //
    PRODIGAL ( ch_contigs, 'gff' )
    ch_versions = ch_versions.mix(PRODIGAL.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
