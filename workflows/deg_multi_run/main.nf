#!/usr/bin/env nextflow


annotated_gem_ch = Channel.fromPath( params.annotated_gem_path )
formula_design_ch = Channel.value( params.formula_design )
formula_variables_ch = Channel.value( params.formula_variables )



process Create_DGE_Input {

    input:
        file ( gem_file ) from annotated_gem_ch
        val  ( formula_variables ) from formula_variables_ch

    output:
        file ( "complete_counts.csv" ) into complete_gem_csv_ch
        file ( "dropped_counts.csv" ) into dropped_gem_csv_ch
        file ( "annotations.csv" ) into sample_annotations_ch

    // beforeScript 'conda activate py38'

    script:
    """
    #!/usr/bin/env python

    import GSForge as gsf

    agem = gsf.AnnotatedGEM("${gem_file}")

    dropped_counts, labels = gsf.get_data(
      agem,
      count_mask="dropped",
      annotation_variables=${formula_variables}
    )

    complete_counts, _ = gsf.get_data(
      agem,
      count_mask="complete",
    )

    dropped_counts = gsf.utils.R_interface.Py_counts_to_R(dropped_counts)
    dropped_counts = dropped_counts.round()

    complete_counts = gsf.utils.R_interface.Py_counts_to_R(complete_counts)
    complete_counts = complete_counts.round()

    labels = labels.to_dataframe()

    complete_counts.to_csv("complete_counts.csv", index_label=False)
    dropped_counts.to_csv("dropped_counts.csv", index_label=False)
    labels.to_csv("annotations.csv")

    """
}



process Run_DESeq2 {

  container = 'genomicpariscentre/deseq2:latest'

    input:
        file ( gem_csv ) from dropped_gem_csv_ch
        file ( sample_annotations ) from sample_annotations_ch
        val  ( formula_design ) from formula_design_ch

    output:
        file( "*.csv" ) into results_DESeq2_ch

    script:
    """
    #!/usr/bin/env Rscript

    library("DESeq2")

    dropped_counts <- read.csv(file = "${gem_csv}")
    labels <- read.csv(file = "${sample_annotations}")

    dds <- DESeqDataSetFromMatrix(countData = dropped_counts,
                                  colData = labels,
                                  design = ${formula_design})
    dds <- DESeq(dds)

    full_deseq_df <- data.frame(results(dds))
    write.csv(full_deseq_df, "deseq2_results.csv")
    """
}
