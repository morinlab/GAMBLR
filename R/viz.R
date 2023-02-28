
# global set of aliases for finding specific sets of colours
colour_aliases = list("COO_consensus" = "coo", "COO" = "coo", "DHITsig_consensus" = "coo",
                      "pathology" = "pathology", "analysis_cohort" = "pathology", "group" = "pathology",
                      "FL_group" = "FL", "lymphgen" = "lymphgen", "lymphgen_with_cnv" = "lymphgen",
                      "bcl2_ba" = "pos_neg", "BCL2_status" = "pos_neg", "myc_ba" = "pos_neg",
                      "bcl6_ba" = "pos_neg", "EBV_status_inf"="EBV_status",
                      "manta_BCL2_sv" = "pos_neg", "manual_BCL2_sv" = "pos_neg", "manta_MYC_sv" = "pos_neg")


#' Plot a rainfall plot for one sample. This function takes in MAF data frame, or path to custom MAF file.
#' If non are specified, the SSM will be obtained though GAMBLR directly.
#'
#' @param this_sample_id Sample id for the sample to display. This is argument is not required if you want a multi-sample plot but is otherwise needed.
#' @param label_ashm_genes Boolean argument indicating whether the aSHM regions will be labeled or not.
#' @param projection Specify projection (grch37 or hg38) of mutations. Default is grch37.
#' @param chromosome Provide one or more chromosomes to plot. The chr prefix can be inconsistent with projection and will be handled.
#' @param this_maf Specify custom MAF data frame of mutations.
#' @param maf_path Specify path to MAF file if it is not already loaded into data frame.
#' @param zoom_in_region Provide a specific region in the format "chromosome:start-end" to zoom in to a specific region.
#' @param label_sv Boolean argument to specify whether label SVs or not. Only supported if a specific chromosome or zoom in region are specified.
#' @param seq_type Specify one of "genome" or "capture" when relying on the function to obtain mutations from a region (i.e. if you haven't provided a MAF or single sample_id)
#'
#' @return a ggplot2 plot. Print it using print() or save it using ggsave()
#' @export
#' @import ggplot2 dplyr ggrepel
#'
#' @examples
#' prettyRainfallPlot("Raji")
#' prettyRainfallPlot("Raji", chromosome = c(3,9,"chr14",22,"X"))
#' prettyRainfallPlot("Raji", chromosome = c(3,9), projection = "hg38", label_ashm_genes = FALSE)
#' prettyRainfallPlot("Raji", zoom_in_region = "8:125252796-135253201", label_sv = TRUE)
#' prettyRainfallPlot("Raji", chromosome = 6, label_sv = TRUE)
#' prettyRainfallPlot( zoom_in_region = "chr3:5,221,286-5,269,723", seq_type="genome") #multi-sample rainfall plot for one gene region

prettyRainfallPlot = function(this_sample_id,
                              label_ashm_genes = TRUE,
                              projection = "grch37",
                              chromosome,
                              this_maf,
                              maf_path,
                              zoom_in_region,
                              seq_type,
                              label_sv = FALSE) {
  if (missing(this_sample_id)) {
    warning("No sample_id was provided. Using all mutations in the MAF within your region!")
    if(missing(zoom_in_region)){
      stop("Must provide a zoom_in_region to plot when showing data from more than one patient")
    }
  }

  # allow user to specify chromosome prefix inconsistent with chromosome names
  if (!missing(chromosome)) {
    chromosome = standardize_chr_prefix(incoming_vector = chromosome, projection = projection)
  }

  # allow to zoom in to a specific region
  if (!missing(zoom_in_region)) {
    region = zoom_in_region
    zoom_in_region = region_to_chunks(zoom_in_region)
    zoom_in_region$chromosome = standardize_chr_prefix(incoming_vector = zoom_in_region$chromosome,
                                                       projection = projection)
    zoom_in_region$start = as.numeric(zoom_in_region$start)
    zoom_in_region$end = as.numeric(zoom_in_region$end)
  }

  if (label_ashm_genes) {
    if (projection == "grch37") {
      ashm_regions = grch37_ashm_regions %>%
        dplyr::rename("start" = "hg19_start",
                      "end" = "hg19_end",
                      "Chromosome" = "chr_name") %>%
        dplyr::mutate(Chromosome = str_remove(Chromosome, pattern = "chr"))
    } else if (projection == "hg38") {
      ashm_regions = hg38_ashm_regions %>%
        rename("start" = "hg38_start",
               "end" = "hg38_end",
               "Chromosome" = "chr_name")
    } else {
      stop("Please specify one of grch37 or hg38 projections")
    }
    if (!missing(chromosome)) {
      ashm_regions = dplyr::filter(ashm_regions, Chromosome %in% chromosome)
    }
    if (!missing(zoom_in_region)) {
      ashm_regions = dplyr::filter(
        ashm_regions,
        (
          Chromosome %in% zoom_in_region$chromosome &
            start >= zoom_in_region$start &
            end <= zoom_in_region$end
        )
      )
    }
    ashm_regions = ashm_regions %>%
      group_by(gene) %>%
      slice_head() %>%
      ungroup()

    # this will be needed for consistent labeling with rainfall plots
    ashm_regions = ashm_regions %>%
      arrange(match(
        Chromosome,
        str_sort(ashm_regions$Chromosome, numeric = TRUE)
      ))
    ashm_regions = ashm_regions %>%
      mutate(Chromosome_f = factor(Chromosome, levels = unique(ashm_regions$Chromosome)))
  }

  # if user is subsetting by chromosome or zooming in to a specific region, it is possible there are no aSHM features to show
  # handle this case separately
  if (nrow(ashm_regions) == 0) {
    message(
      "Warning: after subsetting to a regions you requested to plot, there are no aSHM features to overlap on the final graph."
    )
    label_ashm_genes = FALSE
  }

  # get ssm for the requested sample
  if (!missing(this_maf)) {
    if(missing(this_sample_id)){
      these_ssm=this_maf
      this_sample_id = "all samples"
    }else{
    message ("Using the suppplied MAF df to obrain ser of SSM for the specified sample ...")
    these_ssm = this_maf %>%
      dplyr::filter(Tumor_Sample_Barcode %in% this_sample_id)
    }
  } else if (!missing (maf_path)) {
    message ("Path to custom MAF file was provided, reading SSM using the custom path ...")

    this_maf = suppressMessages(read_tsv(maf_path))
    if(!missing(this_sample_id)){
      this_maf = this_maf %>% dplyr::filter(Tumor_Sample_Barcode %in% this_sample_id)
    }else{
      this_sample_id = "all samples"
    }
  } else if(!missing(this_sample_id)) {
    message ("MAF df or path to custom MAF file was not provided, getting SSM using GAMBLR ...")
    these_ssm = get_ssm_by_sample(this_sample_id,
                                  projection = projection, this_seq_type = seq_type)
  }else if(!missing(seq_type)){
    if(missing(this_sample_id)){
      this_sample_id = "all samples"
    }
    message(paste("Will use all mutations for",seq_type, "in this region:",zoom_in_region))
    these_ssm = get_ssm_by_region(region = region,seq_type = seq_type,projection=projection,
                                  maf_columns = c("Hugo_Symbol","Chromosome", "Start_Position", "End_Position", "Tumor_Sample_Barcode", "t_alt_count","Reference_Allele","Tumor_Seq_Allele2"),
                                  maf_column_types = c("c","c","i","i","c","i","c","c"))
  }

  # do rainfall calculation using lag
  rainfall_points = dplyr::select(
    these_ssm,
    Hugo_Symbol,
    Chromosome,
    Start_Position,
    End_Position,
    Reference_Allele,
    Tumor_Seq_Allele2
  ) %>%
    arrange(Chromosome, Start_Position) %>%
    group_by(Chromosome) %>%  # group by chromosome to calculate lag per chromosome
    dplyr::mutate(
      IMD = Start_Position - dplyr::lag(Start_Position),
      # used for coloring
      # all indels are squished to the same color
      Substitution = ifelse((
        Reference_Allele %in% c("A", "T",  "C", "G") &
          Tumor_Seq_Allele2 %in% c("A", "T",  "C", "G")
      ),
      paste(Reference_Allele, Tumor_Seq_Allele2, sep = '>'),
      "InDel"
      )
    ) %>%
    dplyr::mutate(IMD = log(IMD)) %>%
    ungroup() %>%
    drop_na(IMD) # for the first point of each chromosome, NAs are produced generating a warning message

  # collapse substitutions into classes
  rainfall_points$Substitution = rainfall_conv[as.character(rainfall_points$Substitution)]

  # ensure order of grids in the plot is sorted
  rainfall_points = rainfall_points %>%
    arrange(match(
      Chromosome,
      str_sort(rainfall_points$Chromosome, numeric = TRUE)
    ))
  rainfall_points = rainfall_points %>%
    mutate(Chromosome_f = factor(Chromosome, levels = unique(rainfall_points$Chromosome)))
  if (!missing(chromosome)) {
    rainfall_points = dplyr::filter(rainfall_points, Chromosome %in% chromosome)
  }
  if (!missing(zoom_in_region)) {
    rainfall_points = dplyr::filter(
      rainfall_points,
      (
        Chromosome %in% zoom_in_region$chromosome &
          Start_Position >= zoom_in_region$start &
          End_Position <= zoom_in_region$end
      )
    )
  }

  # if user is subsetting by chromosome or zooming in to a specific region, are there any SSM left to plot?
  if (nrow(rainfall_points) == 0) {
    stop("After subsetting to a regions you requested to plot, there are no SSM to display.")
  }

  # label SVs if user wants to overlap this data
  if (!missing(chromosome) & label_sv) {
    sv_chromosome = chromosome
  } else if (!missing(zoom_in_region) & label_sv) {
    sv_chromosome = zoom_in_region$chromosome
  } else if (label_sv) {
    stop(
      "Labeling SV is only supported when a particular chromosome or zoomed region is plotted."
    )
  }

  if (label_sv) {
    message("Getting combined manta + GRIDSS SVs using GAMBLR ...")
    these_sv = get_combined_sv(sample_ids = this_sample_id)
    if ("SCORE" %in% colnames(these_sv)) {
      these_sv = these_sv %>%
        rename("SOMATIC_SCORE" = "SCORE")
    }
    # annotate SV
    these_sv = annotate_sv(these_sv)

    # make SVs a long df with 1 record per SV corresponding to the strand
    sv_to_label =
      melt(
        these_sv %>% select(
          chrom1,
          start1,
          end1,
          chrom2,
          start2,
          end2,
          tumour_sample_id,
          gene,
          partner,
          fusion
        ),
        id.vars = c(
          "tumour_sample_id",
          "gene",
          "partner",
          "fusion",
          "start1",
          "end1",
          "start2",
          "end2"
        ),
        variable.name = "chromosomeN",
        value.name = "Chromosome"
      ) %>%
      dplyr::filter(Chromosome %in% sv_chromosome)

    # are there any SVs on this chromosome/region?
    if (nrow(sv_to_label) > 0) {
      sv_to_label =
        sv_to_label %>%
        melt(
          .,
          id.vars = c(
            "tumour_sample_id",
            "gene",
            "partner",
            "fusion",
            "chromosomeN",
            "Chromosome"
          )
        ) %>%
        group_by(fusion, chromosomeN) %>%
        dplyr::filter(if (grepl("1", chromosomeN))
          variable %in% c("start1", "end1")
          else
            variable %in% c("start2", "end2")) %>%
        dplyr::mutate(variable = gsub("1|2", "", variable)) %>%
        distinct(fusion, Chromosome, variable, .keep_all = TRUE) %>%
        spread(., variable, value) %>%
        dplyr::rename("End_Position" = "end",
                      "Start_Position" = "start") %>%
        ungroup
    } else {
      message(
        "Warning: after subsetting to a regions you requested to plot, there are no SV features to overlap on the final graph."
      )
      label_sv = FALSE
    }

    # when we are plotting region and not whole chromosome, ensure SV is within that region
    if (!missing(zoom_in_region) & label_sv) {
      sv_to_label = dplyr::filter(
        sv_to_label,
        (
          Start_Position >= zoom_in_region$start &
            End_Position <= zoom_in_region$end
        )
      )
      # When we did filtering to start/end for a region, are there any SV to plot?
      if (nrow(sv_to_label) == 0) {
        message(
          "Warning: after subsetting to a regions you requested to plot, there are no SV features to overlap on the final graph."
        )
        label_sv = FALSE
      }
    }

    sv_to_label = sv_to_label %>%
      mutate(Chromosome_f = factor(Chromosome))
  }

  p = ggplot(rainfall_points) +
    geom_point(aes(x = Start_Position, y = IMD, color = Substitution)) +
    scale_color_manual(values = get_gambl_colours("rainfall")) +
    ylab("log(IMD)") +
    theme_Morons() +
    facet_wrap( ~ Chromosome_f, scales = "free_x") +
    ggtitle(this_sample_id) +
    theme(plot.title = element_text(hjust = 0)) # left-align title plot

  if (label_ashm_genes) {
    p = p +
      ggrepel::geom_text_repel(
        data = ashm_regions,
        aes(start, 1, label = gene),
        size = 4,
        segment.size = 0.5,
        segment.color = "#000000",
        force_pull = 0,
        arrow = arrow(length = unit(0.05, "inches"), type = "closed"),
        max.overlaps = 15,
        segment.curvature = 0.25,
        segment.ncp = 4,
        segment.angle = 25
      )
  }

  if (label_sv) {
    p = p +
      geom_vline(
        data = sv_to_label,
        aes(xintercept = Start_Position),
        color = "lightgreen",
        alpha = .7
      ) +
      geom_text(data = sv_to_label,
                aes(End_Position, 15, label = fusion, color = "lightgreen"))
  }

  # show x-axis coordinates if zooming in to a specific region, but not if looking chromosome/genome-wide
  if (missing(zoom_in_region)) {
    p = p + guides(x = "none")
  }

  return(p)
}

gene_mutation_tally = function(maf_df,these_samples_metadata,these_genes,grouping_variable="cohort"){
  meta = dplyr::select(these_samples_metadata,sample_id,{{grouping_variable}})
  maf_filt = dplyr::filter(maf_df,Hugo_Symbol %in% these_genes, Variant_Classification %in% coding_class) %>%
    dplyr::filter(Variant_Classification !="Silent")
  meta_anno = left_join(maf_filt,meta,by=c("Tumor_Sample_Barcode"="sample_id")) %>%
    group_by(Hugo_Symbol,Tumor_Sample_Barcode) %>%
    slice_head() %>%
    ungroup()
  denom = meta %>% group_by(!!sym(grouping_variable)) %>% tally() %>% dplyr::rename(c("total"="n"))
  meta_anno_tally = group_by(meta_anno,Hugo_Symbol,!!sym(grouping_variable)) %>% tally()
  meta_anno_tally = full_join(meta_anno_tally,denom) %>% mutate(frequency=100*n/total)
  return(meta_anno_tally)
}

#' Make a word cloud of gene names from a MAF file based on mutation frequency
#'
#' @param maf_df A MAF-format data frame containing the mutations you want to summarize in a gene word cloud
#' @param these_genes An optional vector of gene symbols (defaults to all lymphoma genes)
#' @param other_genes An optional second vector of gene symbols to include in your cloud in a second colour
#' @param these_genes_colour Specify the hex code of a colour to use for the first set of genes
#' @param other_genes_colour Specify another hex code of a colour to use for the second set of genes
#' @param colour_index Optional named character vector with a name for each gene in these_genes and a colour as the value
#'
#' @return data frame with counts for each gene
#' @export
#' @import wordcloud RColorBrewer
#'
#'
#' @examples
prettyGeneCloud = function(maf_df,these_genes,other_genes,
                           these_genes_colour="#B2DF8A",
                           other_genes_colour="#bc42f5",
                           colour_index){
  if(missing(these_genes)){
    these_genes = pull(lymphoma_genes,Gene)
  }
  #drop genes not in the list then tally only coding variants (by default).
  # TODO: eventually allow an option to collapse samples from the same patient
  if(missing(other_genes)){
    these_genes_maf = dplyr::filter(maf_df,Hugo_Symbol %in% these_genes)
  }else{
    these_genes_maf = dplyr::filter(maf_df,Hugo_Symbol %in% c(these_genes,other_genes))
  }

  #drop non-coding
  these_genes_maf = dplyr::filter(these_genes_maf,Variant_Classification %in% coding_vc)
  these_genes_unique = group_by(these_genes_maf,Hugo_Symbol,Tumor_Sample_Barcode) %>%
    slice_head() %>% ungroup() %>% group_by(Hugo_Symbol) %>% tally()
  print(as.data.frame(head(these_genes_unique,25)))
  if(!missing(other_genes)){
    #assign a colour to each gene list
    these_genes_unique = these_genes_unique %>%
      mutate(this_col=ifelse(Hugo_Symbol %in% these_genes,these_genes_colour,other_genes_colour)) %>% arrange(desc(n))
    wordcloud::wordcloud(these_genes_unique$Hugo_Symbol,these_genes_unique$n,colors=these_genes_unique$this_col,
              ordered.colors = T,scale=c(8,0.3),random.order = F)
  }else{
    if(!missing(colour_index)){
      #use the colours in colour_index to colour the gene names
      if(any(!these_genes %in% names(colour_index))){
        stop("all genes in these_genes must be among the names of colour_index if you specify this variable")
      }
      these_genes_unique$color=colour_index[these_genes_unique$Hugo_Symbol]
      #make cloud with the user-specified colours mapped to the genes
      wordcloud::wordcloud(these_genes_unique$Hugo_Symbol,these_genes_unique$n,random.order=F,ordered.colors=T,colors=these_genes_unique$color)
    }else{
      wordcloud::wordcloud(these_genes_unique$Hugo_Symbol,these_genes_unique$n,random.color=TRUE,colors=RColorBrewer::brewer.pal(12,"Set3"))
    }
  }
  these_genes_unique = arrange(these_genes_unique,n)
  these_genes_unique$Hugo_Symbol = factor(these_genes_unique$Hugo_Symbol,levels=these_genes_unique$Hugo_Symbol)
  return(these_genes_unique)
}


#' Generate a plot of all CN segments.
#'
#' @param region Genomic region for plotting in bed format.
#' @param gene Optional variable, converts gene to region if region not supplied.
#' @param these_samples_metadata GAMBL metadata subset to the cases you want to process (or full metadata)
#' @param type Type of CN segment to be ploitted. Default is gain (CN > 2).
#' @param crop_segments Boolean statement that crops segment by first checking if crop segment is smaller than lef/right distance, then adds or subtracts  crop distance to end/start coordiantes. Default is TRUE.
#' @param sort_by_annotation Sort CN by annotation, default is "pathology".
#' @param crop_distance Crop distance for cropping segments. Default value is 10000000 bp.
#'
#' @return Nothing
#' @import tidyverse
#' @export
#'
#' @examples
#' plot = focal_cn_plot(gene = "BCL2", type = "gain", segment_size = 2, crop_distance = 1000000)
#' plot = focal_cn_plot(region = "chr4:100154645-5488465512", type = "loss", crop_distance = 100000000)
#'
focal_cn_plot = function(region,
                         gene,
                         these_samples_metadata,
                         type = "gain",
                         segment_size = 1,
                         crop_segments = TRUE,
                         sort_by_annotation = c('pathology'),
                         crop_distance = 100000000){

  if(!missing(gene)){
    region = gene_to_region(gene)
    chunks = region_to_chunks(region)
  }else{
    chunks = region_to_chunks(region)
  }
  if(type == "gain"){
    all_not_dip = get_cn_segments(region = region) %>%
      mutate(size = end - start) %>%
      dplyr::filter(CN>2)
  }else{
    all_not_dip = get_cn_segments(region = region) %>%
      mutate(size = end - start) %>%
      dplyr::filter(CN<2)
  }

  #crop start and end if they're further than crop_distance from your region
  all_not_dip = mutate(all_not_dip, left_distance = as.numeric(chunks$start) - start)
  all_not_dip = mutate(all_not_dip, right_distance = end - as.numeric(chunks$end))
  if(crop_segments){
    all_not_dip = mutate(all_not_dip, end = ifelse(right_distance > crop_distance, as.numeric(chunks$end) + crop_distance, end))
    all_not_dip = mutate(all_not_dip, start = ifelse(left_distance > crop_distance, as.numeric(chunks$start) - crop_distance, start))
  }
  all_not_dip = left_join(all_not_dip, these_samples_metadata, by = c("ID" = "sample_id")) %>%
    dplyr::filter(!is.na(pathology))

  all_not_dip = all_not_dip %>%
    arrange(across(all_of(c(sort_by_annotation, "size"))))

  all_not_dip$ID = factor(all_not_dip$ID, levels = unique(all_not_dip$ID))

  ggplot(all_not_dip, aes(x = start, xend = end, y = ID, yend = ID, colour = lymphgen)) +
    geom_vline(aes(xintercept = as.numeric(chunks$start)), alpha = 0.5, colour = get_gambl_colours()[type]) +
    geom_segment(size = segment_size) + theme_cowplot() +
    theme(axis.text.y = element_blank())
}


#' Generate a more visually appealing and flexible lollipop plot.
#'
#' @param maf_df A data frame containing the mutation data (from a MAF).
#' @param gene The gene symbol to plot.
#' @param plot_title Optional (defaults to gene name).
#' @param plot_theme Options: cbioportal(default), blue, simple, nature, nature2, ggplot2, and dark.
#'
#' @return Nothing
#' @export
#' @import g3viz tidyverse
#'
#' @examples
#' pretty_lollipop_plot = (mutation_df, "MYC", "Mutation data for MYC", "blue")
#' pretty_lollipop_plot = (mutation_df, "BCL2")
#'
pretty_lollipop_plot = function(maf_df,
                                gene,
                                plot_title,
                                plot_theme = "cbioportal"){
  if(missing(plot_title)){
    plot_title = gene
  }
  maf_df = maf_df %>%
    dplyr::filter(Hugo_Symbol == gene)

  #use the readMAF function (modified by Ryan) to parse/convert
  maf_df = g3viz::readMAF(maf.df = maf_df)
  chart.options = g3Lollipop.theme(theme.name = plot_theme, title.text = plot_title)
  g3Lollipop(maf_df,
             gene.symbol = gene,
             plot.options = chart.options,
             output.filename = "default_theme")
}


#' Count hypermutated bins and generate heatmaps/cluster the data.
#'
#' @param regions Vector of regions in the format "chr:start-end".
#' @param regions_df Data frame of regions with four columns (chrom,start,end,gene_name).
#' @param these_samples_metadata GAMBL metadata subset to the cases you want to process (or full metadata).
#' @param region_padding How many bases will be added on the left and right of the regions to ensure any small regions are sufficiently covered by bins. Default is  1000.
#' @param metadataColumns What metadata will be shown in the visualization.
#' @param sortByColumns Which of the metadata to sort on for the heatmap.
#' @param expressionColumns Optional variable for retreiving expression values for a specific gene(s).
#' @param orientation Specify the sample orientation, default is sample_rows.
#' @param customColour Optional named list of named vectors for specifying all colours for metadata. Can be generated with map_metadata_to_colours. Default is NULL.
#' @param slide_by How far to shift before starting the next window.
#' @param window_size The width of your sliding window.
#' @param min_count_per_bin Minimum counts per bin, default is 3.
#' @param min_bin_recurrence How many samples a bin must be mutated in to retain in the visualization.
#' @param min_bin_patient How many bins must a patient mutated in to retain in the visualization.
#' @param region_fontsize Fontsize of regions in plot, default is 8ppt.
#' @param cluster_rows_heatmap Optional parameter to enable/disable clustering of each dimension of the heatmap. Default is FALSE.
#' @param cluster_cols_heatmap Optional parameter to enable/disable clustering of each dimension of the heatmap. Default is FALSE.
#' @param show_gene_colours Optional logical argument indicating whether regions should have associated colours plotted as annotation track of heatmap.
#' @param legend_row Fiddle with these to widen or narrow your legend.
#' @param legend_col Fiddle with these to widen or narrow your legend.
#' @param legend_direction Accepts one of "horizontal" (default) or "vertical" to indicate in which direction the legend will be drawn.
#' @param legendFontSize Fontsize of legend in plot, defualt is 10ppt.
#' @param from_indexed_flatfile Set to TRUE to avoid using the database and instead rely on flatfiles (only works for streamlined data, not full MAF details).
#' @param mode Only works with indexed flatfiles. Accepts 2 options of "slms-3" and "strelka2" to indicate which variant caller to use. Default is "slms-3".
#'
#'
#' @return Nothing
#' @import tidyverse ComplexHeatmap
#' @export
#'
#' @examples
#' roi = c("chr1:102502-130210")
#' mut_freq = get_mutation_frequency_bin_matrix(regions = roi, region_padding = 1500, show_gene_colours = TRUE, legendFontSize = 12)
#'
get_mutation_frequency_bin_matrix = function(regions,
                                             regions_df,
                                             these_samples_metadata,
                                             seq_type="genome",
                                             region_padding = 1000,
                                             metadataColumns = c("pathology"),
                                             sortByColumns = c("pathology"),
                                             expressionColumns = c(),
                                             orientation = "sample_rows",
                                             skip_regions=c("MYC", "BCL2", "IGLL5"),
                                             customColour = NULL,
                                             slide_by = 100,
                                             window_size = 500,
                                             min_count_per_bin = 3,
                                             min_bin_recurrence = 5,
                                             min_bin_patient = 0,
                                             region_fontsize = 8,
                                             cluster_rows_heatmap = FALSE,
                                             cluster_cols_heatmap = FALSE,
                                             show_gene_colours = FALSE,
                                             legend_row = 3,
                                             legend_col = 3,
                                             legend_direction = "horizontal",
                                             legendFontSize = 10,
                                             from_indexed_flatfile = TRUE,
                                             mode = "slms-3"){

  if(missing(regions)){
    if(missing(regions_df)){
      regions_df = grch37_ashm_regions #drop MYC and BCL2
      regions_df = grch37_ashm_regions %>%
        dplyr::filter(!gene %in% skip_regions)
    }
    regions = unlist(apply(regions_df, 1, function(x){paste0(x[1], ":", as.numeric(x[2]) - region_padding, "-", as.numeric(x[3]) + region_padding)})) #add some buffer around each
  }
  dfs = lapply(regions, function(x){calc_mutation_frequency_sliding_windows(
    this_region = x, drop_unmutated = TRUE,seq_type=seq_type,
    slide_by = slide_by, plot_type = "none", window_size = window_size,
    min_count_per_bin = min_count_per_bin, return_count = TRUE,
    metadata = these_samples_metadata,
    from_indexed_flatfile = from_indexed_flatfile, mode = mode)})

  all= do.call("rbind", dfs)

  #add a fake bin for one gene and make every patient not mutated in it (to fill gaps)
  fake = these_samples_metadata %>%
    dplyr::select(sample_id) %>%
    mutate(bin = "1_chrN") %>%
    mutate(mutated = 0)

  all = bind_rows(all, fake)
  completed = complete(all, sample_id, bin, fill = list(mutated = 0))
  widened = pivot_wider(completed, names_from = sample_id, values_from = mutated)
  widened_df = column_to_rownames(widened, var = "bin")

  if(length(expressionColumns)>0){
    these_samples_metadata = these_samples_metadata %>%
      mutate(across(all_of(expressionColumns), ~ trim_scale_expression(.x)))
  }
  meta_show = these_samples_metadata %>%
    select(sample_id, all_of(metadataColumns)) %>%
    arrange(across(all_of(sortByColumns))) %>%
    dplyr::filter(sample_id %in% colnames(widened_df)) %>%
    column_to_rownames(var = "sample_id")

  message(paste("starting with", length(colnames(widened_df)), "patients"))
  patients_show = colnames(widened_df)[which(colSums(widened_df)>= min_bin_patient)]
  message(paste("returning matrix with", length(patients_show), "patients"))
  meta_show = dplyr::filter(meta_show, rownames(meta_show) %in% patients_show)
  to_show = widened_df[which(rowSums(widened_df) > min_bin_recurrence), patients_show]

  bin_col_fun = colorRamp2(c(0, 3, 6, 9), c("white", "orange", "red", "purple"))
  to_show_t = t(to_show)
  meta_show_t = meta_show[rownames(to_show_t),]
  lg_cols = get_gambl_colours("lymphgen")
  path_fun = function(x){
    path_cols = get_gambl_colours("pathology")
    lg_cols = get_gambl_colours("lymphgen")

    return(unname(path_cols[x]))
  }
  path_cols = get_gambl_colours("pathology")

  #assign bins back to regions for better annotation
  assign_bins_to_region = function(bin_names, rdf){
    bin_df = data.frame(bin_name = bin_names)

    separated = bin_df %>%
      separate(bin_name, into = c("start", "chrom")) %>%
      mutate(start = as.integer(start)) %>%
      mutate(end = start + 1)

    separated$bin_name = bin_names
    colnames(rdf)[c(1:3)] = c("chrom", "start", "end")
    rdf = mutate(rdf, start = start - 1500) %>%
      mutate(end = end + 1500)

    regions.dt = as.data.table(rdf)

    setkey(regions.dt, chrom, start, end)
    bin.dt = as.data.table(separated)
    setkey(bin.dt, chrom, start, end)
    bin_overlapped = foverlaps(bin.dt, regions.dt, mult = "first") %>%
      as.data.frame() %>%
      select(bin_name, gene) %>%
      column_to_rownames(var = "bin_name")

    return(bin_overlapped)
  }

  #regions_df = grch37_ashm_regions
  if(is.null(customColour)){
    meta_cols = map_metadata_to_colours(metadataColumns, these_samples_metadata = meta_show, as_vector = F)

  }else{
    meta_cols = customColour
  }
  col_fun = circlize::colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  for(exp in expressionColumns){
    meta_cols[[exp]] = col_fun
  }
  bin_annot = assign_bins_to_region(bin_names = colnames(to_show_t), rdf = regions_df)
  heatmap_legend_param = list(title = "Bin value", nrow = legend_row, ncol = legend_row, legend_direction = legend_direction, labels_gp = gpar(fontsize = legendFontSize))

  annotation_legend_param = list(nrow = legend_row, ncol = legend_col, direction = legend_direction, labels_gp = gpar(fontsize = legendFontSize))

  if(orientation == "sample_rows"){
    row_annot = HeatmapAnnotation(df = meta_show, show_legend = T, which = 'row', col = meta_cols, annotation_legend_param = annotation_legend_param)
    if(show_gene_colours){
      col_annot = HeatmapAnnotation(df = bin_annot, show_legend = F, which = 'col', annotation_legend_param = annotation_legend_param)
    }else{
      col_annot = HeatmapAnnotation(value = anno_empty(border = FALSE))
    }
    Heatmap(to_show_t[rownames(meta_show), rownames(bin_annot)],
            cluster_columns = cluster_cols_heatmap,
            cluster_rows = cluster_rows_heatmap,
            col = bin_col_fun,
            bottom_annotation = col_annot,
            left_annotation = row_annot,
            show_row_names = F,
            show_column_names = F,
            column_split = factor(bin_annot$gene),
            column_title_gp = gpar(fontsize = region_fontsize),
            column_title_rot = 90,
            row_title_gp = gpar(fontsize = 10),
            heatmap_legend_param = heatmap_legend_param)
  }else{
    col_annot = HeatmapAnnotation(df = meta_show, show_legend = T, which = 'col', col = meta_cols, annotation_legend_param = annotation_legend_param)
    if(show_gene_colours){
      row_annot = HeatmapAnnotation(df = bin_annot,show_legend = F, which = 'row', annotation_legend_param = annotation_legend_param)
    }else{
      row_annot = rowAnnotation(value = anno_empty(border = FALSE))
    }
    Heatmap(to_show[rownames(bin_annot),rownames(meta_show)],
            show_heatmap_legend = F,
            cluster_columns = cluster_rows_heatmap,
            cluster_rows = cluster_cols_heatmap,
            col = bin_col_fun,
            bottom_annotation = col_annot,
            left_annotation = row_annot,
            show_row_names = F,
            show_column_names = F,
            row_split = factor(bin_annot$gene),
            row_title_gp = gpar(fontsize = region_fontsize),
            row_title_rot = 0,
            column_title_gp = gpar(fontsize = 8),
            heatmap_legend_param = heatmap_legend_param)
  }
}


#' Plot a heatmap comparing the VAF of mutations in T1/T2 pairs.
#'
#' @param maf1 Data frame of simple somatic mutations at time point A.
#' @param maf2 Data frame of simple somatic mutations at time point B.
#' @param vafcolname Name of variable that holds VAF in maf. If not supplied, vaf will be calcualted.
#' @param patient_colname Column name that holds patient name (default is "patient_id").
#' @param these_samples_metadata GAMBL metadata subset to the cases you want to process (or full metadata).
#' @param sortByColumns Which of the metadata to sort on for the heatmap.
#' @param metadata_columns A vector containing the categorical column names you want to plot below.
#' @param gene_orientation Where genes would be plotted. Default is "bottom".
#' @param annotate_zero Indicate a variant that had VAF = 0 in one of the two time points. Default is FALSE.
#' @param genes An optional list of genes to restrict your plot to.
#' @param top_n_genes How many genes to be added to the plot.
#' @param drop_unless_lowvaf Will drop some genes where VAF is low, default is FALSE.
#' @param vaf_cutoff_to_drop Which VAF cut-off value to use when dropping variants before plotting.
#' @param cluster_columns Boolean statement for clustering by columns, defaults to FALSE.
#' @param cluster_rows Boolean statement for clustering by rows, defaults to FALSE.
#'
#' @return Nothing
#' @import tidyverse ComplexHeatmap
#' @export
#'
#' @examples
#' plot = plot_mutation_dynamics_heatmap(maf_df1, maf_df2, "patient_id", gene_orientation = "bottom", annotate_zero = TRUE, top_genes = 100, drop_inless_lowvaf = TRUE, vaf_cutoff_to_drop = 0.04, cluster_rows = TRUE)
#'
plot_mutation_dynamics_heatmap = function(maf1,
                                          maf2,
                                          vafcolname,
                                          patient_colname = "patient_id",
                                          these_samples_metadata,
                                          sortByColumns,
                                          metadata_columns = c("sample_id"),
                                          gene_orientation = "bottom",
                                          annotate_zero = FALSE,
                                          genes,
                                          top_n_genes,
                                          drop_unless_lowvaf = FALSE,
                                          vaf_cutoff_to_drop = 0.05,
                                          cluster_columns = FALSE,
                                          cluster_rows = FALSE){

  if(missing(vafcolname)){
    t1_pair_mafdat = mutate(maf1, vaf = t_alt_count / (t_alt_count + t_ref_count))
    t2_pair_mafdat = mutate(maf2, vaf = t_alt_count / (t_alt_count + t_ref_count))
  }else{
    t1_pair_mafdat[,"vaf"] = t1_pair_mafdat[,vafcolname]
    t2_pair_mafdat[,"vaf"] = t2_pair_mafdat[,vafcolname]
  }
  if(!missing(sortByColumns)){
    these_samples_metadata = arrange(these_samples_metadata, across(sortByColumns))
  }

  t1_pair_mafdat = t1_pair_mafdat %>%
    dplyr::rename("patient_id" = patient_colname)

  t2_pair_mafdat = t2_pair_mafdat %>%
    dplyr::rename("patient_id" = patient_colname)

  these_samples_metadata = these_samples_metadata %>%
    dplyr::rename("patient_id" = patient_colname)

  both_vaf_all = full_join(t1_pair_mafdat, t2_pair_mafdat, by = c("patient_id", "Start_Position")) %>%
    dplyr::select(patient_id, HGVSp_Short.x, Hugo_Symbol.x, Hugo_Symbol.y, Tumor_Sample_Barcode.x,
                  Tumor_Sample_Barcode.y, HGVSp_Short.y, vaf.x, vaf.y, hot_spot.x, hot_spot.y) %>%
    mutate(ANNO = ifelse(is.na(vaf.x), HGVSp_Short.y, HGVSp_Short.x)) %>%
    mutate(GENE = ifelse(is.na(vaf.x), Hugo_Symbol.y, Hugo_Symbol.x)) %>%
    mutate(VAF1 = ifelse(is.na(vaf.x), 0, vaf.x)) %>%
    mutate(VAF2 = ifelse(is.na(vaf.y), 0, vaf.y)) %>%
    mutate(hot_spot.x = ifelse(is.na(hot_spot.x), 0, 1)) %>%
    mutate(hot_spot.y = ifelse(is.na(hot_spot.y), 0, 1)) %>%
    mutate(HOTSPOT = ifelse(hot_spot.x == 1 | hot_spot.y == 1, 1, 0)) %>%
    dplyr::filter(ANNO != "") %>%
    mutate(Mutation = paste(GENE, ANNO, sep = "_")) %>%
    dplyr::select(patient_id, GENE, Mutation, VAF1, VAF2, ANNO, HOTSPOT)

  if(drop_unless_lowvaf){
    both_vaf_all = dplyr::filter(both_vaf_all, VAF1 < vaf_cutoff_to_drop | VAF2 < vaf_cutoff_to_drop)
  }
  if(!missing(genes)){
    both_vaf_all = dplyr::filter(both_vaf_all, GENE %in% genes)
  }
  both_vaf_all = mutate(both_vaf_all, unique_id = paste(patient_id, Mutation, sep = "_")) %>%
    mutate(fold_change = log(VAF2 + 0.1)-log(VAF1 + 0.1)) %>%
    group_by(patient_id, GENE) %>%
    mutate(Number = paste(GENE, row_number(), sep = "_"))

  print(head(both_vaf_all))
  h = both_vaf_all %>%
    select(patient_id, Number, fold_change) %>%
    pivot_wider(id_cols = patient_id, names_from = Number, values_from = fold_change) %>%
    column_to_rownames("patient_id")

  both_vaf_all = mutate(both_vaf_all, minVAF = ifelse(VAF1 < VAF2, VAF1, VAF2))
  zeroes = both_vaf_all %>%
    select(patient_id, Number, minVAF) %>%
    pivot_wider(id_cols = patient_id, names_from = Number, values_from = minVAF) %>%
    column_to_rownames("patient_id")

  hotspots = both_vaf_all %>%
    select(patient_id, Number, HOTSPOT) %>%
    pivot_wider(id_cols = patient_id, names_from = Number, values_from = HOTSPOT) %>%
    column_to_rownames("patient_id")

  zeroes[is.na(zeroes)] = 0.001
  hotspots[is.na(hotspots)] = -1
  hotspots[hotspots == 0] = -1
  hotspots[hotspots == 1] = 0
  zeroes[zeroes > 0] = 1
  print(head(hotspots))

  these_samples_metadata_rn = dplyr::filter(these_samples_metadata, patient_id %in% rownames(h)) %>%
    select(all_of(c("patient_id", metadata_columns))) %>%
    column_to_rownames("patient_id")

  la = HeatmapAnnotation(df = as.data.frame(these_samples_metadata_rn), which = "row")
  ta = HeatmapAnnotation(df = as.data.frame(these_samples_metadata_rn), which = "column")

  h[is.na(h)] = 0
  cs = colSums(zeroes)
  ordered = names(cs[order(cs)])
  if(!missing(top_n_genes)){
    print(head(ordered, top_n_genes))
    genes = ordered[c(1:top_n_genes)]
  }
  col_fun = colorRamp2(c(0, 1), c("white", "red"))

  if(gene_orientation == "bottom"){
    if(annotate_zero){

    }else{
      H = Heatmap(h[rownames(these_samples_metadata_rn),], cluster_rows = F, cluster_columns = F, left_annotation = la)
    }
  }else{
    if(!missing(top_n_genes)){
      these_zeroes = t(zeroes[rownames(these_samples_metadata_rn), genes])
      these_zeroes = t(hotspots[rownames(these_samples_metadata_rn), genes])
      to_show = t(h[rownames(these_samples_metadata_rn), genes])

    }else{
      these_zeroes = t(zeroes[rownames(these_samples_metadata_rn),])
      these_zeroes = t(hotspots[rownames(these_samples_metadata_rn),])
      to_show = t(h[rownames(these_samples_metadata_rn),])
    }
    if(annotate_zero){
      H = Heatmap(to_show, layer_fun = function(j, i, x, y, width, height, fill) {
                                                v = pindex(these_zeroes, i, j)
                                                l = v == 0
                                                grid.points(x[l], y[l], pch = 16, size = unit(1, "mm"), gp = gpar(col = "white"))
    },
      cluster_columns = cluster_columns, cluster_rows = cluster_rows, bottom_annotation = ta)
      print("HERE")
    }else{
      H = Heatmap(t(h[rownames(these_samples_metadata_rn),]), cluster_rows = cluster_rows, cluster_columns = cluster_columns, bottom_annotation = ta)
    }
  }
  return(H)
}


#' INTERNAL FUNCTION called by plot_sample_circos and get_mutation_frequency_bin_matrix, not meant for out-of-package usage.
#' Assign a colour palette to metadata columns automatically and consistently.
#'
#' @param metadataColumns Names of the metadata columns to assign colours for.
#' @param these_samples_metadata Metadata for just the samples you need colours for.
#' @param column_alias A list of column_names with values indicating what gambl colour scheme they should use (e.g. pos_neg, pathology, lymphgen).
#' @param as_vector Boolean statement that is set to TRUE per default.
#' @param verbose Set to TRUE to enable verbose mode (debugging messages.
#' @param annoAlpha Optional alpha to apply to annotation colours.
#'
#' @return Either a vector or list of colours.
#' @import dplyr ggsci
#'
#' @examples
#' all_cols=map_metadata_to_colours(legend_metadata_columns,these_meta,verbose=T)
#' all_cols=map_metadata_to_colours(c("lymphgen","pathology","genetic_subgroup"),these_samples_metadata = all_meta,column_alias=list("nothing"="FL"),as_vector = F)
#'
map_metadata_to_colours = function(metadataColumns,
                                   these_samples_metadata,
                                   column_alias = list(),
                                   as_vector = TRUE,
                                   verbose = FALSE,
                                   annoAlpha = 1){

  clinical_colours = ggsci::get_ash("clinical")
  all_gambl_colours = get_gambl_colours()
  colours = list()
  colvec = c()


  aliases = c(colour_aliases, column_alias)
  for(column in metadataColumns){
    this_value = these_samples_metadata[[column]]
  options = this_value
    if(verbose){
      print(">>>>>>>")
      message("finding colour for", this_value)
      print("<<<<<<<")
    }
    if(column %in% names(aliases)){
      key = aliases[[column]]
      if(verbose){
        print(paste("using alias to look up colours for", column))
        message(paste("using", key, "for", column))
      }
      these = get_gambl_colours(classification = key)
      colours[[column]] = these
      colvec = c(colvec, these[this_value])
      if(verbose){
        message("adding:", these[this_value])
      }
    }else if(column == "sex"){
      these = get_gambl_colours("sex", alpha = annoAlpha)
      if(!"NA" %in% names(these)){
        these= c(these, "NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec, these[this_value])
      message("adding:", these[this_value])
    }else if(sum(levels(options) %in% names(clinical_colours)) == length(levels(options))){

      #we have a way to map these all to colours!
      if(verbose){
        message(paste("found colours for", column, "in clinical"))
      }
      these = clinical_colours[levels(options)]
      if(!"NA" %in% names(these)){
        these = c(these,"NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec, these[this_value])
    }else if(("positive" %in% options | "POS" %in% options) & length(options)<4){
      if(verbose){
        print("using pos_neg")
      }

      these = get_gambl_colours("pos_neg", alpha = annoAlpha)
      these = these[levels(options)]

      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec, these[this_value])
    }else if("GCB" %in% options){
      these = get_gambl_colours("COO", alpha = annoAlpha)
      if(!"NA" %in% names(these)){
        these = c(these,"NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec,these)
    }else if(column %in% c("pathology")){
      these = get_gambl_colours(column, alpha = annoAlpha)

      if(!"NA" %in% names(these)){
        these = c(these,"NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec, these)
    }else if(grepl("lymphgen", column, ignore.case = TRUE)){
      these = get_gambl_colours("lymphgen", alpha = annoAlpha)
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec, these)
    }else if(column == "HMRN"){
      these = get_gambl_colours("hmrn", alpha = annoAlpha)
      if(!"NA" %in% names(these)){
        these= c(these,"NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec, these)
    }else if(sum(levels(options) %in% names(all_gambl_colours)) == length(levels(options))){
      if(verbose){
        message(paste("found colours for", column, "in all_gambl_colours"))
      }
      these = all_gambl_colours[levels(options)]
      if(!"NA" %in% names(these)){
        these= c(these,"NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec, these)
    }else if(length(levels(options)) > 15){

      these = rainbow(length(levels(options)), alpha = annoAlpha)
      names(these) = levels(options)

      colours[[column]] = these
      colvec = c(colvec, these)
    }else{
      these = blood_cols[sample(c(1:length(blood_cols)), size = length(levels(options)))]
      names(these) = levels(options)
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
      colvec = c(colvec, these)
    }
  }
  if(as_vector){
    return(colvec)
  }
  return(colours)
}


#' Plot a sample-centric circos overview.
#'
#' @param this_sample_id Sample ID for the sample to plot.
#' @param sv_df Optional data frame of SVs (default is to use the database).
#' @param cnv_df Optional data frame of CNVs (default is to use the database).
#' @param ssm_df Optional data frame of SSMs (default is to use the database).
#' @param include_sv Default TRUE.
#' @param include_cnv Default TRUE.
#' @param include_ssm Defaul FALSE.
#' @param legend_metadata_columns Column names from meta data
#' @param legend_metadata_names List of meta data names to be plotted.
#' @param chrom_list List of chromosomes to be plotted. If not stated, chr1-22+X will bes used.
#' @param label_genes Gene labels (df, list or what type?)
#' @param auto_label_sv Default is FALSE
#'
#' @return Nothing
#' @export
#' @import circlize ComplexHeatmap
#'
#' @examples
#' this_samp = "13-38657_tumorB"
#' GAMBLR::plot_sample_circos(this_sample_id=this_samp,legend_metadata_columns = c("pathology","lymphgen","COO_consensus","DHITsig_consensus","bcl2_ba","myc_ba"),legend_metadata_names = c("pathology","LymphGen","COO","DHITsig","BCL2","MYC"),auto_label_sv = TRUE,chrom_list = c("chr2","chr3","chr8","chr14","chr18"))
#'
plot_sample_circos = function(this_sample_id,
                              sv_df,
                              cnv_df,
                              ssm_df,
                              include_sv = TRUE,
                              include_ssm = FALSE,
                              legend_metadata_columns,
                              legend_metadata_names = c(),
                              include_cnv = TRUE,
                              chrom_list,
                              label_genes,
                              auto_label_sv = FALSE){

  add_cnv = function(cnv_df){
    bed = data.frame(cnv_df[,c("chrom", "start", "end", "log.ratio")])
    colnames(bed) = c("chr", "start", "end", "value1")
    col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    circos.genomicTrackPlotRegion(bed, stack = TRUE, panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = col_fun(value), border = NA, posTransform = NULL, ...)
      i = getI(...)
      cell.xlim = get.cell.meta.data("cell.xlim")
    }, bg.border = NA)
  }
  if(missing(cnv_df)){
    cnv_df = get_sample_cn_segments(this_sample_id = this_sample_id, with_chr_prefix = TRUE)
  }
  if(missing(chrom_list)){

  #should we add chr list for males as well? Sex could then also be added as a paramter for the function were the appropiate chr list is called if not stated?
   chrom_list = paste0("chr", c(1:22,"X"))
  }
  if(!missing(label_genes)){
    gene_bed = grch37_gene_coordinates %>%
      dplyr::filter(gene_name %in% label_genes) %>%
      dplyr::select(chromosome, start, end, gene_name) %>%
      dplyr::mutate(chromosome = paste0("chr", chromosome))
  }
  if(missing(sv_df)){
    sv_df = get_manta_sv(with_chr_prefix = TRUE) %>%
      dplyr::filter(tumour_sample_id == this_sample_id)

  }
  sv_df = sv_df %>%
    dplyr::filter(CHROM_A %in% chrom_list) %>%
    dplyr::filter(CHROM_B %in% chrom_list)

  if(auto_label_sv){

    #annotate oncogene SVs and label them
    annotated_sv  = annotate_sv(sv_df, with_chr_prefix = TRUE) %>%
      dplyr::filter(!is.na(partner)) %>%
      dplyr::filter(tumour_sample_id == this_sample_id)

    these_oncogenes = unique(pull(annotated_sv, gene))
    these_partners = unique(pull(annotated_sv, partner))
    if("IGH" %in% these_partners){
      these_partners = c(these_partners, "IGHV3-62")
    }
    anno_bed1 = annotated_sv %>%
      dplyr::select(chrom1, start1, end1, tumour_sample_id)

    anno_bed2 = annotated_sv %>%
      dplyr::select(chrom2, start2, end2, tumour_sample_id)

    colnames(anno_bed1) = c("chrom", "start", "end", "sample_id")
    colnames(anno_bed2) = c("chrom", "start", "end", "sample_id")

    bed_mut_partner = grch37_partners %>%
      dplyr::filter(gene %in% these_partners) %>%
      mutate(chrom = paste0("chr", chrom))

    bed_mut_onco = grch37_oncogene %>%
      dplyr::filter(gene %in% these_oncogenes) %>%
      mutate(chrom = paste0("chr", chrom))

    bed_mut = bind_rows(bed_mut_partner, bed_mut_onco)
    print(bed_mut)
  }
  bed1 = sv_df %>%
    dplyr::select(CHROM_A, START_A, END_A, tumour_sample_id)

  bed2 = sv_df %>%
    dplyr::select(CHROM_B, START_B, END_B, tumour_sample_id)

  colnames(bed1) = c("chrom", "start", "end", "sample_id")
  colnames(bed2) = c("chrom", "start", "end", "sample_id")
  circos.clear()
  circos.initializeWithIdeogram(chromosome.index = chrom_list)
  add_cnv(cnv_df)
  circos.genomicLink(bed1, bed2, col = "#bdbdc1")
  if(!missing(label_genes)){
    circos.genomicLabels(gene_bed, labels.column = "gene_name")
  }
  if(auto_label_sv){
    circos.genomicLink(anno_bed1, anno_bed2,col = 'red')
    circos.genomicLabels(bed_mut, labels.column = "gene")
  }
  text(c(0.75, 0.75), this_sample_id, cex = 0.8)
  if(!missing(legend_metadata_columns)){
    samp_meta = get_gambl_metadata() %>%
      dplyr::filter(sample_id == this_sample_id)

    these_meta = samp_meta[legend_metadata_columns]
    these_cols = get_gambl_colours()
    vals = as.character(these_meta)
    names = colnames(these_meta)

    all_cols = map_metadata_to_colours(legend_metadata_columns, these_meta, verbose = T)

    cols = all_cols[vals]
    print(cols)
    if(length(legend_metadata_names) == length(legend_metadata_columns)){
      for(i in c(1:length(vals))){
        if(!legend_metadata_names[i] == ""){
          vals[i] = paste(legend_metadata_names[i], vals[i])
        }
      }
    }

    lgd_discrete = Legend(labels = vals, title_position = "topleft", legend_gp = gpar(fill = cols))
    draw(lgd_discrete, x = unit(3, "mm"), y = unit(1, "npc") - unit(5, "mm"), just = c("left", "top"))
  }

  #continuous
  lgd_cnv = Legend(at = c(-2, -1, 0, 1, 2),
                   col_fun = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                   title_position = "topleft",
                   title = "log\nratio")

  draw(lgd_cnv, x = unit(1, "npc") - unit(2, "mm"), y = unit(4, "mm"), just = c("right", "bottom"))
}


#' Make an oncoplot that is pretty using ComplexHeatmap. The metadata is expected to follow the structure and column naming used in GAMBL.
#' If you provide your own non-GAMBL samples and metadata, you must include at least the following columns with these names.
#' The first one should match the Tumor_Sample_Barcode in the MAF object or onco_matrix you provide.
#' sample_id, pathology
#'
#' @param maftools_obj A maftools object containing the mutations you want to plot.
#' @param onco_matrix_path Provide a path to an onco_matrix file instead of a MAF object if the former is unavailable (this limits functionality a bit).
#' @param genes An optional list of genes to restrict your plot to.
#' @param include_noncoding List of non-coding regions to be included, default is NULL. Specify like this: include_noncoding=list("NFKBIZ" = c("3'UTR"), "HNRNPH1" = "Splice_Region")
#' @param keepGeneOrder Set to TRUE if you want to preserve the gene order specified.
#' @param keepSampleOrder Set to TRUE if you want to preserve the sample order specified.
#' @param highlightHotspots Set to TRUE to highlight hot spots. Default is FALSE.
#' @param these_samples_metadata Data frame containing metadata for your samples.
#' @param metadataColumns A vector containing the categorical column names you want to plot below.
#' @param numericMetadataColumns A vector containing the numeric columns you want to plot below.
#' @param expressionColumns Optional variable for retreiving expression values for a specific gene(s).
#' @param numericMetadataMax A numeric vector of cutoffs to apply to numeric columns above.
#' @param sortByColumns A vector containing the column names you want to sort columns (patients) on.
#' @param removeNonMutated Set to TRUE to drop unmutated cases.
#' @param minMutationPercent Only genes mutated in more than minMutationPercent % patients will be included.
#' @param fontSizeGene Font size for gene labels (default 6).
#' @param annoAlpha Optional alpha to apply to annotation colours.
#' @param mutAlpha Optional alpha to apply to mutation colours.
#' @param recycleOncomatrix Set to TRUE most of the time to reuse the oncomatrix saved by maftools.
#' @param box_col Colour of boxes for outlining mutations (can be problematic with larger oncoprints).
#' @param metadataBarHeight Optional argument to adjust the height of bar with annotations. The default is 1.5.
#' @param metadataBarFontsize Optional argument to control for the font size of metadata annotations. The default is 5.
#' @param hideTopBarplot Optional argument for removing top bar plot. Default value is TRUE.
#' @param hideSideBarplot Optional argument for removing side bar plot. Default value is FALSE.
#' @param splitColumnName Optional argument to indicate which metadata column to split on. Default is set to pathology.
#' @param splitGeneGroups Split genes into groups for better seperation (between different gene-groups) in prettyOncoplot.
#' @param legend_row Fiddle with these to widen or narrow your legend.
#' @param legend_col Fiddle with these to widen or narrow your legend.
#' @param showTumorSampleBarcode Optional argument for showing tumor barcode. Default is FALSE.
#' @param groupNames optional vector of group names to be displayed above heatmap. Should be the same length as the number of groups that will be shown. Default is NULL (no labels).
#' @param verbose Set to TRUE to enable verbose mode (debugging messages.
#' @param hide_annotations Hide annotations for specifc ashms. argument takes a list with annotations.
#' @param annotate_specific_genes Optional argument, specifying whether the features should be labelled according to their significance in one of the pathologies. Default is FALSE (no annotation).
#' @param this_forest_object If annotate_specific_genes is specified, this arguments takes the output of GAMBLR::prettyForestPlot directly to determine the annotations.
#' @param custom_colours Provide named vector (or named list of vectors) containing custom annotation colours if you do not want to use standartized pallette.
#' @param legend_direction Direction of lgend, defualt is "horizontal".
#' @param ylim Limit for y-axis.
#' @param legend_position Position of legend, default is "bottom".
#' @param annotation_row Row for annotations, default is 2.
#' @param annotation_col Column for annotations, default is 1.
#' @param legendFontSize Font size for legend, default is 10.
#'
#' @return Nothing
#' @export
#' @import ComplexHeatmap grid
#'
#' @examples
#' prettyOncoplot(maftools_obj = maf_obj,genes = bl_genes,
#' these_samples_metadata = extra_meta,
#' metadataColumns = c("pathology","COO_consensus", "cluster_name", "lymphgen","EBV_status_inf", "manta_BCL6_sv"),
#' hide_annotations = c(some_ashm,"lymphgen","COO_consensus", "pathology","manta_BCL6_sv"),
#' expressionColumns = c("IRF4",some_ashm),
#' sortByColumns = c("IRF4"),
#' keepGeneOrder = FALSE,splitGeneGroups = groups,
#' splitColumnName = "cluster_name",
#' metadataBarHeight = 2.5,metadataBarFontsize = 6,fontSizeGene = 8,
#' recycleOncomatrix = TRUE,removeNonMutated = FALSE)
#'
prettyOncoplot = function(maftools_obj,
                          onco_matrix_path,
                          genes,
                          include_noncoding = NULL,
                          keepGeneOrder = FALSE,
                          keepSampleOrder = TRUE,
                          highlightHotspots = FALSE,
                          these_samples_metadata,
                          metadataColumns,
                          numericMetadataColumns,
                          expressionColumns = c(),
                          numericMetadataMax,
                          sortByColumns,
                          arrange_descending = FALSE,
                          removeNonMutated = FALSE,
                          minMutationPercent,
                          mutAlpha = 1,
                          recycleOncomatrix = FALSE,
                          splitColumnName,
                          splitGeneGroups,
                          showTumorSampleBarcode = FALSE,
                          groupNames,
                          hide_annotations,
                          annotate_specific_genes = FALSE,
                          this_forest_object = NULL,
                          custom_colours = NULL,
                          hideTopBarplot = TRUE,
                          tally_all_mutations = FALSE,
                          tally_all_mutations_max = 1000,
                          hideSideBarplot = FALSE,
                          box_col = NA,
                          annoAlpha = 1,
                          legend_direction = "horizontal",
                          ylim = NULL,
                          legend_position = "bottom",
                          legend_row = 3,
                          legend_col = 3,
                          metadataBarHeight = 1.5,
                          metadataBarFontsize = 5,
                          legendFontSize = 10,
                          fontSizeGene = 6,
                          annotation_row = 2,
                          annotation_col = 1,
                          verbose = FALSE){

  patients = pull(these_samples_metadata, sample_id)
  #ensure patients not in metadata get dropped up-front to ensure mutation frequencies are accurate
  if(!recycleOncomatrix & missing(onco_matrix_path)){
    onco_matrix_path = "onco_matrix.txt"
  #order the data frame the way you want the patients shown
    maf_patients = unique(as.character(maftools_obj@data$Tumor_Sample_Barcode))
    if(any(!maf_patients %in% patients)){
      extra = maf_patients[which(!maf_patients %in% patients)]
      patients = maf_patients[which(maf_patients %in% patients)]
      n_drop = length(extra)
      message(paste(n_drop, "patients are not in your metadata, will drop them from the data before displaying"))
      maftools_obj = subsetMaf(maf = maftools_obj, tsb = patients)
    }
    if(missing(genes)){
      #check that our MAFtools object only contains samples in the supplied metadata
      genes = maftools::getGeneSummary(x = maftools_obj)[order(MutatedSamples, decreasing = TRUE)][,.(Hugo_Symbol, MutatedSamples)]
      colnames(genes)[2] = "mutload"
      totSamps = as.numeric(maftools_obj@summary[3, summary])
      genes[,fractMutated := mutload / totSamps]
      genes = genes[fractMutated * 100 >= minMutationPercent, Hugo_Symbol]
      lg = length(genes)
      message(paste("creating oncomatrix with", lg, "genes"))
      om = maftools:::createOncoMatrix(m = maftools_obj, g = genes, add_missing = TRUE)
      mat_origin = om$oncoMatrix
      tsbs = levels(maftools:::getSampleSummary(x = maftools_obj)[,Tumor_Sample_Barcode])
      print(paste("numcases:", length(tsbs)))

      if(!removeNonMutated){
        tsb.include = matrix(data = 0, nrow = nrow(mat_origin), ncol = length(tsbs[!tsbs %in% colnames(mat_origin)]))
        colnames(tsb.include) = tsbs[!tsbs %in% colnames(mat_origin)]
        rownames(tsb.include) = rownames(mat_origin)
        mat_origin = cbind(mat_origin, tsb.include)
      }
      write.table(mat_origin, file = onco_matrix_path, quote = F, sep = "\t")
    }else{
      if(any(duplicated(genes))){
        stop("There are duplicated elements in the provided gene list (@param genes). Please ensure only unique entries are present in this list.")
      }
      om = maftools:::createOncoMatrix(m = maftools_obj, g = genes, add_missing = TRUE)
      mat_origin = om$oncoMatrix
      tsbs = levels(maftools:::getSampleSummary(x = maftools_obj)[,Tumor_Sample_Barcode])
      print(paste("numcases:",length(tsbs)))
      print(paste("numgenes:",length(mat_origin[,1])))
      if(!removeNonMutated){
        tsb.include = matrix(data = 0, nrow = nrow(mat_origin), ncol = length(tsbs[!tsbs %in% colnames(mat_origin)]))
        colnames(tsb.include) = tsbs[!tsbs %in% colnames(mat_origin)]
        rownames(tsb.include) = rownames(mat_origin)
        mat_origin = cbind(mat_origin, tsb.include)
      }
      write.table(mat_origin, file = onco_matrix_path, quote = F, sep = "\t")
    }
  }
  if(missing(onco_matrix_path)){
    onco_matrix_path = "onco_matrix.txt"
  }
  if(!missing(numericMetadataColumns)){
    message(paste0("The column(s) ", numericMetadataColumns, " specified both in metadata and numeric metadata. Plotting as numeric values..."))
    metadataColumns = metadataColumns[!metadataColumns %in% numericMetadataColumns]
  }
  patients = pull(these_samples_metadata, sample_id)
  #because the way MAFtools writes this file out is the absolute worst for compatability
  old_style_mat = read.table(onco_matrix_path, sep = "\t", stringsAsFactors = FALSE)
  mat = read.table(onco_matrix_path, sep = "\t", header = TRUE, check.names = FALSE, row.names = 1, fill = TRUE, stringsAsFactors = F, na.strings = c("NA", ""))
  colnames(old_style_mat) = colnames(mat)
  mat = old_style_mat
  mat[mat==0]=""
  #add the noncoding mutations to this if requested (only for genes and types specified)
  if(length(include_noncoding) > 0){
    all_genes_df = data.frame(Hugo_Symbol = rownames(mat))
    all_samples_df = data.frame(Tumor_Sample_Barcode = colnames(mat))
    for(gene in names(include_noncoding)){
      for(this_vc in unname(include_noncoding[[gene]])){
        message(paste(gene, "and", this_vc))
        these_samples = dplyr::filter(maftools_obj@maf.silent,Hugo_Symbol == gene & Variant_Classification == this_vc) %>%
          dplyr::select(Tumor_Sample_Barcode, Variant_Classification) %>%
          unique() %>%
          pull(Tumor_Sample_Barcode)
        for(samp in these_samples){
          if(samp %in% colnames(mat)){
            if(mat[gene, samp] == ""){
              mat[gene, samp] = this_vc
            }else{
              mat[gene, samp] = paste0(this_vc, ";", mat[gene, samp])
            }
          }
        }
      }
    }
  }
  #annotate hot spots if necessary
  if(missing(metadataColumns)){
    message("you should name at least one metadata column to show as an annotation. Defaulting to pathology")
    metadataColumns = c("pathology")
  }
  if(missing(genes)){
    genes = rownames(mat)
  }
  col = get_gambl_colours("mutation", alpha = mutAlpha)
  mat[mat == 0]=""
  patients_kept = patients[which(patients %in% colnames(mat))]
  patients_dropped = patients[which(!patients %in% colnames(mat))]
  if(verbose){
    message("====DROPPED=====")
    message(patients_dropped)
  }
  genes_kept = genes[which(genes %in% rownames(mat))]
  genes_dropped = genes[which(!genes %in% maftools_obj@gene.summary$Hugo_Symbol)]
  for (g in genes_dropped) {
    maftools_obj@gene.summary = dplyr::add_row(maftools_obj@gene.summary, Hugo_Symbol = g)
  }
  maftools_obj@gene.summary <- maftools_obj@gene.summary %>% replace(is.na(.), 0)
  if(!missing(minMutationPercent)){
    if(! onco_matrix_path == "onco_matrix.txt"){

      warning("mintMutationPercent option is not available when you provide your own oncomatrix. Feel free to implement this if you need it")
      return()
    }
    mutation_counts <- maftools_obj@gene.summary %>%
      select(Hugo_Symbol, MutatedSamples)

    numpat = length(patients)
    mutation_counts = mutate(mutation_counts, percent_mutated = 100 * MutatedSamples / numpat)
    genes_keep = mutation_counts %>%
      dplyr::filter(percent_mutated >= minMutationPercent) %>%
      pull(Hugo_Symbol)

    genes_kept = genes[genes %in% genes_keep]
  }
  mat = mat[,patients_kept]
  mat = mat[which(rownames(mat) %in% genes_kept),]
  spacing = 0
  height_scaling = 1
  alter_fun = list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w - unit(spacing, "pt"), h * height_scaling,
                gp = gpar(fill = "#e6e6e6", col = box_col))
    },
    RNA = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = "#F2ED36", col = box_col))
    },
    `3'UTR` = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = "#F2ED36", col = box_col))
    },
    Intron = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), 0.75* h*height_scaling,
                gp = gpar(fill = col["Nonsense_Mutation"], col = box_col))
    },
    #big blue
    Nonsense_Mutation = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = "#D8A7CA", col = box_col))
    },
    Splice_Site = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["Splice_Site"], col = box_col))
    },
    Splice_Region = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["Splice_Region"], col = box_col))
    },
    Nonstop_Mutation = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["Nonstop_Mutation"], col = box_col))
    },
    Translation_Start_Site = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["Translation_Start_Site"], col = box_col))
    },
    In_Frame_Ins = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["In_Frame_Ins"], col = box_col))
    },
    In_Frame_Del = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["In_Frame_Del"], col = box_col))
    },
    #all frame shifts will be the same colour, magenta
    Frame_Shift_Del = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["Frame_Shift_Del"], col = box_col))
    },
    Frame_Shift_Ins = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["Frame_Shift_Ins"], col = box_col))
    },
    #big red
    Multi_Hit = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["Multi_Hit"], col = box_col))
    },
    #small green
    Missense_Mutation = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), h*height_scaling,
                gp = gpar(fill = col["Missense_Mutation"], col = box_col))
    },
    hot_spot = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), (height_scaling/5)*h,
                gp = gpar(fill = "white", col = box_col))
    },
    Silent = function(x, y, w, h) {
      grid.rect(x, y, w-unit(spacing, "pt"), (height_scaling/5)*h,
                gp = gpar(fill = col["Silent"], col = box_col))
    }
  )
  #automagically assign colours for other metadata columns.
  #TO DO: convert the loop below into a "map_metadata_to_colours" function HAS THIS BEEN RESOLVED?
  blood_cols = get_gambl_colours("blood", alpha = annoAlpha)
  colours = list()
  clinical_colours = ggsci::get_ash("clinical")
  all_gambl_colours = get_gambl_colours()
  for(column in metadataColumns){
    these_samples_metadata[[column]] = factor(these_samples_metadata[[column]], levels = unique(these_samples_metadata[[column]]))
    options = these_samples_metadata %>%
      arrange(column) %>%
      dplyr::filter(!is.na(column)) %>%
      pull(column) %>%
      unique()

    options = options[!is.na(options)]
    if(verbose){
      print(">>>>>>>")
      print(levels(options))
      print("<<<<<<<")
    }
    if(column == "sex"){
      these = get_gambl_colours("sex", alpha = annoAlpha)
      these = these[levels(options)]
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }else if(sum(levels(options) %in% names(clinical_colours)) == length(levels(options))){
      #we have a way to map these all to colours!
      if(verbose){
        message(paste("found colours for", column, "here"))
      }
      these = clinical_colours[levels(options)]
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }else if(("positive" %in% options | "POS" %in% options | "yes" %in% options) & length(options) < 4){
      if(verbose){
        print("using pos_neg")
      }
      these = get_gambl_colours("pos_neg", alpha = annoAlpha)
      these = these[levels(options)]
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }else if("GCB" %in% options){
      these = get_gambl_colours("COO", alpha = annoAlpha)
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }else if(column %in% c("pathology")){
      these = get_gambl_colours(column, alpha = annoAlpha)
      if(!"NA" %in% names(these)){
        these = c(these,"NA" = "white")
      }
      colours[[column]] = these
    }else if(grepl("lymphgen", column, ignore.case = TRUE)){
      these = get_gambl_colours("lymphgen", alpha = annoAlpha)
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }else if(column == "HMRN"){
      these = get_gambl_colours("hmrn", alpha = annoAlpha)
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }else if(sum(levels(options) %in% names(all_gambl_colours)) == length(levels(options))){
      if(verbose){
        message(paste("found colours for", column, "in all_gambl_colours"))
      }
      these = all_gambl_colours[levels(options)]
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }else if(length(levels(options)) > 15){
      these = rainbow(length(levels(options)), alpha = annoAlpha)
      names(these) = levels(options)
        colours[[column]] = these
    }else{
      these = blood_cols[sample(c(1:length(blood_cols)), size = length(levels(options)))]
      names(these) = levels(options)
      if(!"NA" %in% names(these)){
        these = c(these, "NA" = "white")
      }
      colours[[column]] = these
    }
  }
  if (! is.null(custom_colours)){
    colours = custom_colours
  }

  if(highlightHotspots){
    hot_samples = dplyr::filter(maftools_obj@data, hot_spot == TRUE & Hugo_Symbol %in% genes) %>%
      dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode) %>%
      mutate(mutated = "hot_spot") %>%
      unique()

    all_genes_df = data.frame(Hugo_Symbol = rownames(mat))
    all_samples_df = data.frame(Tumor_Sample_Barcode = colnames(mat))
    hs = left_join(all_samples_df, hot_samples)
    hot_mat = hs %>%
      pivot_wider(names_from = "Tumor_Sample_Barcode", values_from = "mutated") %>%
      left_join(all_genes_df,.) %>%
      column_to_rownames("Hugo_Symbol") %>%
      as.matrix()

    #annotate hotspots in matrix
    for (i in colnames(mat)){
      mat[genes, i][!is.na(hot_mat[genes, i])] = paste0(mat[genes, i][!is.na(hot_mat[genes, i])], ";", hot_mat[genes, i][!is.na(hot_mat[genes, i])])
    }
    colours[["hot_spots"]] = c("hot_spot" = "magenta")
  }
  if(verbose){
    print(colours) #eventually get rid of this once the bugs are gone
  }
  if(!missing(numericMetadataColumns)){
    metadata_df = dplyr::filter(these_samples_metadata, sample_id %in% patients_kept) %>%
      column_to_rownames("sample_id") %>%
      dplyr::select(all_of(c(metadataColumns, numericMetadataColumns, expressionColumns)))

    if(!missing(numericMetadataMax)){
      max_list = setNames(numericMetadataMax, numericMetadataColumns)

      metadata_df = metadata_df %>%
        mutate(across(names(max_list), ~ ifelse(.x > max_list[[cur_column()]], max_list[[cur_column()]], .x)))
    }
  }else{
    metadata_df = dplyr::filter(these_samples_metadata, sample_id %in% patients_kept) %>%
      column_to_rownames("sample_id") %>%
      dplyr::select(all_of(c(metadataColumns, expressionColumns)))
  }
  metadata_df = metadata_df %>%
    mutate(across(all_of(expressionColumns), ~ trim_scale_expression(.x)))
  if(!missing(sortByColumns)){
    if (arrange_descending) {
      metadata_df = arrange(metadata_df, across(sortByColumns, desc))
    } else {
      metadata_df = arrange(metadata_df, across(sortByColumns))
    }
    patients_kept = rownames(metadata_df)
  }
  if(verbose){
    print(genes_kept)
  }
  col_fun = circlize::colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  for(exp in expressionColumns){
    colours[[exp]] = col_fun
  }
  if(missing(splitColumnName)){
    column_split = rep("", length(patients_kept))
  }else{
    column_split = factor(metadata_df[patients_kept, splitColumnName])
  }
  if(missing(splitGeneGroups)){
    row_split = rep("", length(genes))
  }else{
    row_split = factor(splitGeneGroups[genes], levels = unique(splitGeneGroups[genes]))
  }
  if(!missing(groupNames)){
    column_title = groupNames
  }else{
    column_title = NULL
  }
  if(keepGeneOrder){
    gene_order = genes
  }else{
    gene_order = NULL
  }
  if(missing(hide_annotations)){
    show_legend = rep(TRUE, length(colnames(metadata_df)))
  }else{
    show_legend = rep(TRUE, length(colnames(metadata_df)))
    names(show_legend) = colnames(metadata_df)
    show_legend[hide_annotations] = FALSE
  }
  if(missing(sortByColumns)){
    column_order = NULL
  }else{
    column_order = patients_kept
  }
  heatmap_legend_param = list(title = "Alterations",
                         at = c("RNA", "3'UTR" , "Nonsense_Mutation", "Splice_Site","Splice_Region", "Nonstop_Mutation", "Translation_Start_Site",
                         "In_Frame_Ins", "In_Frame_Del", "Frame_Shift_Ins", "Frame_Shift_Del", "Multi_Hit", "Missense_Mutation", "Silent", "hot_spot"),
                         labels = c("RNA", "3'UTR", "Nonsense Mutation", "Splice Site","Splice Region", "Nonstop Mutation", "Translation Start Site",
                         "In Frame Insertion", "In Frame Deletion", "Frame Shift Insertion", "Frame Shift Deletion",
                         "Multi Hit", "Missense Mutation", "Silent", "Hotspot"),
                         nrow = annotation_row, ncol = annotation_col,
                         legend_direction = legend_direction,
                         labels_gp = gpar(fontsize = legendFontSize))
  if(hideTopBarplot){
    top_annotation = NULL
  }else{
    tally_mutations = maftools_obj@data %>%
      dplyr::filter(Tumor_Sample_Barcode %in% patients_kept) %>%
      group_by(Tumor_Sample_Barcode) %>%
      summarize(n_mutations = n()) %>%
      ungroup %>%
      arrange(match(Tumor_Sample_Barcode, patients_kept)) %>%
      select(n_mutations) %>%
      mutate(n_mutations = ifelse(n_mutations > tally_all_mutations_max,
                                  tally_all_mutations_max,
                                  n_mutations))

    if(is.null(ylim) & ! tally_all_mutations){
      top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot())

    }else if (!is.null(ylim) & ! tally_all_mutations){
      top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(ylim=ylim))

    } else if (is.null(ylim) & tally_all_mutations) {
      top_annotation = columnAnnotation(" " = anno_barplot(tally_mutations))
    } else if (! is.null(ylim) & tally_all_mutations) {
      top_annotation = columnAnnotation(" " = anno_barplot(tally_mutations, ylim=ylim))
    }
  }

  # Handle right annotation for specific genes
  if (annotate_specific_genes & is.null(this_forest_object)) {
    message("WARNING: You requested right annotation, but forgot to provide output of GAMBLR::prettyForestPlot")
    message("No right annotation will be drawn.")
    right_annotation = NULL
  } else if (annotate_specific_genes) {

    these_comparisons = this_forest_object$mutmat$comparison %>% levels

    enrichment_label =
      mat[intersect(genes, genes_kept),patients_kept] %>%
      rownames_to_column("gene") %>%
      select(gene) %>%
      left_join(this_forest_object$fisher %>% select(gene, estimate, q.value)) %>%
      mutate("Enriched in" = case_when(
        estimate == "Inf" & q.value <= 0.1 ~ these_comparisons[1],
        estimate == "-Inf" & q.value <= 0.1 ~ these_comparisons[2],
        is.na(estimate) ~ "NA",
        estimate<=1 & q.value <= 0.1 ~ these_comparisons[2],
        estimate > 1 & q.value <= 0.1 ~ these_comparisons[1],
        TRUE ~ "Both"
      )) %>%
      pull("Enriched in")

    right_annotation = rowAnnotation(" " = enrichment_label,
                             col = list(" " = c(get_gambl_colours()[these_comparisons], Both = "#ACADAF", "NA" = "#000000")),
                             simple_anno_size = unit(metadataBarHeight, "mm"),
                             annotation_legend_param =
                               list(title = "Enriched in",
                                    nrow=legend_row,
                                    ncol = legend_col,
                                    direction=legend_direction,
                                    labels_gp = gpar(fontsize = legendFontSize)))
  } else {
    right_annotation = NULL
  }

  ch = ComplexHeatmap::oncoPrint(mat[intersect(genes, genes_kept),patients_kept],
                                 alter_fun = alter_fun,
                                 top_annotation = top_annotation,
                                 right_annotation = right_annotation,
                                 col = col,
                                 row_order = gene_order,
                                 column_order = column_order,
                                 column_labels = NULL,
                                 show_column_names = showTumorSampleBarcode,
                                 column_split = column_split,
                                 column_title = column_title,
                                 row_title = NULL,
                                 row_split = row_split[intersect(genes, genes_kept)],
                                 heatmap_legend_param = heatmap_legend_param,
                                 row_names_gp = gpar(fontsize = fontSizeGene),
                                 pct_gp = gpar(fontsize = fontSizeGene),
                                 bottom_annotation = ComplexHeatmap::HeatmapAnnotation(df = metadata_df,
                                                                                       show_legend = show_legend,
                                                                                       col = colours,
                                                                                       simple_anno_size = unit(metadataBarHeight, "mm"),
                                                                                       gap = unit(0.25 * metadataBarHeight, "mm"),
                                                                                       annotation_name_gp = gpar(fontsize = metadataBarFontsize),
                                                                                       annotation_legend_param = list(nrow = legend_row,
                                                                                                                      col_fun = col_fun,
                                                                                                                      ncol = legend_col,
                                                                                                                      direction = legend_direction,
                                                                                                                      labels_gp = gpar(fontsize = legendFontSize))))

    draw(ch, heatmap_legend_side = legend_position, annotation_legend_side = legend_position)
}


#' Display 2 prettyOncoplots side-by-side.
#'
#' `prettyCoOncoplot` returns ggplot-compatible figure of 2 prettyOncoplots side-by-side.
#'
#' This function will generate a graphic displaying 2 oncoplots side-by-side. Optionally user can
#' annotate each oncoplot with it's own title that will be displayed at the top. All the arguments
#' recognized by prettyOncoplot are supported and can be specified when calling this function.
#' For both oncoplots the same specified parameters will be applied (e.g. genes to display, split columns,
#' font size, top annotation etc). If the provided argument is not recognized by prettyOncoplot,
#' it will be discarded. If you want a specific order of oncoplots on the left and right, please
#' ensure the argument `comparison_column` is a factor with first level being the group
#' you want to be plotted on the left side. For developers: new arguments added to prettyOncoplot in the future
#' are expected to be out-of-the-box compatible with this function nd would not need code modifications.
#'
#' @param maf Required argument. A maftools object containing the mutations you want to plot on both oncoplots.
#' @param metadata Required argument. A data.frame with metadata for both oncoplots.
#' @param comparison_values Optional: If the comparison column contains more than two values or is not a factor, specify a character vector of length two in the order you would like the factor levels to be set, reference group first.
#' @param comparison_column Required: the name of the metadata column containing the comparison values.
#' @param label1 Optional argument. Label to be shown as a title for the oncoplot #1.
#' @param label2 Optional argument. Label to be shown as a title for the oncoplot #2.
#'
#' @return A ggplot object with 2 oncoplots side-by-side.
#' @export
#' @import ComplexHeatmap ggpubr maftools
#'
#' @examples
#' ssm=get_coding_ssm(limit_cohort = c("BL_Adult", "BL_Pediatric"))
#' ssm=maftools::read.maf(ssm)
#' meta=get_gambl_metadata() %>% dplyr::filter(cohort %in% c("BL_Adult", "BL_Pediatric"))
#' prettyCoOncoplot(maf=ssm,
#'     metadata = meta,
#'     comparison_column = "cohort",
#'     include_noncoding = NULL,
#'     minMutationPercent = 0,
#'     genes=c("MYC", "TET2", "TP53", "DDX3X", "ID3"),
#'     metadataColumns=c("pathology", "EBV_status_inf", "pairing_status", "cohort"),
#'     splitColumnName="EBV_status_inf",
#'     metadataBarHeight=10,
#'     fontSizeGene=12,
#'     metadataBarFontsize=10,
#'     label1="Adult",
#'     label2="Pediatric")
#'
prettyCoOncoplot <-   function(maf,
                               metadata,
                               comparison_column,
                               comparison_values,
                               label1,
                               label2,
                               ...) {
    # check for required arguments
    required = c("maf", "metadata", "comparison_column")

    defined = names(as.list(match.call())[-1])

    if (any(!required %in% defined)) {
      stop("Please provide mutation data and metadata for 2 pretty Oncoplots with specified comparison_column.")
    }

    #If no comparison_values are specified, derive the comparison_values from the specified comparison_column
    if(missing(comparison_values)){
      if(class(metadata[[comparison_column]]) == "factor"){
        comparison_values = levels(metadata[[comparison_column]])
      } else {
        comparison_values = unique(metadata[[comparison_column]])
      }
    }

    #Ensure there are only two comparison_values
    {
      if(length(comparison_values) != 2)
        stop(paste0("Your comparison must have two values. \nEither specify comparison_values as a vector of length 2 or subset your metadata so your comparison_column has only two unique values or factor levels."))
    }

    #Subset the metadata to the specified comparison_values and the maf to the remaining sample_ids
    meta1 = metadata[metadata[[comparison_column]] %in% comparison_values[1], ]
    meta2 = metadata[metadata[[comparison_column]] %in% comparison_values[2], ]

    # Subset maf to only samples in the comparison values
    ssm1 = maftools::subsetMaf(maf,tsb=pull(meta1, Tumor_Sample_Barcode))
    ssm2 = maftools::subsetMaf(maf,tsb=pull(meta2, Tumor_Sample_Barcode))

    # Arguments to pass into prettyOncoplot
    oncoplot_args = list(...)
    # Discard any arguments not supported by prettyOncoplot
    oncoplot_args = oncoplot_args[names(oncoplot_args) %in% intersect(names(oncoplot_args),
                                                                      formalArgs(prettyOncoplot))]
    # Build oncoplot No1
    op1 = do.call(prettyOncoplot, c(
      list(
        maftools_obj = ssm1,
        these_samples_metadata = meta1
      ),
      oncoplot_args
    ))
    # convert it to ggplot object
    op1 = grid.grabExpr(draw(op1), width = 10, height = 17)
    # if user provided annotation label, place it as a name for oncoplot No1
    if (!missing(label1)) {
      op1 = annotate_figure(op1,
                            top = text_grob(label1,
                                            face = "bold",
                                            size = 20))
    }
    # Build oncoplot No2
    op2 = do.call(prettyOncoplot, c(
      list(
        maftools_obj = ssm2,
        these_samples_metadata = meta2
      ),
      oncoplot_args
    ))
    # convert it to ggplot object
    op2 = grid.grabExpr(draw(op2), width = 10, height = 17)
    # if user provided annotation label, place it as a name for oncoplot No2
    if (!missing(label2)) {
      op2 = annotate_figure(op2,
                            top = text_grob(label2,
                                            face = "bold",
                                            size = 20))
    }
    # arrange 2 oncoplots together side by side
    p = ggarrange(op1,
                  op2,
                  ncol = 2,
                  nrow = 1)

    return(p)
  }


#' Generate a colourful multi-panel overview of hypermutation in regions of interest across many samples.
#'
#' @param regions_bed Bed file with chromosome coordinates, should contain columns chr, start, end, name (with these exact names).
#' @param regions_to_display Optional vector of names from default regions_bed to use.
#' @param exclude_classifications Optional argument for excluding specific classifications from a metadeta file.
#' @param metadata A metadata file already subsetted and arranged on the order you want the samples vertically displayed.
#' @param custom_colours Provide named vector (or named list of vectors) containing custom annotation colours if you do not want to use standartized pallette.
#' @param classification_column Optional. Override default column for assigning the labels used for colouring in the figure.
#' @param maf_data Either a maf loaded from disk or from the database using a get_ssm function.
#'
#' @return Nothing
#' @export
#' @import tidyverse DBI RMariaDB
#'
#' @examples
#' my_plot = ashm_multi_rainbow_plot(regions_bed = "my_bed.bed",
#'                                   regions_to_display = "chr3",
#'                                   exclude_classification = "some-classification",
#'                                   metadata = "my_metadata",
#'                                   custom_colours = c("#onecolour, "#anothercolour", "athirdcolour"),
#'                                   classification_column = "lymphgen",
#'                                   maf_data = my_maf)
#'
ashm_multi_rainbow_plot = function(regions_bed,
                                   regions_to_display,
                                   exclude_classifications,
                                   metadata,
                                   seq_type,
                                   custom_colours,
                                   classification_column = "lymphgen",
                                   maf_data){

  table_name = config::get("results_tables")$ssm
  db = config::get("database_name")
  #get the mutations for each region and combine
  #regions_bed should contain chr, start, end, name (with these exact names)
  if(missing(metadata)){
    metadata = get_gambl_metadata()
    meta_arranged = arrange(metadata, pathology_rank, lymphgen)
  }else{
    meta_arranged = metadata #assume the user already arranged it the way they wanted
  }
  if(!missing(exclude_classifications)){
    meta_arranged = dplyr::filter(meta_arranged, !get(classification_column) %in% exclude_classifications)
  }
  if(missing(regions_bed)){
    regions_bed = grch37_ashm_regions
    regions_bed = mutate(regions_bed, regions = paste0(chr_name, ":", hg19_start, "-", hg19_end))
    regions_bed = mutate(regions_bed, name = paste0(gene, "-", region))
  }else{
    regions_bed = mutate(regions_bed,regions=paste0(chr,":",start,"-",end))
    #if name column is missing, add it
    if(!"name" %in% colnames(regions_bed))
    {
      regions_bed$name = regions_bed$regions
    }
  }
  print(regions_bed)
  names = pull(regions_bed, name)
  names = c(names, "NFKBIZ-UTR", "MAF", "PAX5", "WHSC1", "CCND1",
                   "FOXP1-TSS1", "FOXP1-TSS2", "FOXP1-TSS3", "FOXP1-TSS4",
                   "FOXP1-TSS5", "BCL6", "IGH", "IGL", "IGK", "PVT1", "BCL2") #add some additional regions of interest
  regions = pull(regions_bed, regions)
  regions = c(regions,"chr3:101578214-101578365", "chr16:79627745-79634622", "chr9:36898851-37448583", "chr4:1867076-1977887", "chr11:69451233-69460334", "chr3:71623481-71641671",
                      "chr3:71532613-71559445", "chr3:71343345-71363145", "chr3:71167050-71193679", "chr3:71105715-71118362", "chr3:187406804-188522799","chr14:106144562-106344765",
                      "chr22:23217074-23250428","chr2:89073691-89320640", "chr8:128774985-128876311","chr18:60982124-60990180")
  regions_bed = data.frame(regions = regions, names = names)
  regions_bed = dplyr::filter(regions_bed, names %in% regions_to_display)
  regions = pull(regions_bed, regions)
  names = pull(regions_bed, names)
  if(missing(maf_data)){
    region_mafs = lapply(regions, function(x){get_ssm_by_region(region = x, streamlined = TRUE,seq_type = seq_type)})
  }else{
    region_mafs = lapply(regions, function(x){get_ssm_by_region(region = x, streamlined = TRUE, maf_data = maf_data)})
  }
  tibbled_data = tibble(region_mafs, region_name = names)
  unnested_df = tibbled_data %>%
    unnest_longer(region_mafs)

  unlisted_df = mutate(unnested_df, start = region_mafs$Start_Position, sample_id = region_mafs$Tumor_Sample_Barcode) %>%
    dplyr::select(start,sample_id, region_name)

  meta_arranged = meta_arranged %>% mutate_if(is.factor, as.character)
  meta_arranged = meta_arranged %>% mutate(classification = factor(!!sym(classification_column)))


  muts_anno = left_join(unlisted_df, meta_arranged)
  muts_first = dplyr::select(muts_anno, start, region_name) %>%
    group_by(region_name) %>%
    arrange(start) %>%
    dplyr::filter(row_number() == 1)

  eg = expand_grid(start = pull(muts_first, start), sample_id = pull(meta_arranged, sample_id))
  eg = left_join(eg, muts_first)

  #concatenate expanded frame of points with original mutation data
  real_and_fake = bind_rows(unlisted_df, eg)
  muts_anno = left_join(real_and_fake, meta_arranged)

  muts_anno$sample_id = factor(muts_anno$sample_id, levels = meta_arranged$sample_id)

  if(!missing(regions_to_display)){
    muts_anno = dplyr::filter(muts_anno, region_name %in% regions_to_display)
  }
  #make the plot
  p = muts_anno %>%
        ggplot() +
        geom_point(aes(x = start, y = sample_id, colour = classification), alpha = 0.4, size = 0.6) +
        labs(title = "", subtitle = "", x = "", y = "Sample") +
        theme_Morons() +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(1,1,1,1, "cm"), title = element_blank(), plot.subtitle = element_blank(), axis.title.x = element_blank()) +
        facet_wrap(~region_name, scales = "free_x") +
        guides(color = guide_legend(reverse = TRUE,
                                    override.aes = list(size = 3),
                                    title={{ classification_column }}))

  if(! missing(custom_colours)){
    # ensure only relevant color keys are present
    custom_colours = custom_colours[intersect(names(custom_colours), pull(meta_arranged[,classification_column]))]
    p = p +
        scale_colour_manual(values = custom_colours)
  }

  print(p)

}


#' Create a genome-wide copy number plot for one sample and (optionally) display mutation VAF.
#'
#' @param this_sample The sample_id for the sample to plot.
#' @param just_segments Specify whether only the segments will be plotted (instead of mutation VAF).
#' @param coding_only Optional. Set to TRUE to restrict to plotting only coding mutations.
#' @param one_chrom Subset plot to one chromosome.
#' @param genes_to_label Optional. Provide a list of genes to label (if mutated). Can only be used with coding_only (see above).
#' @param from_flatfile If set to true the function will use flatfiles instead of the database.
#' @param use_augmented_maf Boolean statement if to use augmented maf, default is FALSE.
#' @param add_chr_prefix If TRUE, "chr" prefix will be added to chr column. Default is FALSE.
#'
#' @return Nothing
#' @export
#' @import tidyverse DBI RMariaDB
#'
#' @examples
#' cnv_vaf_plot = copy_number_vaf_plot(this_sample = "some-sample-name",
#'                                     just_segments = TRUE,
#'                                     coding_only = FALSE,
#'                                     one_chrom = "chr2",
#'                                     genes_to_label = "MYC",
#'                                     from_flatfile = FALSE,
#'                                     use_augmented_maf = FALSE,
#'                                     add_chr_prefix = TRUE)
#'
copy_number_vaf_plot = function(this_sample,
                                just_segments = FALSE,
                                coding_only = FALSE,
                                one_chrom,
                                genes_to_label,
                                from_flatfile = FALSE,
                                use_augmented_maf = FALSE,
                                add_chr_prefix = FALSE){

  chrom_order = factor(c(1:22, "X"))
  if(add_chr_prefix){
    chrom_order = c(1:22, "X")
    chrom_order = factor(unlist(lapply(chrom_order, function(x){paste0("chr", x)})))
  }
  cn_colours = get_gambl_colours(classification = "copy_number")
  maf_and_seg = assign_cn_to_ssm(this_sample = this_sample, coding_only = coding_only, from_flatfile = from_flatfile, use_augmented_maf = use_augmented_maf)
  vaf_cn_maf = maf_and_seg[["maf"]]
  vaf_cn_maf = mutate(vaf_cn_maf, CN = case_when(LOH == "1" & CN == 2 ~ "nLOH", TRUE ~ as.character(CN)))
  if(!missing(one_chrom)){
    vaf_cn_maf = dplyr::filter(vaf_cn_maf, Chromosome == one_chrom)
  }
  if(just_segments){
    cn_seg = maf_and_seg[["seg"]]
    cn_seg = mutate(cn_seg, CN_segment = as.numeric(CN), CN = as.character(CN))
    print(head(cn_seg))
    if(!missing(one_chrom)){
      cn_seg = dplyr::filter(cn_seg, Chromosome %in% one_chrom)
    }
    ggplot(cn_seg) +
      geom_segment(data = cn_seg, aes(x = Start_Position, xend = End_Position, y = CN_segment, yend = CN_segment, colour = CN)) +
      facet_wrap(~Chromosome, scales = "free_x") +
      scale_colour_manual(values = cn_colours) +
      theme_minimal() +
      guides(color = guide_legend(reverse = TRUE, override.aes = list(size = 3)))
  }else{
    if(coding_only){
      if(missing(genes_to_label)){
        p = mutate(vaf_cn_maf, vaf = t_alt_count / (t_ref_count + t_alt_count)) %>%
        ggplot() +
          geom_point(aes(x = Start_Position, y = vaf, colour = CN), alpha = 0.6, size = 2) +
          scale_colour_manual(values = cn_colours) +
          facet_wrap(~factor(Chromosome, levels = chrom_order), scales = "free_x") +
          theme_minimal() +
          guides(color = guide_legend(reverse = TRUE, override.aes = list(size = 3)))
        p + ggtitle(this_sample)
      }else{
        #label any mutations that intersect with our gene list
        plot_genes = vaf_cn_maf %>%
          dplyr::filter(Hugo_Symbol %in% my_genes)

        p = mutate(vaf_cn_maf, vaf = t_alt_count / (t_ref_count + t_alt_count)) %>%
        ggplot() +
          geom_point(aes(x = Start_Position, y = vaf, colour = CN), size = 2) +
          geom_text(data = plot_genes, aes(x = Start_Position, y = 0.8, label = Hugo_Symbol), size = 3, angle = 90) +
          scale_colour_manual(values = cn_colours) +
          facet_wrap(~factor(Chromosome, levels = chrom_order), scales = "free_x") +
          ylim(c(0,1)) +
          theme_minimal() +
          guides(color = guide_legend(reverse = TRUE, override.aes = list(size = 3)))
        p + ggtitle(this_sample)
      }
    }else{
      p = mutate(vaf_cn_maf, vaf = t_alt_count / (t_ref_count + t_alt_count)) %>%
      ggplot() +
        geom_point(aes(x = Start_Position, y = vaf, colour = CN), alpha = 0.6, size = 0.2) +
        scale_colour_manual(values = cn_colours) +
        facet_wrap(~factor(Chromosome, levels = chrom_order), scales = "free_x") +
        theme_minimal() +
        guides(color = guide_legend(reverse = TRUE, override.aes = list(size = 3)))
      p + ggtitle(this_sample)
    }
  }
}


#TODO migrate viz/plotting functions that don't directly rely on the database to a separate file, DONE?
#' Make a rainbow plot of all mutations in a region, ordered and coloured by metadata.
#'
#' @param mutations_maf A data frame containing mutations (MAF format) within a region of interest (i.e. use the get_ssm_by_region).
#' @param metadata should be a data frame with sample_id as a column that should match Tumor_Sample_Barcode in the database.
#' @param exclude_classifications Optional argument for excluding specific classifications from a metadeta file.
#' @param drop_unmutated Boolean argument for removing unmutated sample ids in mutated cases.
#' @param classification_column The name of the metadata column to use for ordering and colouring samples.
#' @param bed Optional data frame specifying the regions to annotate (required columns: start, end, name).
#' @param region Genomic region for plotting in bed format.
#' @param custom_colours Provide named vector (or named list of vectors) containing custom annotation colours if you do not want to use standartized pallette.
#' @param hide_ids Boolean argument, if TRUE, ids will be removed.
#'
#' @return ggplot2 object
#' @export
#' @import tidyverse DBI RMariaDB dbplyr
#'
#' @examples
#' #basic usage
#' region = "chr6:90975034-91066134"
#' metadata = get_gambl_metadata()
#' plot = ashm_rainbow_plot(metadata = metadata, region = region)
#' #advanced usages
#' mybed = data.frame(start=c(128806578,128805652,128748315), end=c(128806992,128809822,128748880), name=c("TSS","enhancer","MYC-e1"))
#' ashm_rainbow_plot(mutations_maf=my_mutations,metadata=my_metadata,bed=mybed)
#'
ashm_rainbow_plot = function(mutations_maf,
                             metadata,
                             exclude_classifications,
                             drop_unmutated = FALSE,
                             classification_column,
                             bed,
                             region,
                             custom_colours,
                             hide_ids = TRUE){

  table_name = config::get("results_tables")$ssm
  db=config::get("database_name")
  if(!missing(region)){
    region = gsub(",", "", region)
    split_chunks = unlist(strsplit(region, ":"))
    chromosome = split_chunks[1]
    startend = unlist(strsplit(split_chunks[2], "-"))
    qstart = as.numeric(startend[1])
    qend = as.numeric(startend[2])
    if(missing(mutations_maf)){
      mutations_maf = get_ssm_by_region(region = region, streamlined = TRUE,from_indexed_flatfile = T)
    }else{
      #ensure it only contains mutations in the region specified
      mutations_maf = get_ssm_by_region(region = region, streamlined = TRUE, maf_data = mutations_maf)
    }
  }
  if(!missing(classification_column)){
    meta_arranged = arrange(metadata, pathology_rank, lymphgen)
    if(!missing(exclude_classifications)){
      meta_arranged = dplyr::filter(meta_arranged,!get(classification_column) %in% exclude_classifications)
    }
  }else{
    classification_column = "lymphgen"
    meta_arranged = metadata
  }
  mutation_positions = mutations_maf %>%
    dplyr::select(Tumor_Sample_Barcode, Start_Position) %>%
    as.data.frame()

  mutated_cases = pull(mutation_positions, Tumor_Sample_Barcode) %>%
    unique()

  if(drop_unmutated){
    meta_arranged = meta_arranged %>%
      dplyr::filter(sample_id %in% mutated_cases)
  }
  #add a fake mutation at the start position for each sample to ensure every sample shows up
  fake_mutations = data.frame(Tumor_Sample_Barcode = pull(metadata, sample_id), Start_Position = qstart - 1000)
  mutation_positions = rbind(mutation_positions, fake_mutations)

  meta_arranged$classification = meta_arranged[[classification_column]] %>%
    as.factor()

  muts_anno = dplyr::left_join(mutation_positions, meta_arranged, by = c("Tumor_Sample_Barcode" = "sample_id")) %>%
    subset(!is.na(classification))

  muts_anno$sample_id = factor(muts_anno$Tumor_Sample_Barcode, levels = unique(meta_arranged$sample_id))

  if(missing(custom_colours)){
    p = ggplot(muts_anno) +
      geom_point(aes(x = Start_Position, y = sample_id, colour = classification), alpha = 0.4)
  }else{
    p = ggplot(muts_anno) +
      geom_point(aes(x = Start_Position, y = sample_id, colour = classification), alpha = 0.4) +
      scale_colour_manual(values = custom_colours)
  }
  if(missing(bed)){
    p + guides(color = guide_legend(reverse = TRUE, override.aes = list(size = 3)))
  }else{
    bed = bed %>%
      mutate(size = end - start) %>%
      mutate(midpoint = start + size / 2)
    height = length(unique(meta_arranged$sample_id)) + 10
    p = p + geom_rect(data = bed, aes(xmin = start, xmax = end, ymin = 0, ymax = height + 20), alpha = 0.1) +
      geom_text(data = bed, aes(x = midpoint, y = height, label = name), size = 2.5, angle = 90) +
      guides(color = guide_legend(reverse = TRUE, override.aes = list(size = 3)))
  }

  p = p +
    labs(y = "Sample") +
    theme_Morons() +
    theme(plot.margin = margin(1,1,1,1, "cm"), title = element_blank(), plot.subtitle = element_blank(), axis.title.x = element_blank())

  if(hide_ids){
    p = p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }else{
    p = p + theme(axis.text.y = element_text(size = 5))
  }
  return(p)
}


  #' This function doesn't do anything yet
#'
#' @param mafs TODO
#' @param sample_id TODO
#' @param genes TODO
#' @param show_noncoding TODO
#' @param detail TODO
#'
#' @return
#' @export
#' @import tidyverse
#'
plot_multi_timepoint = function(mafs,
                                sample_id,
                                genes,
                                show_noncoding = FALSE,
                                detail){

  tp = c("A","B","C")
  title = paste(sample_id, detail, sep = "\n")
  i = 1
  for (i in c(1:length(mafs))){
    maf_file = mafs[i]
    time_point = tp[i]
    print(paste(maf_file, time_point))
  }
  if(length(mafs)==2){
    A.maf = fread_maf(mafs[1])
    B.maf = fread_maf(mafs[2])

    A.maf = A.maf %>%
      dplyr::select(c(Hugo_Symbol, Chromosome, Start_Position, End_Position, Variant_Classification, Tumor_Sample_Barcode, HGVSp_Short, t_ref_count, t_alt_count)) %>%
      mutate(VAF = t_alt_count/(t_ref_count + t_alt_count)) %>%
      mutate(time_point = 1) %>%
      mutate(coord = paste(Chromosome, Start_Position, sep = ":"))

    B.maf = B.maf %>%
      dplyr::select(c(Hugo_Symbol, Chromosome, Start_Position, End_Position, Variant_Classification, Tumor_Sample_Barcode, HGVSp_Short, t_ref_count, t_alt_count)) %>%
      mutate(VAF = t_alt_count/(t_ref_count + t_alt_count)) %>%
      mutate(time_point = 2) %>%
      mutate(coord = paste(Chromosome, Start_Position, sep = ":"))

    all.maf = rbind(A.maf, B.maf)
    if(show_noncoding){
      coding.maf = subset_regions(all.maf, shm_regions)
    }
    else{
      coding.maf = dplyr::filter(all.maf, !Variant_Classification %in% c("Silent", "RNA", "IGR", "Intron", "5'Flank", "3'Flank", "5'UTR"))
      #keep certain 3' UTR mutations, toss the rest
      coding.maf = dplyr::filter(coding.maf, (Variant_Classification == "3'UTR" & Hugo_Symbol == "NFKBIZ") | (Variant_Classification != "3'UTR"))
    }
    A.rows = which(coding.maf$time_point==1)
    A.zero.coords = pull(unique(coding.maf[which(coding.maf$time_point==1 & VAF == 0), "coord"]))
    A.zero = dplyr::filter(coding.maf, coord %in% A.zero.coords)

    B.rows = which(coding.maf$time_point==2)
    B.zero.coords = pull(unique(coding.maf[which(coding.maf$time_point==2 & VAF == 0), "coord"]))
    B.zero = dplyr::filter(coding.maf, coord %in% B.zero.coords)

    coding.maf$category = "trunk"
    coding.maf[which(coord %in% A.zero.coords), "category"] = "not-A"
    coding.maf[which(coord %in% B.zero.coords), "category"] = "not-B"

    #actually this is changed in eitehr direction, not just gained
    just_gained_lg_all = dplyr::filter(coding.maf, Hugo_Symbol %in% genes & category != "trunk" )

    #just_gained_lg = filter(coding.maf, Hugo_Symbol %in% genes & category != "trunk" & time_point !=2) %>%
    #  mutate(time_point = time_point +0.4)

    just_gained_lg = dplyr::filter(coding.maf, Hugo_Symbol %in% genes & category != "trunk" & time_point == 2 & VAF > 0) %>%
      mutate(time_point = time_point +0.4)

    print(just_gained_lg)

    just_trunk = dplyr::filter(coding.maf, Hugo_Symbol %in% genes & category == "trunk" & time_point ==1) %>%
      mutate(time_point = time_point-0.4)

    ggplot(coding.maf, aes(x = time_point, y = VAF, group = coord, colour = category)) +
      geom_point() +
      geom_line(alpha = 0.5) +
      geom_text_repel(data = just_gained_lg, aes(label = Hugo_Symbol), size = 4,segment.linetype = 0) +
      geom_text_repel(data = just_trunk, aes(label = Hugo_Symbol), size = 4, segment.linetype = 0) +
      ggtitle(title) +
      theme_minimal()
    return(TRUE)
  }
  if(length(mafs)==3){
    A.maf = fread_maf(mafs[1])
    B.maf = fread_maf(mafs[2])
    C.maf = fread_maf(mafs[3])

    A.maf = A.maf %>%
      dplyr::select(c(Hugo_Symbol, Chromosome, Start_Position, End_Position, Variant_Classification, Tumor_Sample_Barcode, HGVSp_Short, t_ref_count, t_alt_count))  %>%
      mutate(VAF = t_alt_count/(t_ref_count + t_alt_count)) %>%
      mutate(time_point = 1) %>%
      mutate(coord = paste(Chromosome, Start_Position, sep = ":"))

    B.maf = B.maf %>%
      dplyr::select(c(Hugo_Symbol, Chromosome, Start_Position, End_Position, Variant_Classification, Tumor_Sample_Barcode, HGVSp_Short, t_ref_count, t_alt_count))  %>%
      mutate(VAF = t_alt_count/(t_ref_count + t_alt_count)) %>%
      mutate(time_point = 2) %>%
      mutate(coord = paste(Chromosome, Start_Position,sep = ":"))

    C.maf = C.maf %>%
      dplyr::select(c(Hugo_Symbol, Chromosome, Start_Position, End_Position, Variant_Classification, Tumor_Sample_Barcode, HGVSp_Short, t_ref_count, t_alt_count))  %>%
      mutate(VAF = t_alt_count/(t_ref_count + t_alt_count)) %>%
      mutate(time_point = 3) %>%
      mutate(coord = paste(Chromosome, Start_Position, sep = ":"))

    all.maf = rbind(A.maf, B.maf, C.maf)
    if(show_noncoding){
      coding.maf = subset_regions(all.maf, shm_regions)
    }
    else{
      coding.maf = dplyr::filter(all.maf, !Variant_Classification %in% c("Silent", "RNA", "IGR", "Intron", "5'Flank", "3'Flank", "5'UTR"))
      #keep certain 3' UTR mutations, toss the rest
      coding.maf = dplyr::filter(coding.maf, (Variant_Classification == "3'UTR" & Hugo_Symbol == "NFKBIZ") | (Variant_Classification != "3'UTR"))
    }

    A.rows = which(coding.maf$time_point==1)
    A.zero.coords = pull(unique(coding.maf[which(coding.maf$time_point==1 & VAF == 0), "coord"]))
    A.zero = dplyr::filter(coding.maf, coord %in% A.zero.coords)

    B.rows = which(coding.maf$time_point==2)
    B.zero.coords = pull(unique(coding.maf[which(coding.maf$time_point==2 & VAF == 0), "coord"]))
    B.zero = dplyr::filter(coding.maf, coord %in% B.zero.coords)

    C.rows = which(coding.maf$time_point==3)
    C.zero.coords = pull(unique(coding.maf[which(coding.maf$time_point==3 & VAF == 0), "coord"]))
    C.zero = dplyr::filter(coding.maf, coord %in% C.zero.coords)

    coding.maf$category = "trunk"
    coding.maf[which(coord %in% A.zero.coords), "category"] = "not-A"
    coding.maf[which(coord %in% B.zero.coords), "category"] = "not-B"
    coding.maf[which(coord %in% C.zero.coords), "category"] = "not-C"

    just_gained_lg_all = dplyr::filter(coding.maf, Hugo_Symbol %in% lg & category != "trunk" )

    just_gained_lg = dplyr::filter(coding.maf, Hugo_Symbol %in% lg & category != "trunk" & time_point ==3) %>%
      mutate(time_point = time_point +0.4)

    just_trunk = dplyr::filter(coding.maf, Hugo_Symbol %in% lg & category == "trunk" & time_point ==1) %>%
      mutate(time_point = time_point-0.4)

    ggplot(coding.maf, aes(x = time_point, y = VAF, group = coord, colour = category)) +
      geom_point() +
      geom_line(alpha = 0.3) +
      geom_text_repel(data = just_gained_lg, aes(label = Hugo_Symbol), size = 4, segment.linetype = 0) +
      geom_text_repel(data = just_trunk, aes(label = Hugo_Symbol), size = 4, segment.linetype = 0) +
      ggtitle(title) +
      theme_minimal()
    return(TRUE)
  }
  return(FALSE)
}


#' Use GISTIC2.0 scores output to reproduce maftools::chromoplot with more flexibility.
#'
#' @param scores output file scores.gistic from the run of GISTIC2.0
#' @param genes_to_label optional. Provide a data frame of genes to label (if mutated). The first 3 columns must contain chromosome, start, and end coordinates. Another required column must contain gene names and be named `gene`. All other columns are ignored. If no data frame provided, oncogenes from GAMBLR packages are used by default to annotate on the plot.
#' @param cutoff optional. Used to determine which regions to color as aberrant. Must be float in the range [0-1]. The higher the number, the less regions will be considered as aberrant. The default is 0.5.
#' @param adjust_amps optional. The value of G-score for highest amplification peak will be multiplied by this value to determine how far up the gene label will be displayed. Default 0.5.
#' @param adjust_dels optional. The value of G-score for highest deletion peak will be multiplied by this value to determine how far down the gene label will be displayed. Default 2.75.
#' @param label_size optional. The font size for the gene label to be displayed. Default 3.
#' @param force_pull optional. How strong the gene name label will be pulled towards a data point. Default 0 (no pulling).
#' @param segment.curvature optional. Indicates whether arrow to the data point should be curved. Accepts numeric value, where negative is for left-hand and positive for right-hand curves, and 0 for straight lines. Default 0.25
#' @param segment.ncp optional. Indicates number of control points to make a smoother curve. Higher value allows for more flexibility for the curve. Default 4
#' @param segment.angle optional. Numeric value in the range 0-180, where less than 90 skews control points of the arrow from label to data point toward the start point. Default 25
#'
#' @return Nothing
#' @export
#' @import tidyverse ggrepel
#'
#' @examples
#' #basic usage
#' prettyChromoplot("path_to_gistic_results/scores.gistic")
#' #advanced usages
#' prettyChromoplot("path_to_gistic_results/scores.gistic", genes_to_label="path_to_gene_coordinates_table.tsv", cutoff=0.75) +
#' ...(any ggplot options to customize plot appearance)
#'
prettyChromoplot = function(scores,
                            genes_to_label,
                            cutoff = 0.5,
                            adjust_amps = 0.5,
                            adjust_dels = 2.75,
                            label_size = 3,
                            force_pull = 0,
                            segment.curvature = 0.25,
                            segment.ncp = 4,
                            segment.angle = 25){

  #read GISTIC scores file, convert G-score to be negative for deletions, and relocate chromosome, start, and end columns to be the first three
  scores = data.table::fread(scores) %>%
    dplyr::mutate(`G-score` = ifelse(Type == "Amp", `G-score`, - 1 * `G-score`)) %>%
    dplyr::relocate(Type, .after = frequency)

  #annotate each region with direction of changes - used for coloring
  scores$fill = ifelse(scores$Type == "Amp" & scores$`-log10(q-value)` > cutoff, "up",
                ifelse(scores$Type == "Del" & scores$`-log10(q-value)` > cutoff, "down", "neutral"))

  #colors to plot
  cnv_palette = c("up" = "#bd0000", "down" = "#2e5096", "neutral" = "#D2D2D3")

  #if no file is provided, annotate with oncogenes in GAMBLR package
  if(missing(genes_to_label)){
    genes_to_label = GAMBLR::grch37_oncogene %>%
      dplyr::mutate(across(c(chrom, start, end), as.integer)) %>%
      data.table::as.data.table()
  }else{
    genes_to_label = data.table::fread(genes_to_label)
    colnames(genes_to_label)[1:3] = c("chrom", "start", "end")
    genes_to_label = genes_to_label %>%

      #for now, drop the X chromosome since GISTIC runs without sex chromosmes
      dplyr::filter(!grepl("X", chrom)) %>%
      dplyr::mutate(across(c(chrom, start, end), as.integer)) %>%
      data.table::as.data.table()
  }
  #overlap scores with genes to annotate
  scores = data.table::foverlaps(scores %>%
    data.table::setkey(., Chromosome, Start, End), genes_to_label %>%
    data.table::setkey(., chrom, start, end), by.x = c("Chromosome", "Start", "End"), by.y = c("chrom", "start", "end"), type = "within") %>%

    #if gene to annotate is provided, but it is in region with no CNV, do not label it
    dplyr::mutate(gene=ifelse(!is.na(gene) & fill=="neutral", NA, gene)) %>%
    #if gene is covering multiple adjacent regions, label only once
    dplyr::group_by(gene) %>%
    dplyr::mutate(newcol = ifelse(!is.na(gene) & !duplicated(gene), gene, NA), gene = newcol) %>%
    dplyr::select(-newcol)

  #get coordinates to label chromosome numbers
  xses = scores %>%
    dplyr::group_by(Chromosome) %>%
    dplyr::mutate(End = max(End) / 2) %>%
    dplyr::pull(End)

  #main plotting
  ggplot(data = scores, aes(x = Start, y = `G-score`, color = fill, label = gene)) +
    geom_bar(size = 0.2, stat = 'identity', position = "dodge") +
    ylab('G-score') +
    ggrepel::geom_text_repel(data = subset(scores, !is.na(gene) & Type == "Amp"),
                             nudge_y = max(subset(scores, !is.na(gene) & Type == "Amp")$`G-score`) * adjust_amps,
                             size = label_size,
                             segment.size = 0.5,
                             segment.color = "#000000",
                             force_pull = force_pull,
                             arrow = arrow(length = unit(0.05, "inches"), type = "closed"),
                             segment.curvature = segment.curvature,
                             segment.ncp = segment.ncp,
                             segment.angle = segment.angle) +
    ggrepel::geom_text_repel(data = subset(scores, !is.na(gene) & Type == "Del"),
                             nudge_y = min(subset(scores, !is.na(gene) & Type == "Del")$`G-score`) * adjust_dels,
                             nudge_x = subset(scores, !is.na(gene) & Type == "Del")$Start,
                             size = label_size,
                             segment.size = 0.5,
                             segment.color = "#000000",
                             force_pull = force_pull,
                             arrow = arrow(length = unit(0.05, "inches"), type = "closed"),
                             segment.curvature = segment.curvature,
                             segment.ncp = segment.ncp,
                             segment.angle = segment.angle) +
    facet_grid(. ~ Chromosome, scales = "free") +
    scale_color_manual(values = cnv_palette) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 16, colour = "black"),
          axis.ticks.x = element_blank(), axis.ticks.y = element_line(colour = "black"), legend.position = "none",
          panel.spacing.x = unit(0.1, "lines"), panel.border = element_blank(), text = element_text(size = 16, colour = "black"),
          strip.background = element_blank(), strip.text.x = element_blank(), panel.grid = element_blank()) +
    geom_hline(yintercept = 0, size = 7) +
    geom_text(aes(label = Chromosome, x = xses, y = 0), size = 4, color = "white")
}


#' Define function for consistent plot theme.
#'
#' @param base_size Size of the font on the plot. Defaults to 14.
#' @param base_family Font family to be used on the plot. Defaults to Arial. Always use cairo device when saving the resulting plot!
#' @param my_legend_position Where to draw the legend? Defaults to the bottom of the plot.
#' @param my_legend_direction Which direction to draw the legend? Defaults to horizontal.
#'
#'
#' @return Nothing.
#' @export
#' @import ggplot2 ggthemes
#'
#' @examples
#' ggplot(mpg, aes(displ, hwy, colour = class)) +
#' geom_point() +
#' theme_Morons()
#'
theme_Morons = function(base_size = 14,
                        base_family = "Arial",
                        my_legend_position = "bottom",
                        my_legend_direction = "horizontal"){

  library(ggthemes)
  (theme_foundation(base_size = base_size, base_family = base_family) +
   theme(plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
         text = element_text(colour = "black"),
         panel.background = element_rect(colour = NA),
         plot.background = element_rect(colour = NA),
         panel.border = element_rect(colour = NA),
         axis.title = element_text(face = "bold", size = rel(1.2)),
         axis.title.y = element_text(angle = 90, vjust = 2),
         axis.title.x = element_text(vjust = -0.2),
         axis.text = element_text(size = base_size, family = base_family),
         axis.line = element_line(colour = "black", size = rel(0.8)),
         axis.ticks = element_line(),
         panel.grid.major = element_line(colour = "#f0f0f0"),
         panel.grid.minor = element_blank(),
         legend.key = element_rect(colour = NA),
         legend.position = my_legend_position,
         legend.direction = my_legend_direction,
         legend.title = element_text(face = "italic"),
         strip.background = element_rect(color = "black", fill = "white", size = 1, linetype = "solid"),
         strip.text = element_text(face = "bold")
    ))
}


#' Create a forest plot comparing mutation frequencies for a set of genes between two groups.
#'
#' @param maf A maf data frame. Minimum required columns are Hugo_Symbol and Tumor_Sample_Barcode.
#' @param mutmat Optional argument for binary mutation matrix. If not supplied, function will generate this matrix from the file used in argument "maf".
#' @param metadata Metadata for the comparisons. Minimum required columns are Tumor_Sample_Barcode and the column assigning each case to one of two groups.
#' @param genes An optional list of genes to restrict your plot to. If no gene-list is supplied, the function will extract all mutated genes from the incoming maf. See min_mutations parameter for more info.
#' @param min_mutations Optional parameter for when no gene list is provided. This parameter ensures only genes with n mutations are kept in the gene list. Default value is 1, this means all genes in the incoming maf will be plotted.
#' @param comparison_column Mandatory: the name of the metadata column containing the comparison values.
#' @param rm_na_samples Set to TRUE to remove 0 mutation samples. Default is FALSE.
#' @param comparison_values Optional: If the comparison column contains more than two values or is not a factor, specify a character vector of length two in the order you would like the factor levels to be set, reference group first.
#' @param separate_hotspots Optional: If you would like to treat hotspots separately from other mutations in any gene. Requires that the maf file is annotated with GAMBLR::annotate_hotspots.
#' @param comparison_name Optional: Specify the legend title if different from the comparison column name.
#' @param custom_colours Optional: Specify a named vector of colours that match the values in the comparison column.
#' @param custom_labels Optional: Specify custom labels for the legend categories. Must be in the same order as comparison_values.
#' @param max_q cut off for q values to be filtered in fish test
#'
#' @return A convenient list containing all the data frames that were created in making the plot, including the mutation matrix.
#' @return It also produces (and returns) ggplot object with a side-by-side forest plot and bar plot showing mutation incidences across two groups.
#' @export
#' @import dplyr cowplot broom
#'
#' @examples
#' metadata = get_gambl_metadata(case_set = "tFL-study") #%>%
#'   dplyr::filter(pairing_status == "matched") %>%
#'   dplyr::filter(consensus_pathology %in% c("FL", "DLBCL"))
#'
#' maf = get_coding_ssm(limit_samples = metadata$sample_id, basic_columns = TRUE, seq_type = "genome")
#'
#' prettyForestPlot(maf = maf,
#'                  metadata = metadata,
#'                  genes = c("ATP6V1B2", "EZH2", "TNFRSF14", "RRAGC"),
#'                  comparison_column = "consensus_pathology",
#'                  comparison_values = c("DLBCL", "FL"),
#'                  separate_hotspots = FALSE,
#'                  comparison_name = "FL vs DLBCL")
#'
prettyForestPlot = function(maf,
                            mutmat,
                            metadata,
                            genes,
                            min_mutations = 1,
                            comparison_column,
                            rm_na_samples = FALSE,
                            comparison_values = FALSE,
                            separate_hotspots = FALSE,
                            comparison_name = FALSE,
                            custom_colours = FALSE,
                            custom_labels = FALSE,
                            max_q = 1){

  #If no comparison_values are specified, derive the comparison_values from the specified comparison_column
  if(comparison_values[1] == FALSE){
    if(class(metadata[[comparison_column]]) == "factor"){
      comparison_values = levels(metadata[[comparison_column]])
    } else {
      comparison_values = unique(metadata[[comparison_column]])
    }
  }
  #Ensure there are only two comparison_values
  {
    if(length(comparison_values) != 2)
      stop(paste0("Your comparison must have two values. \nEither specify comparison_values as a vector of length 2 or subset your metadata so your comparison_column has only two unique values or factor levels."))
  }
  #Subset the metadata to the specified comparison_values and the maf to the remaining sample_ids
  metadata = metadata[metadata[[comparison_column]] %in% comparison_values, ]

  #Ensure the metadata comparison column is a factor with levels matching the input
  metadata$comparison = factor(metadata[[comparison_column]], levels = comparison_values)

  #read maf into r
  if(!missing(maf)){
    #extract gene symbols from maf with minimum N mutations (if no genes list is provided)
    if(missing(genes)){
      genes = maf %>%
        dplyr::select(Hugo_Symbol) %>%
        add_count(Hugo_Symbol) %>%
        distinct(Hugo_Symbol, .keep_all = TRUE) %>%
        dplyr::filter(n >= min_mutations) %>%
        pull(Hugo_Symbol)
    }
    maf = maf[maf$Hugo_Symbol %in% genes, ]
    maf = maf[maf$Tumor_Sample_Barcode %in% metadata$Tumor_Sample_Barcode, ]
  }

  #If separate_hotspots = true, confirm the input maf is hotspot annotated
  if(!missing(mutmat)){
    #add the required columns from the metadata and make the names consistent
    mutmat = left_join(dplyr::select(metadata, sample_id, comparison),mutmat) %>%
      dplyr::rename("Tumor_Sample_Barcode"="sample_id")
  }else if(!missing(maf)){
    if(separate_hotspots){
      if(!"hot_spot" %in% colnames(maf))
        stop("No \"hot_spot\" column in maf file. Annotate your maf file with GAMBLR::annotate_hot_spots() first. ")
      maf$Hugo_Symbol = ifelse(!is.na(maf$hot_spot), paste0(maf$Hugo_Symbol, "_hotspot"), maf$Hugo_Symbol)
    }
    #Convert the maf file to a binary matrix
    mutmat = maf %>%
      dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode) %>%
      full_join(dplyr::select(metadata, Tumor_Sample_Barcode, comparison), by = "Tumor_Sample_Barcode") %>%
      distinct() %>%
      dplyr::mutate(is_mutated = 1) %>%
      pivot_wider(names_from = Hugo_Symbol, values_from = is_mutated, values_fill = 0) %>%
      dplyr::mutate(across(where(is.numeric), ~replace_na(., 0)))
    if(rm_na_samples && "NA" %in% colnames(mutmat)){
      mutmat = mutmat %>%
        dplyr::filter(`NA` == 0) #filter out all samples that show no mutations in the selected genes (i.e na_samples = 1).
    }
    if("NA" %in% colnames(mutmat)){
      mutmat = dplyr::select(mutmat, -`NA`) #remove NA column (if there).
    }
  }else{
    message("provide a MAF or mutation matrix")
    return()
  }

  fish_test = mutmat %>%
    pivot_longer(-c(Tumor_Sample_Barcode, comparison), names_to = "gene", values_to = "is_mutated") %>%
    dplyr::mutate(is_mutated = factor(is_mutated, levels = c("1", "0"))) %>%
    group_by(gene) %>%
    dplyr::summarise(table = list(table(is_mutated, comparison))) %>%
    dplyr::mutate(test = map(table, fisher.test), tidy = map(test, broom::tidy)) %>%
    unnest(tidy) %>%
    dplyr::mutate(q.value = p.adjust(p.value, "BH")) %>%
    dplyr::select(-c(test, method, alternative)) %>%
    dplyr::filter(q.value <= max_q) %>%
    dplyr::mutate(gene = fct_reorder(gene, estimate))

  flatten_table <- function(Row){

    mut_n <- Row[2] %>%
    as.data.frame %>%
    `colnames<-`(
      gsub(
        "table.",
        "",
        colnames(.)
      )
    ) %>%
    mutate(
      is_mutated = ifelse(
        is_mutated == 1,
        "mutated",
        "non-mutated"
      ),
      group = paste(
        is_mutated,
        comparison,
        sep = "_"
      )
    ) %>%
    select(group, Freq) %>%
    pivot_wider(
        names_from = group,
        values_from = Freq
    ) %>%
    t %>%
    as.data.frame

    Row <- Row[-2] %>%
      do.call(cbind, .) %>%
      as.data.frame %>%
      t %>% as.data.frame

    rbind(Row, mut_n)

  }


  fish_test <- apply(
    fish_test,
    1,
    flatten_table
  ) %>%
  do.call(cbind, .) %>%
  t %>%
  as.data.frame %>%
  `rownames<-`(NULL) %>%
  mutate_at(c(2:10), as.numeric) %>%
  arrange(estimate)


  point_size = 50 / round(length(fish_test$gene))
  if(point_size < 1){
    point_size = 1
  }
  font_size = 360 / round(length(fish_test$gene))
  if(font_size < 4){
    font_size = 4
  }else if(font_size > 20){
    font_size = 20
  }
  message(paste("FONT:", font_size, "POINT:", point_size, length(fish_test$gene)))
  forest = fish_test %>%
    dplyr::mutate(gene = factor(gene, levels = fish_test$gene)) %>%
    ggplot(aes(x = gene, y = log(estimate))) +
    geom_point(size = point_size, shape = "square") +
    geom_hline(yintercept = 0, lty = 2) +
    coord_flip() +
    geom_errorbar(aes(ymin = log(conf.low), ymax = log(conf.high), width = 0.2)) +
    ylab("ln(Odds Ratio)") +
    xlab("Mutated Genes") +
    cowplot::theme_cowplot() +
    theme(axis.text.y = element_text(size = font_size))

  if(comparison_name == FALSE){
    comparison_name = comparison_column
  }

  if(custom_colours[1] == FALSE){
    if(length(levels(metadata$comparison)[levels(metadata$comparison) %in% names(get_gambl_colours())]) == 2){
      colours = get_gambl_colours()[levels(metadata$comparison)]
    } else {
      colours = get_gambl_colours(classification = "blood")[c("Red", "Blue")]
      names(colours) = levels(metadata$comparison)
    }
  } else {
    colours = custom_colours
  }

  if(custom_labels[1] == FALSE){
    labels = levels(metadata$comparison)
    names(labels) = levels(metadata$comparison)
  } else if(length(custom_labels) != 2) {
    labels = levels(metadata$comparison)
    names(labels) = levels(metadata$comparison)
    print("Provided custom labels is not a character vector of length 2. Defaulting to comparison factor levels as labels. ")
  } else {
    labels = custom_labels
    names(labels) = comparison_values
  }

  bar = mutmat %>%
    dplyr::select(-Tumor_Sample_Barcode) %>%
    pivot_longer(
		  !comparison,
		  names_to = "gene",
		  values_to = "is_mutated"
	  ) %>%
    group_by(gene, comparison) %>%
    drop_na() %>%
    summarise(percent_mutated = sum(is_mutated) / n() * 100) %>%
    dplyr::filter(gene %in% fish_test$gene) %>%
    dplyr::mutate(gene = factor(gene, levels = fish_test$gene)) %>%
    ggplot(aes(x = gene, y = percent_mutated, fill = comparison)) +
    geom_col(position = "dodge", width = 0.5) +
    xlab("") + ylab("% Mutated") +
    coord_flip() +
    scale_fill_manual(name = comparison_name, values = colours, labels = labels[levels(metadata$comparison)]) +
    cowplot::theme_cowplot() +
    theme(axis.text.y = element_blank(), legend.position = "bottom")

  legend = cowplot::get_legend(bar)

  plots = plot_grid(forest, bar +
                      theme(legend.position = "none"), rel_widths = c(1, 0.6), nrow = 1)

  arranged_plot = cowplot::plot_grid(plot_grid(NULL, legend, NULL, nrow = 1), plots, nrow = 2, rel_heights = c(0.1, 1))

  return(list(fisher = fish_test, forest = forest, bar = bar, legend = legend, arranged = arranged_plot, mutmat = mutmat))
}


#' Make an heatmap that is looking cute using ComplexHeatmap. The metadata is expected to follow the structure and column naming used in GAMBL.
#' If you provide your own non-GAMBL samples and metadata, you must include at least the columns with names corresponding to annotation tracks and column "Tumor_Sample_Barcode".
#' showing sample ids. The metadata can contain numeric columns, which will be plotted as numeric variables in the annotation. The efature matrix is supplied in this_matrix argument.
#' and is expected to have samples in rows, and features in columns. The argument importance_values is similar to the widths of NMF object or importance values for feature/group from RF models.
#' It is also expected to have column names (having names of the groups that will be shown on heatmap) and rownames (corresponding to feature ids).
#' @param this_matrix A data frame with column Tumor_Sample_Barcode and a column for each feature. Can be binary. Expected to not contain negative values.
#' @param importance_values Provide a data frame of feature (in rows) by group (in columns) with numeric values representative of feature importance. Can be obtained from rf$inportance or basis(NMF).
#' @param these_samples_metadata Data frame containing metadata for your samples.
#' @param max_number_of_features_per_group Optional argument to indicate how many features from each group to be considered for display. Default is 10.
#' @param splitColumnName Optional argument to indicate which metadata column to split on. Default is set to pathology.
#' @param metadataColumns A vector containing the categorical column names you want to plot below.
#' @param numericMetadataColumns A vector containing the numeric columns you want to plot below.
#' @param numericMetadataMax A numeric vector of cutoffs to apply to numeric columns above.
#' @param custom_colours Provide named vector (or named list of vectors) containing custom annotation colours if you do not want to use standartized pallette.
#' @param prioritize_ordering_on_numeric Logical argument specifying whether to sort on numeric metadata first or other metadata columns. Default is TRUE (sort on numeric metadata, then on other columns).
#' @param legend_direction Optional argument to indicate whether legend should be in horizontal (default) or vertical position.
#' @param legend_position Optional argument to indicate where the legend should be drawn. The default is set to bottom, but can also accept top, right, and left.
#' @param legend_row Fiddle with these to widen or narrow your legend (default 3).
#' @param legend_col Fiddle with these to widen or narrow your legend (default 3).
#' @param fontSizeGene Font size for gene labels (default 6).
#' @param metadataBarHeight Optional argument to adjust the height of bar with annotations. The default is 1.5
#' @param leftStackedWidth Optional argument to control how wide should the stacked plot on the left be. The default is 4.
#' @param metadataBarFontsize Optional argument to control for the font size of metadata annotations. The default is 5.
#' @param groupNames optional vector of group names to be displayed above heatmap. Should be the same length as the number of groups that will be shown. Default is NULL (no labels).
#'
#' @return Nothing
#' @export
#' @import ComplexHeatmap grid dplyr circlize
#'
#' @examples
#' splendidHeatmap(
#'  this_matrix = data,
#'  importance_values = rf$importance[,c(1:3)],
#'  these_samples_metadata = MASTER.METADATA,
#'  splitColumnName = "pathology",
#'  metadataColumns = c("cohort", "pathology", "sex", ".", "COO_consensus", "DHITsig_consensus", "seq_type"),
#'  numericMetadataColumns = ".",
#'  numericMetadataMax = 0.7,
#'  custom_colours=custom_colours)
#'
splendidHeatmap = function(this_matrix,
                           importance_values,
                           these_samples_metadata,
                           max_number_of_features_per_group = 10,
                           splitColumnName = "pathology",
                           metadataColumns = c("pathology"),
                           numericMetadataColumns = NULL,
                           numericMetadataMax = NULL,
                           prioritize_ordering_on_numeric = TRUE,
                           custom_colours = NULL,
                           legend_direction = "horizontal",
                           legend_position = "bottom",
                           legend_row = 3,
                           legend_col = 3,
                           fontSizeGene = 6,
                           metadataBarHeight = 1.5,
                           leftStackedWidth = 4,
                           metadataBarFontsize = 5,
                           groupNames = NULL){
  comparison_groups = colnames(importance_values)

  if(!is.null(splitColumnName) & (splitColumnName %in% metadataColumns)){
    metadataColumns = c(splitColumnName, metadataColumns[!metadataColumns == splitColumnName])
  }

  if(!is.null(numericMetadataColumns) & length(intersect(numericMetadataColumns, metadataColumns))>0){
    message(paste0("The column(s) ", numericMetadataColumns, " specified both in metadata and numeric metadata. Plotting as numeric values..."))
    metadataColumns = metadataColumns[!metadataColumns %in% numericMetadataColumns]
  }

  #get which group samples belong to
  metadata_df = these_samples_metadata[,c("Tumor_Sample_Barcode", metadataColumns, numericMetadataColumns)] %>%
    as.data.frame() %>%
    column_to_rownames(., "Tumor_Sample_Barcode")

  if(!is.null(numericMetadataMax)){
      max_list = setNames(numericMetadataMax, numericMetadataColumns)
      metadata_df = metadata_df %>%
        dplyr::mutate(across(names(max_list), ~ ifelse(.x > max_list[[cur_column()]], max_list[[cur_column()]], .x)))
  }

  # count N of features for every dsample and add it to metadata
  metadata_df =
  this_matrix %>%
    as.data.frame %>%
    column_to_rownames("Tumor_Sample_Barcode") %>%
    rowSums %>%
    as.data.frame %>%
    `names<-`("N_features") %>%
    rownames_to_column ("Tumor_Sample_Barcode") %>%
    base::merge(metadata_df %>%
                  rownames_to_column ("Tumor_Sample_Barcode"),
                .) %>%
    column_to_rownames("Tumor_Sample_Barcode")

  my_colours = NULL
  these_names = NULL
  for (i in 1:length(metadataColumns)){
    this_metadata_column = get_gambl_colours(metadataColumns[i])
    if (sum(is.na(names(this_metadata_column[unlist(c(unique(these_samples_metadata[,metadataColumns[i]])))]))) <= 1 &
        nrow(unique(these_samples_metadata[,metadataColumns[i]])) > 1){
      these_names = c(these_names, metadataColumns[i])
      my_colours = append(my_colours, list(c(this_metadata_column, "NA" = "#BDBDC1FF")))
      names(my_colours) = these_names
    }
  }

  my_colours = c(custom_colours, my_colours)

  col_fun=circlize::colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  for(exp in numericMetadataColumns){
    my_colours[[exp]] = col_fun
  }

  #get all features
  w = importance_values[,comparison_groups]
  w = as.data.frame(w) %>%
    dplyr::mutate_if(is.character,as.numeric)

  # extract most important features, while taking the feature with highest weight for a particular cluster if it was seen before for other cluster with lower weight
  FEATURES <- w[,1] %>%
    as.data.frame() %>%
    `rownames<-`(rownames(w)) %>%
    `names<-`("importance") %>%
    dplyr::arrange(desc(importance)) %>%
    head(., max_number_of_features_per_group) %>%
    rownames_to_column(., var = "Feature") %>%
    dplyr::mutate(group = comparison_groups[1])
  for (i in 2:length(comparison_groups)){
    FEATURES = rbind(as.data.frame(FEATURES), w[,i] %>%
      as.data.frame() %>%
      `rownames<-`(rownames(w)) %>%
      `names<-`("importance") %>%
       arrange(desc(importance)) %>%
       head(., max_number_of_features_per_group + 3) %>%
       rownames_to_column(., var = "Feature") %>%
       dplyr::mutate(group = comparison_groups[i])) %>%
       dplyr::mutate(importance=as.numeric(importance)) %>%
       dplyr::group_by(Feature) %>%
       dplyr::filter(importance == max(importance)) %>%
       dplyr::arrange(group)
  }
  FEATURES = as.data.frame(FEATURES)

  mat = this_matrix %>%
    base::merge(., metadata_df %>%
    rownames_to_column(., "Tumor_Sample_Barcode") %>%
    dplyr::select(Tumor_Sample_Barcode, splitColumnName)) %>%
    as.data.frame()
  mat[,splitColumnName] = factor(mat[,splitColumnName])

  #breaks used to display groups with different colors on heatmap
  bk = c(0, seq(0.5, length(comparison_groups) + 0.5, 1))

  #colors used to show on the heatmap body. Starts with white - the color of feature absence
  my_palette = c("white", rev(unlist(my_colours[splitColumnName])))
  my_palette = unname(my_palette)

  #get each group and label the events for each feature with group number
  mat_2 = mat[,-ncol(mat)]
  #subset samples of each group
  MY.LIST = list()
  for (i in 1:length(comparison_groups)){
    MY.LIST[[i]] = assign(comparison_groups[i], mat_2 %>%
      as.data.frame(.) %>%
      column_to_rownames(., var = "Tumor_Sample_Barcode") %>%
      t(.) %>%
      as.data.frame(.) %>%
      dplyr::select(metadata_df %>%
      dplyr::filter(base::get(splitColumnName) == comparison_groups[i]) %>%
      rownames))
  }

  #assign numbers - used for coloring of heatmap body
  for(i in 1:length(comparison_groups)){
    MY.LIST[[i]][MY.LIST[[i]] > 0] = i
  }

  #bind them all together for plotting
  mat_2 = do.call(cbind, MY.LIST) %>%
    as.data.frame(.) %>%
    t(.) %>%
    as.data.frame(.) %>%
    rownames_to_column(., var = "Tumor_Sample_Barcode") %>%
    base::merge(., metadata_df %>%
    rownames_to_column(., "Tumor_Sample_Barcode") %>%
    dplyr::select(Tumor_Sample_Barcode, splitColumnName), by = "Tumor_Sample_Barcode")

  #specify where row breaks should be on heatmap
  FEATURES = FEATURES %>%
    arrange(match(
      group,
      str_sort(FEATURES$group, numeric = TRUE)
    ))
  breaks = 0
  for (this_group in comparison_groups){
    N = (nrow(FEATURES %>%
      dplyr::filter(group == this_group)))
    breaks = c(breaks, N)
  }

  #second, make a vector that will be supplied to ComplexHeatmap
  my_vector = NULL
  for (i in 1:(length(breaks))){
    my_vector = c(my_vector, rep(i - 1, breaks[i]))
  }

  #prepare matrix for stacked barplots on the left
  STACKED = data.frame(matrix(NA, ncol = 1, nrow = nrow(FEATURES)))[-1]
  rownames(STACKED) = FEATURES$Feature
  for (i in 1:length(comparison_groups)) {
  STACKED = cbind(STACKED, mat_2[,c("Tumor_Sample_Barcode", FEATURES$Feature)] %>%
    base::merge(., metadata_df %>%
    rownames_to_column(., "Tumor_Sample_Barcode") %>%
    dplyr::select(Tumor_Sample_Barcode, splitColumnName), by = "Tumor_Sample_Barcode") %>%
    dplyr::arrange(!!sym(splitColumnName)) %>%
    dplyr::filter(base::get(splitColumnName) == comparison_groups[i]) %>%
    dplyr::select(-Tumor_Sample_Barcode, -splitColumnName) %>%
    dplyr::summarise_all(funs(sum)) %>%
    t(.) %>%
    `colnames<-`(comparison_groups[i]) %>%
    as.data.frame(.) %>%
    dplyr::mutate_all(~(./i) / nrow(metadata_df)))
  }

  m = t(apply(STACKED, 1, function(x) x/sum(x)))
  m[is.na(m)] <- 0

  if(prioritize_ordering_on_numeric & ! is.null(numericMetadataColumns)){ # numeric metadata is provided and is prioritized for column sorting
    used_for_ordering_df = t(base::merge(mat_2 %>%
    dplyr::select(-splitColumnName), metadata_df %>%
    rownames_to_column(., "Tumor_Sample_Barcode"), by = "Tumor_Sample_Barcode") %>%
    column_to_rownames(., var = "Tumor_Sample_Barcode") %>%
      dplyr::arrange(desc(!!!syms(numericMetadataColumns)),
        !!!syms(metadataColumns)) %>%
      dplyr::select(FEATURES$Feature))

    this_is_ordered_df = metadata_df[ (order(match(rownames(metadata_df), colnames(used_for_ordering_df)))), ] %>%
      dplyr::arrange(desc(!!!syms(numericMetadataColumns)),
        !!!syms(metadataColumns))
  }else if(! is.null(numericMetadataColumns)){ # numeric metadata is provided, but is not prioritized for column sorting
    used_for_ordering_df = t(base::merge(mat_2 %>%
    dplyr::select(-splitColumnName), metadata_df %>%
    rownames_to_column(., "Tumor_Sample_Barcode"), by = "Tumor_Sample_Barcode") %>%
    column_to_rownames(., var = "Tumor_Sample_Barcode") %>%
      dplyr::arrange(!!!syms(metadataColumns),
        desc(!!!syms(numericMetadataColumns))) %>%
      dplyr::select(FEATURES$Feature))

    this_is_ordered_df = metadata_df[ (order(match(rownames(metadata_df), colnames(used_for_ordering_df)))), ] %>%
      dplyr::arrange(!!!syms(metadataColumns),
        desc(!!!syms(numericMetadataColumns)))
  }else{ # no numeric metadata is proveded to plot
    used_for_ordering_df = t(base::merge(mat_2 %>%
    dplyr::select(-splitColumnName), metadata_df %>%
    rownames_to_column(., "Tumor_Sample_Barcode"), by = "Tumor_Sample_Barcode") %>%
    column_to_rownames(., var = "Tumor_Sample_Barcode") %>%
      dplyr::arrange(!!!syms(metadataColumns)) %>%
      dplyr::select(FEATURES$Feature))

    this_is_ordered_df = metadata_df[ (order(match(rownames(metadata_df), colnames(used_for_ordering_df)))), ] %>%
      dplyr::arrange(!!!syms(metadataColumns))
  }

  # left annotation: stacked feature weights
  ha = rowAnnotation(`feature abundance` = anno_barplot(m, gp = gpar(fill = my_palette[1:length(comparison_groups)+1]),
                                                      bar_width = 1, width = unit(leftStackedWidth, "cm"),
                                                      axis_param = list(side = legend_position, labels_rot = 0)))

  #bottom annotation: tracks indicating metadata
  ha_bottom = HeatmapAnnotation(
    df = this_is_ordered_df %>% dplyr::select(-c(splitColumnName, N_features)),
    col = my_colours,
    simple_anno_size = unit(metadataBarHeight, "mm"),
    gap = unit(0.25 * metadataBarHeight, "mm"),
    annotation_name_gp = gpar(fontsize = metadataBarFontsize),
    annotation_legend_param = list(nrow = legend_row, ncol = legend_col, direction = legend_direction)
  )

  #top annotation: groups of interest to split on
  top_bar_colors = list(my_colours[[splitColumnName]] %>% rev)
  names(top_bar_colors) = splitColumnName
  names(top_bar_colors[[splitColumnName]]) = names(top_bar_colors[[splitColumnName]]) %>% rev()

  ha_top = HeatmapAnnotation(
    group = anno_block(gp = gpar(fill = top_bar_colors[[1]], fontsize = fontSizeGene * 1.5), labels = groupNames),
    "N of features" = anno_barplot(this_is_ordered_df$N_features)
  )

  splendidHM = ComplexHeatmap::Heatmap(used_for_ordering_df,
                                       col = my_palette,
                                       show_column_names = FALSE,
                                       cluster_columns = FALSE,
                                       cluster_rows = FALSE,
                                       row_names_gp = gpar(fontsize = fontSizeGene),
                                       show_heatmap_legend = FALSE,
                                       row_split = my_vector,
                                       row_title = NULL,
                                       left_annotation = ha,
                                       bottom_annotation = ha_bottom,
                                       top_annotation = ha_top,
                                       column_split = dplyr::pull(metadata_df[(order(match(rownames(metadata_df), colnames(used_for_ordering_df)))), ], splitColumnName),
                                       column_title = NULL)

  draw(splendidHM, heatmap_legend_side = legend_position, annotation_legend_side = legend_position)
}


#' Visualizing variant (SSM or SVs) counts per chromosome
#'
#' @param this_sample Sample to be plotted.
#' @param maf_data Optional parameter with maf like df already loaded into R.
#' @param maf_path Optional parameter with path to external maf like file.
#' @param ssm Set to FALSE to get plotting data from get_combined_sv (SVs). Default value is TRUE (plots SSM retrieved from annotate_cn_by_ssm$maf)
#' @param projection Genome build for returned variants (only applicable for ssm = FALSE)
#' @param min_vaf The minimum tumour VAF for a SV to be returned. Recommended: 0 (only applicable for ssm = FALSE).
#' @param variant_type_col Index of column holding Variant Type (to be used with either maf_data or maf_path).
#' @param chromosome_col Index of column holding Chromosome (to be used with either maf_data or maf_path).
#' @param plot_title Title of plot (default to sample ID).
#' @param y_interval Optional parameter for specifying intervals on y-axis.
#' @param hide_legend Set to True to remove legend from plot, default is FALSE.
#' @param plot_subtitle Subtitle for created plot.
#' @param chr_select vector of chromosomes to be included in plot, defaults to autosomes.
#' @param coding_only Optional. Set to TRUE to restrict to plotting only coding mutations.
#' @param from_flatfile If set to true the function will use flat files instead of the database.
#' @param use_augmented Boolean statement if to use augmented maf, default is FALSE.
#' @param add_qc_metric Boolean statement, if set to TRUE specified QC metric will be added (second y-axis).
#' @param seq_type Default is "genome".
#'
#' @return plot as ggplot object.
#' @import tidyverse cowplot
#' @export
#'
#' @examples
#' ssm = fancy_v_chrcount(this_sample = "HTMCP-01-06-00422-01A-01D", ssm = TRUE)
#' svs = fancy_v_chrcount(this_sample = "HTMCP-01-06-00422-01A-01D", ssm = FALSE,
#'                     min_vaf = 0,
#'                     projection = "grch37",
#'                     chr_select = paste0("chr", c(1:5)),
#'                     plot_subtitle = "SV Count Distribution (chr1-5)")
#'
fancy_v_chrcount = function(this_sample,
                            maf_data,
                            maf_path = NULL,
                            ssm = TRUE,
                            projection = "grch37",
                            min_vaf = 0,
                            variant_type_col = 10,
                            chromosome_col = 5,
                            plot_title = paste0(this_sample),
                            y_interval = 1,
                            hide_legend = FALSE,
                            plot_subtitle = "Variant Count Distribution Per Chromosome",
                            chr_select = paste0("chr", c(1:22)),
                            coding_only = FALSE,
                            from_flatfile = TRUE,
                            use_augmented_maf = TRUE,
                            add_qc_metric = FALSE,
                            seq_type = "genome"){

  if(!missing(maf_data)){
    maf = maf_data
    maf = as.data.frame(maf)
    colnames(maf)[variant_type_col] = "Variant_Type"
    colnames(maf)[chromosome_col] = "Chromosome"

  }else if(!is.null(maf_path)){
    maf = fread_maf(maf_path)
    maf = as.data.frame(maf)
    colnames(maf)[variant_type_col] = "Variant_Type"
    colnames(maf)[chromosome_col] = "Chromosome"
  }

  #get maf data for a specific sample.
  if(missing(maf_data) && is.null(maf_path)){
    if(ssm){
      maf = assign_cn_to_ssm(this_sample = this_sample, coding_only = coding_only, from_flatfile = from_flatfile, use_augmented_maf = use_augmented_maf, this_seq_type = seq_type)$maf
    }else{
      maf = get_combined_sv(sample_ids = this_sample, projection = projection, min_vaf = min_vaf) %>%
        dplyr::select(CHROM_A, START_A, END_A, manta_name)

      #get manta results in required format
      maf = data.frame(maf$CHROM_A, maf$START_A, maf$END_A, do.call(rbind, strsplit(maf$manta_name, split = ":", fixed = TRUE)))

      #rename variables
      names(maf)[1:4] = c("Chromosome", "Start_Position", "End_Position","Variant_Type")

      #filter out translocations and set order of variables
      maf = dplyr::filter(maf, Variant_Type %in% c("MantaDEL", "MantaDUP")) %>%
        dplyr::select(Chromosome, Start_Position, End_Position, Variant_Type)

      #remove "Manta" from Variant_Type string
      maf$Variant_Type = gsub("^.{0,5}", "", maf$Variant_Type)
    }
  }

  #convert variables to factors
  maf$Variant_Type = as.factor(maf$Variant_Type)
  maf$Chromosome = as.factor(maf$Chromosome)

  #add chr prefix if missing
  if(!str_detect(maf$Chromosome, "chr")){
    maf = mutate(maf, Chromosome = paste0("chr", Chromosome))
  }

  #subset data frame on sv sub type
  maf_del = dplyr::filter(maf, Variant_Type == "DEL") %>%
    add_count(Chromosome) %>%
    distinct(Chromosome, .keep_all = TRUE) %>%
    dplyr::select(Chromosome, Variant_Type, n)

  if(ssm){
    maf_ins = dplyr::filter(maf, Variant_Type == "INS") %>%
      add_count(Chromosome) %>%
      distinct(Chromosome, .keep_all = TRUE) %>%
      dplyr::select(Chromosome, Variant_Type, n)

    #combine data frames
    maf.count = rbind(maf_del, maf_ins)

    #get max number of mutations for the chromosome harboring most variants (for setting y-axis value).
    ymax = max(maf_del$n) + max(maf_ins$n)

  }else{
    maf_dup = dplyr::filter(maf, Variant_Type == "DUP") %>%
      add_count(Chromosome) %>%
      distinct(Chromosome, .keep_all = TRUE) %>%
      dplyr::select(Chromosome, Variant_Type, n)

    #combine data frames
    maf.count = rbind(maf_del, maf_dup)

    #get max number of mutations for the chromosome harboring most variants (for setting y-axis value).
    ymax = max(maf_del$n) + max(maf_dup$n)
  }

  if(add_qc_metric){
    #get qc data for selected samples
    sample_df = data.frame(sample_id = this_sample)
    qc_metrics = collate_results(sample_table = sample_df, seq_type_filter = seq_type) %>%
      dplyr::select(MeanCorrectedCoverage)
    if(nrow(qc_metrics) < 1){
      message("No QC metrics available for selected sample...")
    }
  }

  #plot
  p = ggplot(maf.count, aes(x = Chromosome, y = n, fill = Variant_Type, label = n)) +
        labs(title = plot_title, subtitle = plot_subtitle, x = "", y = "Variants (n)", fill = "") +
        scale_x_discrete(expand = c(0, 0.58), limits = chr_select) +
        geom_bar(position = "stack", stat = "identity") +
        {if(ssm)scale_fill_manual(values = get_gambl_colours("indels"))} +
        {if(!ssm)scale_fill_manual(values = get_gambl_colours("svs"))} +
        {if(add_qc_metric)geom_hline(qc_metrics, mapping = aes(yintercept = MeanCorrectedCoverage / 10), linetype = "dashed", group = 2)} +
        {if(!add_qc_metric)scale_y_continuous(expand = c(0, 0), breaks = seq(0, ymax + 2, by = y_interval))} +
        {if(add_qc_metric)scale_y_continuous(expand = c(0, 0), breaks = seq(0, ymax + 2, by = y_interval), sec.axis = sec_axis(~.*10, name = "Mean Corrected Coverage (X)", breaks = seq(0, 100, by = 10)))} +
        theme_cowplot() +
        {if(hide_legend)theme(legend.position = "none")} +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(p)
}


#' Generate a plot with SNV distribution per chromosome.
#'
#' @param this_sample Sample to be plotted.
#' @param maf_data Optional parameter with maf like df already loaded into R.
#' @param maf_path Optional parameter with path to external maf like file.
#' @param variant_type_col Index of column holding Variant Type (to be used with either maf_data or maf_path).
#' @param chromosome_col Index of column holding Chromosome (to be used with either maf_data or maf_path).
#' @param plot_title Title of plot (default to sample ID).
#' @param plot_subtitle Subtitle for created plot.
#' @param chr_select vector of chromosomes to be included in plot, defaults to autosomes.
#' @param include_dnp Optional argument for including DNPs. Default is FALSE.
#' @param hide_legend Set to True to remove legend from plot, default is FALSE.
#' @param coding_only Optional. Set to TRUE to restrict to plotting only coding mutations.
#' @param from_flatfile If set to true the function will use flat files instead of the database.
#' @param use_augmented Boolean statement if to use augmented maf, default is FALSE.
#'
#' @return Nothing.
#' @import tidyverse cowplot
#' @export
#'
#' @examples
#' snv_plot = fancy_snv_chrdistplot(this_sample = "HTMCP-01-06-00422-01A-01D")
#' snv_dnp_plot = fancy_snv_chrdistplot(this_sample = "HTMCP-01-06-00422-01A-01D", include_dnp = TRUE, plot_subtitle = "SNV + DNP Distribution Per Chromosome")
#'
fancy_snv_chrdistplot = function(this_sample,
                                 maf_data,
                                 maf_path = NULL,
                                 variant_type_col = 10,
                                 chromosome_col = 5,
                                 plot_title = paste0(this_sample),
                                 plot_subtitle = "SNV Distribution Per Chromosome",
                                 chr_select = paste0("chr", c(1:22)),
                                 include_dnp = FALSE,
                                 hide_legend = FALSE,
                                 coding_only = FALSE,
                                 from_flatfile = TRUE,
                                 use_augmented_maf = TRUE){

  if(!missing(maf_data)){
    maf = maf_data
    maf = as.data.frame(maf)
    colnames(maf)[variant_type_col] = "Variant_Type"
    colnames(maf)[chromosome_col] = "Chromosome"

  }else if (!is.null(maf_path)){
    maf = fread_maf(maf_path)
    maf = as.data.frame(maf)
    colnames(maf)[variant_type_col] = "Variant_Type"
    colnames(maf)[chromosome_col] = "Chromosome"
  }

  #get maf data for a specific sample.
  if(missing(maf_data) && is.null(maf_path)){
    maf = assign_cn_to_ssm(this_sample = this_sample, coding_only = coding_only, from_flatfile = from_flatfile, use_augmented_maf = use_augmented_maf)$maf
  }

  #add chr prefix if missing
  if(!str_detect(maf$Chromosome, "chr")[5]){
    maf = mutate(maf, Chromosome = paste0("chr", Chromosome))
  }

  #subset data frame on snp sub type
  maf_snp = dplyr::filter(maf, Variant_Type == "SNP") %>%
    add_count(Chromosome) %>%
    distinct(Chromosome, .keep_all = TRUE) %>%
    dplyr::select(Chromosome, Variant_Type, n)

  if(!include_dnp){
    #get max number of SNP for the chromosome harboring most variants (for setting y-axis value).
    ymax = max(maf_snp$n)

    #plot
    ggplot(maf_snp, aes(x = Chromosome, y = n)) +
      labs(title = plot_title, subtitle = plot_subtitle, x = "", y = "Count (n)", fill = "") +
      scale_x_discrete(expand = c(0, 0.7), limits = chr_select) +
      geom_bar(position = "stack", stat = "identity", fill = "#2B9971", width = 0.75) +
      scale_y_continuous(expand = c(0, 0)) +
      theme_cowplot() +
      coord_flip()
  }

  else{
    #subset data frame on snp sub type
    maf_dnp = dplyr::filter(maf, Variant_Type == "DNP") %>%
      add_count(Chromosome) %>%
      distinct(Chromosome, .keep_all = TRUE) %>%
      dplyr::select(Chromosome, Variant_Type, n)

    #combine data frames
    maf.count = rbind(maf_snp, maf_dnp)

    #get max number of mutations for the chromosome harboring most variants (for setting y-axis value).
    ymax = max(maf_snp$n) + max(maf_dnp$n)

    #plot
    ggplot(maf.count, aes(x = Chromosome, y = n, fill = Variant_Type)) +
      labs(title = plot_title, subtitle = plot_subtitle, x = "", y = "SNV Count (n)", fill = "") +
      scale_x_discrete(expand = c(0, 0.7), limits = chr_select) +
      geom_bar(position = "stack", stat = "identity", width = 0.75) +
      scale_fill_manual("", values = c("SNP" = "#2B9971", "DNP" = "#993F2B")) +
      scale_y_continuous(expand = c(0, 0)) +
      theme_cowplot() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      {if(hide_legend)theme(legend.position = "none")} +
      coord_flip()
  }
}


#' Generate a bar plot visualizing total variant (SSM or SVs) count for selected contigs.
#'
#' @param this_sample Sample to be plotted.
#' @param maf_data Optional parameter with maf like df already loaded into R.
#' @param maf_path Optional parameter with path to external maf like file.
#' @param ssm Set to FALSE to get plotting data from get_combined_sv (SVs). Default value is TRUE (plots SSM retrieved from annotate_cn_by_ssm$maf).
#' @param projection Genome build for returned variants (only applicable for ssm = FALSE).
#' @param min_vaf The minimum tumour VAF for a SV to be returned. Recommended: 0 (only applicable for ssm = FALSE).
#' @param variant_type_col Index of column holding Variant Type (to be used with either maf_data or maf_path).
#' @param chromosome_col Index of column holding Chromosome (to be used with either maf_data or maf_path).
#' @param plot_title Title of plot (default to sample ID).
#' @param plot_subtitle Subtitle for created plot.
#' @param chr_select vector of chromosomes to be included in plot, defaults to autosomes.
#' @param variant_select Subtypes of SVs to be incldued in plot, default is DEL, INS and DUP.
#' @param snp_colours Optional vector with colours for SNPs (DNP and TNP).
#' @param hide_legend Set to True to remove legend from plot, default is FALSE.
#' @param coding_only Optional. Set to TRUE to restrict to plotting only coding mutations.
#' @param log10_y Set to TRUE to force y axis to be in log10.
#' @param from_flatfile If set to true the function will use flat files instead of the database.
#' @param use_augmented Boolean statement if to use augmented maf, default is TRUE.
#'
#' @return plot as ggplot object.
#' @import tidyverse cowplot
#' @export
#'
#' @examples
#' chr1_sv = fancy_v_count(this_sample = "HTMCP-01-06-00422-01A-01D", chr_select = c(1))
#' svs = fancy_v_count(this_sample = "HTMCP-01-06-00422-01A-01D")
#'
fancy_v_count = function(this_sample,
                         maf_data,
                         maf_path = NULL,
                         ssm = TRUE,
                         projection = "grch37",
                         min_vaf = 0,
                         variant_type_col = 10,
                         chromosome_col = 5,
                         plot_title = paste0(this_sample),
                         plot_subtitle = "Variant Count For Selected Contigs",
                         chr_select = paste0("chr", c(1:22)),
                         variant_select = c("DEL", "INS", "DUP"),
                         snp_colours = c("SNP" = "#2B9971", "DNP" = "#993F2B", "TNP" = "#A62656"),
                         hide_legend = FALSE,
                         coding_only = FALSE,
                         log10_y = FALSE,
                         from_flatfile = TRUE,
                         use_augmented_maf = TRUE){

  if(!missing(maf_data)){
    maf = maf_data
    maf = as.data.frame(maf)
    colnames(maf)[variant_type_col] = "Variant_Type"
    colnames(maf)[chromosome_col] = "Chromosome"

  }else if (!is.null(maf_path)){
    maf = fread_maf(maf_path)
    maf = as.data.frame(maf)
    colnames(maf)[variant_type_col] = "Variant_Type"
    colnames(maf)[chromosome_col] = "Chromosome"
  }

  #get maf data for a specific sample.
  if(missing(maf_data) && is.null(maf_path)){
    if(ssm){
      maf = assign_cn_to_ssm(this_sample = this_sample, coding_only = coding_only, from_flatfile = from_flatfile, use_augmented_maf = use_augmented_maf)$maf
    }else{
      maf = get_combined_sv(sample_ids = this_sample, projection = projection, min_vaf = min_vaf) %>%
        dplyr::select(CHROM_A, START_A, END_A, manta_name)

      #get manta results in required format
      maf = data.frame(maf$CHROM_A, maf$START_A, maf$END_A, do.call(rbind, strsplit(maf$manta_name, split = ":", fixed = TRUE)))

      #rename variables
      names(maf)[1] = "Chromosome"
      names(maf)[2] = "Start_Position"
      names(maf)[3] = "End_Position"
      names(maf)[4] = "Variant_Type"

      #filter out translocations and set order of variables
      maf = dplyr::filter(maf, Variant_Type %in% c("MantaDEL", "MantaDUP")) %>%
        dplyr::select(Chromosome, Start_Position, End_Position, Variant_Type)

      #remove "Manta" from Variant_Type string
      maf$Variant_Type = gsub("^.{0,5}", "", maf$Variant_Type)
    }
  }

  #add chr prefix if missing
  if(!str_detect(maf$Chromosome, "chr")[1]){
    maf = mutate(maf, Chromosome = paste0("chr", Chromosome))
  }

  #read maf into R and select relevant variables and transform to factor.
  maf_df = dplyr::select(maf, Chromosome, Start_Position, End_Position, Variant_Type) %>%
    mutate_at(vars(Chromosome, Variant_Type), list(factor))

  #sub-setting maf based on user-defined parameters
  maf_df = maf_df[maf_df$Chromosome %in% chr_select, ]
  maf_df = maf_df[maf_df$Variant_Type %in% variant_select, ]

  #subset variant data
  sv_count = maf_df %>%
    group_by(Variant_Type) %>%
    summarize(count = n())

  #get colours
  indels_cols = get_gambl_colours("indels")
  colours = append(indels_cols, snp_colours)

  #plot
  p = ggplot(sv_count, aes(x = Variant_Type, y = count, fill = Variant_Type, label = count)) +
    geom_bar(position = "stack", stat = "identity") +
    {if(log10_y)labs(title = plot_title, subtitle = plot_subtitle, x = "", y = "Variants n (log10)", fill = "")} +
    {if(!log10_y)labs(title = plot_title, subtitle = plot_subtitle, x = "", y = "Variants (n)", fill = "")} +
    {if(ssm)scale_fill_manual(values = colours)} +
    {if(!ssm)scale_fill_manual(values = get_gambl_colours("svs"))} +
    geom_text(size = 5, position = position_stack(vjust = 0.5)) +
    {if(log10_y)scale_y_log10(expand = c(0, 0))} +
    {if(!log10_y)scale_y_continuous(expand = c(0, 0))} +
    {if(hide_legend)theme(legend.position = "none")} +
    theme_cowplot()

  return(p)
}


#' Generate a bar plot visualizing sample-specific copy number states and affected bases for each CN segment.
#'
#' @param this_sample Sample to be plotted.
#' @param seq_data Optional parameter with copy number df already loaded into R.
#' @param seq_path Optional parameter with path to external cn file.
#' @param chrom_col Index of column annotating Chromosome (to be used with either maf_data or maf_path).
#' @param start_col Index of column with copy number start coordinates (to be used with either maf_data or maf_path).
#' @param end_col Index of column with copy number end coordinates (to be used with either maf_data or maf_path).
#' @param cn_col Index of column holding copy number information (to be used with either maf_data or maf_path).
#' @param plot_title Title of plot (default to sample ID).
#' @param plot_subtitle Subtitle for created plot.
#' @param chr_select Vector of chromosomes to be included in plot, defaults to autosomes.
#' @param cutoff Set threshold for maximum CN state to be retrieved.
#' @param include_cn2 Optional boolean statement for including CN = 2 states in plot.
#'
#' @return Nothing.
#' @import tidyverse cowplot
#' @export
#'
#' @examples
#' chr1_cns = fancy_cnbar(this_sample = "HTMCP-01-06-00422-01A-01D", chr_select = c(1))
#' cns = fancy_cnbar(this_sample = "HTMCP-01-06-00422-01A-01D")
#'
fancy_cnbar = function(this_sample,
                       seq_data,
                       seq_path = NULL,
                       chrom_col = 2,
                       start_col = 3,
                       end_col = 4,
                       cn_col = 7,
                       plot_title = paste0(this_sample),
                       plot_subtitle = "n CNV Segments (barplots, left y-axis), n Affected bases for each CN state",
                       chr_select = paste0("chr", c(1:22)),
                       cutoff = 15,
                       include_cn2 = TRUE) {

  if(!missing(seq_data)){
    seq = seq_data
    seq = as.data.frame(seq)
    colnames(seq)[chrom_col] = "chrom"
    colnames(seq)[start_col] = "start"
    colnames(seq)[end_col] = "end"
    colnames(seq)[cn_col] = "CN"

  }else if (!is.null(seq_path)){
    seq = read.table(seq_path, sep = "\t", header = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
    seq = as.data.frame(seq)
    colnames(seq)[chrom_col] = "chrom"
    colnames(seq)[start_col] = "start"
    colnames(seq)[end_col] = "end"
    colnames(seq)[cn_col] = "CN"
  }

  #get maf data for a specific sample.
  if(missing(seq_data) && is.null(seq_path)){
    seq = get_sample_cn_segments(this_sample = this_sample, multiple_samples = FALSE, streamlined = FALSE, from_flatfile = TRUE)
  }

  #add chr prefix if missing
  if(!str_detect(seq$chrom, "chr")[2]){
    seq = mutate(seq, chrom = paste0("chr", chrom))
  }

  #read maf into R and select relevant variables and transformt to factor.
  seq_df = dplyr::select(seq, chrom, start, end, CN) %>%
    mutate_at(vars(chrom), list(factor))

  #subsetting maf based on user-defined parameters
  seq_df = seq_df[seq_df$chrom %in% chr_select, ]

  #transform data type
  seq_df$CN = as.factor(seq_df$CN)

  #count levels of factor
  if(include_cn2){
    cn_states = c(0:cutoff)
    cns_count = dplyr::filter(seq_df, CN %in% cn_states) %>%
      group_by(CN) %>%
      summarize(count = n())
  }else{
    cn_states = c(0:1, 3:cutoff)
    cns_count = dplyr::filter(seq_df, CN %in% cn_states) %>%
      group_by(CN) %>%
      summarize(count = n())
  }

  cns_count$Type = paste0(cns_count$CN)
  cns_count = dplyr::select(cns_count, count, Type)

  #compute lenght of cn segments and transform for plotting
  l_cn_seg = seq
  l_cn_seg$lenght = l_cn_seg$end - l_cn_seg$start
  l_cn_seg$CN = as.factor(l_cn_seg$CN)
  l_cn_seg = dplyr::filter(l_cn_seg, CN %in% cn_states)

  if(!include_cn2){
    l_cn_seg = dplyr::filter(l_cn_seg, CN != 2)
  }

  cn_seq_lenghts = aggregate(l_cn_seg$lenght, list(l_cn_seg$CN), sum)
  colnames(cn_seq_lenghts) = c("CN", "lenght")

  joined_cn = merge(cns_count, cn_seq_lenghts) %>%
    as.data.frame() %>%
    dplyr::select(CN, count, lenght)

  #get levels of cn states for plotting
  cn_levels = cns_count$Type

  #plot
  p = ggplot(joined_cn, aes(x = CN)) +
    geom_segment(aes(y = 1, yend = lenght/500000, x = CN, xend = CN)) +
    geom_point(aes(y = lenght/500000), colour = "#E6B315", size = 3, group = 2) +
    geom_bar(aes(y = count, fill = CN, label = count), position = "stack", stat = "identity") +
    scale_y_log10(limits = c(1, max(joined_cn$count) + 5000), sec.axis = sec_axis(~.*500000, name = "Nucleotides (n)")) +
    labs(title = plot_title, subtitle = plot_subtitle, x = "CN States", y = "CN Segments (n)", fill = "Legend") +
    scale_fill_manual(values = get_gambl_colours("copy_number")) +
    scale_x_discrete(limits = cn_levels) +
    geom_text(aes(x = CN, y = count, label = count), colour = "#000000", size = 5, position = position_stack(vjust = 0.5)) +
    theme_cowplot() +
    theme(legend.position = "none")

  return(p)
}


#' Generate a violine plot showing variant (SSM or SVs) size distributions for selected contigs.
#'
#' @param this_sample Sample to be plotted.
#' @param maf_data Optional parameter with maf like df already loaded into R.
#' @param maf_path Optional parameter with path to external maf like file.
#' @param ssm Set to FALSE to get plotting data from get_combined_sv (SVs). Default value is TRUE (plots SSM retrieved from annotate_cn_by_ssm$maf).
#' @param projection Genome build for returned variants (only applicable for ssm = FALSE).
#' @param min_vaf The minimum tumour VAF for a SV to be returned. Recommended: 0 (only applicable for ssm = FALSE).
#' @param variant_type_col Index of column holding Variant Type (to be used with either maf_data or maf_path).
#' @param chromosome_col Index of column holding Chromosome (to be used with either maf_data or maf_path).
#' @param start_col Index of column with variant start coordinates (to be used with either maf_data or maf_path).
#' @param end_col Index of column with variant end coordinates (to be used with either maf_data or maf_path).
#' @param plot_title Title of plot (default to sample ID).
#' @param plot_subtitle Subtitle for created plot.
#' @param scale_value Scale type for violin plot, accepted values are "area", "width", and "count", defualt is "count.
#' @param log_10 Boolean statement for yaxis, default is TRUE.
#' @param trim Boolean statment for trimming violin plot. Default is TRUE.
#' @param chr_select vector of chromosomes to be included in plot, defaults to autosomes.
#' @param coding_only Optional. Set to TRUE to restrict to plotting only coding mutations.
#' @param from_flatfile If set to true the function will use flat files instead of the database.
#' @param use_augmented Boolean statement if to use augmented maf, default is FALSE.
#'
#' @return plot as ggplot object.
#' @import tidyverse cowplot
#' @export
#'
#' @examples
#' plot = fancy_v_sizedis(this_sample = "HTMCP-01-06-00422-01A-01D")
#'
fancy_v_sizedis = function(this_sample,
                           maf_data,
                           maf_path = NULL,
                           ssm = TRUE,
                           projection = "grch37",
                           min_vaf = 0,
                           variant_type_col = 10,
                           chromosome_col = 5,
                           start_col = 6,
                           end_col = 7,
                           plot_title = paste0(this_sample),
                           plot_subtitle = "Variant Size Distribution",
                           scale_value = "width",
                           log_10 = TRUE,
                           plot_trim = FALSE,
                           chr_select = paste0("chr", c(1:22)),
                           coding_only = FALSE,
                           from_flatfile = TRUE,
                           use_augmented_maf = TRUE){

  if(!missing(maf_data)){
    maf = maf_data
    maf = as.data.frame(maf)
    colnames(maf)[variant_type_col] = "Variant_Type"
    colnames(maf)[chromosome_col] = "Chromosome"
    colnames(maf)[start_col] = "Start_Position"
    colnames(maf)[end_col] = "End_Position"

  }else if (!is.null(maf_path)){
    maf = fread_maf(maf_path)
    maf = as.data.frame(maf)
    colnames(maf)[variant_type_col] = "Variant_Type"
    colnames(maf)[chromosome_col] = "Chromosome"
    colnames(maf)[start_col] = "Start_Position"
    colnames(maf)[end_col] = "End_Position"
  }

  #get maf data for a specific sample.
  if(missing(maf_data) && is.null(maf_path)){
    if(ssm){
      maf = assign_cn_to_ssm(this_sample = this_sample, coding_only = coding_only, from_flatfile = from_flatfile, use_augmented_maf = use_augmented_maf)$maf
    }else{
      maf = get_combined_sv(sample_ids = this_sample, projection = projection, min_vaf = min_vaf) %>%
        dplyr::select(CHROM_A, START_A, END_A, manta_name)

      #get manta results in required format
      maf = data.frame(maf$CHROM_A, maf$START_A, maf$END_A, do.call(rbind, strsplit(maf$manta_name, split = ":", fixed = TRUE)))

      #rename variables
      names(maf)[1] = "Chromosome"
      names(maf)[2] = "Start_Position"
      names(maf)[3] = "End_Position"
      names(maf)[4] = "Variant_Type"

      #filter out translocations and set order of variables
      maf = dplyr::filter(maf, Variant_Type %in% c("MantaDEL", "MantaDUP")) %>%
        dplyr::select(Chromosome, Start_Position, End_Position, Variant_Type)

      #remove "Manta" from Variant_Type string
      maf$Variant_Type = gsub("^.{0,5}", "", maf$Variant_Type)
    }
  }

  #add chr prefix if missing
  if(!str_detect(maf$Chromosome, "chr")[1]){
    maf = mutate(maf, Chromosome = paste0("chr", Chromosome))
  }

  #read maf into R and select relevant variables and transform to factor.
  maf_df = dplyr::select(maf, Chromosome, Start_Position, End_Position, Variant_Type) %>%
    mutate_at(vars(Chromosome, Variant_Type), list(factor))

  #calculate variant size
  maf_df$Size = maf_df$End_Position - maf_df$Start_Position

  if(ssm){
    maf_df = maf_df[maf_df$Variant_Type %in% c("DEL", "INS"), ]
    levels(maf_df$Size)[levels(maf_df$Size) == "0"] = "1"
    maf_df[,5][maf_df[,5] == 0] <- 1
  }else{
    maf_df = maf_df[maf_df$Variant_Type %in% c("DEL", "DUP"), ]
  }

  maf_df$Size = as.integer(maf_df$Size)

  #sub-setting maf based on user-defined parameters
  maf_df = maf_df[maf_df$Chromosome %in% chr_select, ]

  p = ggplot(maf_df, aes(x = Variant_Type, y = Size, fill = Variant_Type)) +
        labs(title = plot_title, subtitle = plot_subtitle, x = "", y = "Variant Size (bp)") +
        geom_violin(trim = plot_trim, scale = scale_value, color = NA) +
        stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "black") +
        {if(ssm)scale_fill_manual(values = get_gambl_colours("indels"))} +
        {if(!ssm)scale_fill_manual(values = get_gambl_colours("svs"))} +
        {if(log_10)scale_y_log10()} +
        theme_cowplot() +
        theme(legend.position = "none")

  return(p)
}


#' Generate sample-level ideogram with copy number information, ssm and gene annotations, etc.
#'
#' @param this_sample Sample to be plotted (for multiple samples, see fancy_multisample_ideogram.
#' @param gene_annotation Annotate ideogram with a set of genes. These genes can either be specified as a vector of characters or a data frame.
#' @param seq_data Optional parameter with copy number df already loaded into R.
#' @param seq_path Optional parameter with path to external cn file.
#' @param maf_data Optional parameter with maf like df already loaded into R.
#' @param maf_path Optional parameter with path to external maf like file.
#' @param variant_type_col_maf Index of column holding Variant Type (to be used with either maf_data or maf_path).
#' @param chromosome_col_maf Index of column holding Chromosome (to be used with either maf_data or maf_path).
#' @param start_col_maf Index of column with variant start coordinates (to be used with either maf_data or maf_path).
#' @param end_col_maf Index of column with variant end coordinates (to be used with either maf_data or maf_path).
#' @param chrom_col_seq Index of column annotating Chromosome (to be used with either maf_data or maf_path).
#' @param start_col_seq Index of column with copy number start coordinates (to be used with either maf_data or maf_path).
#' @param end_col_seq Index of column with copy number end coordinates (to be used with either maf_data or maf_path).
#' @param cn_col Index of column holding copy number information (to be used with either maf_data or maf_path).
#' @param plot_title Title of plot (default to sample ID).
#' @param plot_subtitle Optional argument for plot subtitle.
#' @param intersect_regions Optional parameter for subset variant calls to specific regions. Should be either a vector of characters (chr:start-end) or data frame with regions.
#' @param include_ssm Set to TRUE to plot ssms (dels and ins).
#' @param ssm_count Optional parameter to summarize n variants per chromosome, inlcude_ssm must be set to TRUE.
#' @param coding_only Optional. Set to TRUE to restrict to plotting only coding mutations.
#' @param from_flatfile If set to true the function will use flat files instead of the database.
#' @param use_augmented Boolean statement if to use augmented maf, default is FALSE.
#'
#' @import data.table cowplot
#' @return Nothing.
#' @export
#'
#' @examples
#' #
#' fl_genes = dplyr::filter(lymphoma_genes, FL == TRUE) %>%
#'   dplyr::select(Gene) %>%
#'   pull(Gene)
#'
#' fl_genes_chr1 = gene_to_region(gene_symbol = fl_genes, return_as = "df") %>%
#'   dplyr::filter(chromosome == "1") %>%
#'   pull(hugo_symbol)
#'
#' ideogram_fl_chr1 = fancy_ideogram(this_sample = "HTMCP-01-06-00422-01A-01D",
#'                                   gene_annotation = fl_genes_chr1,
#'                                   intersect_regions = "chr1:10000-249250621",
#'                                   include_ssm = TRUE,
#'                                   ssm_count = TRUE,
#'                                   coding_only = FALSE,
#'                                   from_flatfile = FALSE,
#'                                   use_augmented = FALSE)
#'
#'  fl_regions = gene_to_region(gene_symbol = fl_genes, return_as = "df")
#'  ideogram_fl = fancy_ideogram(this_sample = "HTMCP-01-06-00422-01A-01D",
#'                               gene_annotation = fl_genes,
#'                               intersect_regions = fl_regions,
#'                               include_ssm = TRUE,
#'                               ssm_count = TRUE,
#'                               coding_only = FALSE,
#'                               from_flatfile = FALSE,
#'                               use_augmented = FALSE)
#'
fancy_ideogram = function(this_sample,
                          gene_annotation,
                          seq_data,
                          seq_path = NULL,
                          maf_data,
                          maf_path = NULL,
                          variant_type_col_maf = 10,
                          chromosome_col_maf = 5,
                          start_col_maf = 6,
                          end_col_maf = 7,
                          chrom_col_seq = 2,
                          start_col_seq = 3,
                          end_col_seq = 4,
                          cn_col_seq = 7,
                          plot_title = paste0(this_sample),
                          plot_subtitle = "Genome-wide Ideogram (grch37).",
                          intersect_regions,
                          include_ssm = TRUE,
                          ssm_count = TRUE,
                          coding_only = FALSE,
                          from_flatfile = TRUE,
                          use_augmented_maf = TRUE){

  #plot theme
  ideogram_theme = function(){
    theme(legend.position = "bottom", axis.title.x = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank())
  }

  #grch37 coordinates
  grch37_end = GAMBLR::chromosome_arms_grch37[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44),3]
  grch37_cent_start = GAMBLR::chromosome_arms_grch37[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43),3]
  grch37_cent_end = GAMBLR::chromosome_arms_grch37[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44),2]

  #additional regions to plot
  if(!missing(gene_annotation)){
    gene = gene_to_region(gene_symbol = gene_annotation, genome_build = "grch37", return_as = "df")
    gene.annotate = gene[gene$chr %in% paste0(c(1:22)), ]
    cols.int = c("chromosome", "start", "end")
    gene.annotate[cols.int] = sapply(gene.annotate[cols.int], as.integer)
  }

  #build chr table for segment plotting
  chr = paste0("chr", c(1:22))
  chr_start = c(0)
  chr_end = grch37_end
  cent_start = grch37_cent_start
  cent_end =  grch37_cent_end
  y = c(1:22)
  yend = c(1:22)

  #transform to data frame
  segment_data = data.frame(chr, chr_start, chr_end, cent_start, cent_end, y, yend)

  #load CN data
  if(!missing(seq_data)){
    cn_states = seq_data
    cn_states = as.data.frame(cn_states)
    colnames(cn_states)[chrom_col_seq] = "chrom"
    colnames(cn_states)[start_col_seq] = "start"
    colnames(cn_states)[end_col_seq] = "end"
    colnames(cn_states)[cn_col_seq] = "CN"

  }else if(!is.null(seq_path)){
    cn_states = read.table(seq_path, sep = "\t", header = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
    cn_states = as.data.frame(cn_states)
    colnames(cn_states)[chrom_col_seq] = "chrom"
    colnames(cn_states)[start_col_seq] = "start"
    colnames(cn_states)[end_col_seq] = "end"
    colnames(cn_states)[cn_col_seq] = "CN"
  }

  #get maf data for a specific sample.
  if(missing(seq_data) && is.null(seq_path)){
    cn_states = get_sample_cn_segments(this_sample_id = this_sample, multiple_samples = FALSE, with_chr_prefix = FALSE, streamlined = FALSE)
  }

  #convert chr into y coordinates
  cn_states$ycoord = cn_states$chrom

  #paste chr in chromosomecolumn, if not there
  if(!str_detect(cn_states$chrom, "chr")){
    cn_states = mutate(cn_states, chrom = paste0("chr", chrom))
  }

  if(!missing(intersect_regions)){
    #filter CN states on intersecting regions
    #transform regions to data tables
    #convenience function for converting intersect regions to a df (if it's a string) and renaming columns to match required format.
    if(is.list(intersect_regions)){
      colnames(intersect_regions)[1] = "chrom"
      colnames(intersect_regions)[2] = "start"
      colnames(intersect_regions)[3] = "end"
      if(!str_detect(intersect_regions$chrom, "chr")){
        intersect_regions = mutate(intersect_regions, chrom = paste0("chr", chrom))
      }
      intersect_regions = as.data.table(intersect_regions)
      intersect_regions$start = as.numeric(intersect_regions$start)
      intersect_regions$end = as.numeric(intersect_regions$end)
    }

    if(is.character(intersect_regions)){
      if(length(intersect_regions) > 1){
        message("Please only enter one region, only first region will be regarded. For mutiple regions, kindly provide a data frame with regions of interest")
      }

      split_chunks = unlist(strsplit(intersect_regions, ":"))
      split_chunks = unlist(strsplit(split_chunks, "-"))
      chrom = split_chunks[1]
      start = split_chunks[2]
      end = split_chunks[3]
      intersect_regions = cbind(chrom, start, end) %>%
        as.data.frame()

      intersect_regions$start = as.numeric(intersect_regions$start)
      intersect_regions$end = as.numeric(intersect_regions$end)

      if(!str_detect(intersect_regions$chrom, "chr")){
        intersect_regions = mutate(intersect_regions, chrom = paste0("chr", chrom))
      }
    }

    incoming_cn = as.data.table(cn_states)
    regions_sub = as.data.table(intersect_regions)

    #set keys
    data.table::setkey(incoming_cn, chrom, start, end)
    data.table::setkey(regions_sub, chrom, start, end)

    #intersect regions
    intersect = data.table::foverlaps(regions_sub, incoming_cn, nomatch = 0)

    #transform object to data frame
    inter_df = as.data.frame(intersect)

    #organize columns to match the expected format
    cn_states = select(inter_df, ID, chrom, start, end, LOH_flag, log.ratio, CN, ycoord)
  }

  #convert data types
  cols.int = c("start", "end", "ycoord")
  cn_states[cols.int] = sapply(cn_states[cols.int], as.integer)
  cn_states$chrom = as.factor(cn_states$chrom)
  cn_states$CN = as.factor(cn_states$CN)

  #subset on CN state
  cn_states$CN[cn_states$CN > 6] = 6
  cn_states = subset(cn_states, CN != 2)
  cn_states$CN = paste0("cn_", cn_states$CN)
  cn_states$CN = as.factor(cn_states$CN)
  cn_states = droplevels(cn_states)
  l = split(cn_states, cn_states$CN)
  list2env(l, envir = .GlobalEnv)

  #load maf data
  if(include_ssm){
    if(!missing(maf_data)){
      maf = maf_data
      maf = as.data.frame(maf)
      colnames(maf)[variant_type_col_maf] = "Variant_Type"
      colnames(maf)[chromosome_col_maf] = "Chromosome"
      colnames(maf)[start_col_maf] = "Start_Position"
      colnames(maf)[end_col_maf] = "End_Position"
    }else if (!is.null(maf_path)){
      maf = fread_maf(maf_path)
      maf = as.data.frame(maf)
      colnames(maf)[variant_type_col_maf] = "Variant_Type"
      colnames(maf)[chromosome_col_maf] = "Chromosome"
      colnames(maf)[start_col_maf] = "Start_Position"
      colnames(maf)[end_col_maf] = "End_Position"
    }

    if(missing(maf_data) && is.null(maf_path)){
      maf = assign_cn_to_ssm(this_sample = this_sample, coding_only = coding_only, from_flatfile = from_flatfile, use_augmented_maf = use_augmented_maf)$maf
    }

    #transform maf data
    maf_trans = dplyr::select(maf, Chromosome, Start_Position, End_Position, Variant_Type)

    #calculate mid points of variants
    maf_trans$mid = ((maf_trans$End_Position - maf_trans$Start_Position) / 2) + maf_trans$Start_Position

    #convert chr into y coordinates
    maf_trans$ystart = maf_trans$Chromosome
    maf_trans$yend = maf_trans$Chromosome

    #paste chr in maf, if not there
    if(!str_detect(maf_trans$Chromosome, "chr")){
      maf_trans = mutate(maf_trans, Chromosome = paste0("chr", Chromosome))
    }

    #convert data types
    maf_trans$Start_Position = as.double(maf_trans$Start_Position)
    maf_trans$End_Position = as.double(maf_trans$End_Position)
    maf_trans$ystart = as.integer(maf_trans$ystart)
    maf_trans$yend = as.integer(maf_trans$yend)

    if(!missing(intersect_regions)){
      #filter CN states on intersecting regions
      maf_tmp = maf_trans
      colnames(maf_tmp)[1] = "chrom"
      colnames(maf_tmp)[2] = "start"
      colnames(maf_tmp)[3] = "end"
      maf_tmp = dplyr::select(maf_tmp, chrom, start, end)
      maf.table = as.data.table(maf_tmp)
      data.table::setkey(maf.table, chrom, start, end)

      #intersect regions
      intersect_maf = data.table::foverlaps(regions_sub, maf.table, nomatch = 0)

      #transform object to data frame
      inter_maf_df = as.data.frame(intersect_maf)

      #rename columns
      colnames(inter_maf_df)[1] = "Chromosome"
      colnames(inter_maf_df)[2] = "Start_Position"
      colnames(inter_maf_df)[3] = "End_Position"

      #subset
      inter_maf_df = dplyr::select(inter_maf_df, Chromosome, Start_Position, End_Position)

      #perform a semi join with all cn states (to retain necessary columns)
      maf_trans = dplyr::semi_join(maf_trans, inter_maf_df)
    }

    #subset on variant type
    if(nrow(maf_trans > 0)){
      maf_del = dplyr::filter(maf_trans, Variant_Type == "DEL")
      maf_ins = dplyr::filter(maf_trans, Variant_Type == "INS")

      if(ssm_count){
        del_count = dplyr::filter(maf_trans, Variant_Type == "DEL") %>%
          add_count(Chromosome) %>%
          distinct(Chromosome, .keep_all = TRUE) %>%
          dplyr::select(Chromosome, Variant_Type, yend, n)

        ins_count = dplyr::filter(maf_trans, Variant_Type == "INS") %>%
          add_count(Chromosome) %>%
          distinct(Chromosome, .keep_all = TRUE) %>%
          dplyr::select(Chromosome, Variant_Type, yend, n)
      }
    }
  }

  #get colours and combine palette for indels and cn states
  ideogram_palette = c(get_gambl_colours("indels"), get_gambl_colours("copy_number"))
  selected_colours = ideogram_palette[c(1,2,19,18,16:13)]
  names(selected_colours)[c(3:8)] = c("CN0", "CN1", "CN3", "CN4", "CN5", "CN6+")

  #plot
  p = ggplot() +
    {if(include_ssm && nrow(maf_del > 0)) geom_segment(data = maf_del, aes(x = mid - 100000, xend = mid + 100000, y = ystart - 0.27, yend = yend - 0.27), color = "#53B1FC", size = 5, stat = "identity", position = position_dodge())} + #del
    {if(include_ssm && nrow(maf_del > 0)) geom_segment(data = maf_del, aes(x = mid, xend = mid, y = ystart - 0.35, yend = yend - 0.35), color = "black", lineend = "round", size = 3.5, stat = "identity", position = position_dodge())} + #del
    {if(include_ssm && nrow(maf_del > 0)) geom_segment(data = maf_del, aes(x = mid, xend = mid, y = ystart - 0.35, yend = yend - 0.35, color = "DEL"), lineend = "round", size = 3, stat = "identity", position = position_dodge())} + #del
    {if(include_ssm && nrow(maf_ins > 0)) geom_segment(data = maf_ins, aes(x = mid - 100000, xend = mid + 100000, y = ystart - 0.27, yend = yend - 0.27), color = "#FC9C6D", size = 5, stat = "identity", position = position_dodge())} + #ins
    {if(include_ssm && nrow(maf_ins > 0)) geom_segment(data = maf_ins, aes(x = mid, xend = mid, y = ystart - 0.35, yend = yend - 0.35), color = "black", lineend = "round", size = 3.5, stat = "identity", position = position_dodge())} + #ins
    {if(include_ssm && nrow(maf_ins > 0)) geom_segment(data = maf_ins, aes(x = mid, xend = mid, y = ystart - 0.35, yend = yend - 0.35, color = "INS"), lineend = "round", size = 3, stat = "identity", position = position_dodge())} + #ins
    {if(ssm_count && nrow(maf_del > 0)) annotate(geom = "text", x = -4000000, y = del_count$yend, label = del_count$n, color = "#3A8799", size = 3)} + #count del
    {if(ssm_count && nrow(maf_trans > 0)) annotate(geom = "text", x = -2300000, y = segment_data$y, label = " | ", color = "black", size = 3)} + #count sep
    {if(ssm_count && nrow(maf_ins > 0)) annotate(geom = "text", x = -1000000, y = ins_count$yend, label = ins_count$n, color = "#E6856F", size = 3)} + #count ins
    geom_segment(data = segment_data, aes(x = chr_start, xend = chr_end, y = y, yend = yend, label = chr), color = "#99A1A6", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) + #chr contigs
    geom_segment(data = segment_data, aes(x = cent_start, xend = cent_end, y = y, yend = yend), color = "white", size = 6, stat = "identity", position = position_dodge()) + #centromeres
    {if("cn_0" %in% levels(cn_states$CN)) geom_segment(data = cn_0, aes(x = start, xend = end, y = ycoord, yend = ycoord, color = "CN0"), size = 4.7, stat = "identity", position = position_dodge())} + #cn3
    {if("cn_1" %in% levels(cn_states$CN)) geom_segment(data = cn_1, aes(x = start, xend = end, y = ycoord, yend = ycoord, color = "CN1"), size = 4.7, stat = "identity", position = position_dodge())} + #cn3
    {if("cn_3" %in% levels(cn_states$CN)) geom_segment(data = cn_3, aes(x = start, xend = end, y = ycoord, yend = ycoord, color = "CN3"), size = 4.7, stat = "identity", position = position_dodge())} + #cn3
    {if("cn_4" %in% levels(cn_states$CN)) geom_segment(data = cn_4, aes(x = start, xend = end, y = ycoord, yend = ycoord, color = "CN4"), size = 4.7, stat = "identity", position = position_dodge())} + #cn4
    {if("cn_5" %in% levels(cn_states$CN)) geom_segment(data = cn_5, aes(x = start, xend = end, y = ycoord, yend = ycoord, color = "CN5"), size = 4.7, stat = "identity", position = position_dodge())} + #cn5
    {if("cn_6" %in% levels(cn_states$CN)) geom_segment(data = cn_6, aes(x = start, xend = end, y = ycoord, yend = ycoord, color = "CN6+"), size = 4.7, stat = "identity", position = position_dodge())} + #cn6 and more
    {if(!missing(gene_annotation)) geom_point(data = gene.annotate, aes(x = ((end - start) / 2) + start, y = chromosome - 0.28), shape = 25, color = "#A63932", fill = "#A63932", stat = "identity", position = position_dodge())} + #gene annotation
    {if(!missing(gene_annotation)) geom_label(data = gene.annotate, aes((x = end - start) / 2 + start, y = chromosome - 0.52, label = hugo_symbol), fontface = "bold", color = "white", fill = "#A63932", size = 3, check_overlap = TRUE)} + #gene annotation text
    geom_text(aes(x = -10000000 , y = yend, label = segment_data$chr), color = "black", size = 5) + #chr labels
    labs(title = plot_title, subtitle = plot_subtitle) + #plot titles
    scale_colour_manual(name = "", values = selected_colours) + #legend/colours
    scale_x_continuous(breaks = seq(0, max(segment_data$chr_end), by = 30000000)) + #set x-axis boundaries
    scale_y_reverse() + #reverse y axis
    theme_cowplot() + #themes
    ideogram_theme() #themes

  return(p)
}


#' Generate ideograms for selected sample, visualizing copy number variation segments. Also possible to only plot concordant (or discordant) cn segments between two samples. i.e how two samples differ, or are a like.
#'
#' @param these_sample_ids Sample to be plotted (accepts 2, 3 or 4 samples).
#' @param plot_title Main title of plot.
#' @param plot_sub Subtitle of plot.
#' @param chr_anno_dist Optional parameter to adjust chromosome annotations, default value is 3, increase to adjust annotations to left.
#' @param chr_select Optional parameter to subset plot to specific chromosomes. Default value is chr1-22.
#' @param include_cn2 Set to TRUE for plotting CN states == 2.
#' @param kompare Boolean statement, set to TRUE to call cnvKompare on the selected samples for plotting concordant (or discordant) cn segments across selected chromosomes.
#' @param concordance Boolean parameter to be used when kompare = TRUE. Default is TRUE, to plot discordant segments, set parameter to FALSE.
#' @param coding_only Optional. Set to TRUE to restrict to plotting only coding mutations.
#' @param from_flatfile If set to true the function will use flat files instead of the database.
#' @param use_augmented Boolean statement if to use augmented maf, default is FALSE.
#'
#' @return Nothing.
#' @export
#'
#' @examples
#'
#' two_samples = c("00-15201_tumorA", "00-15201_tumorB")
#' ideo_2_samp = fancy_multisamp_ideogram(these_sample_ids = two_samples,
#'                                        plot_title = "CN Segments Ideogram",
#'                                        plot_sub = "grch37",
#'                                        chr_anno_dist = 4,
#'                                        chr_select = paste0("chr", c(1:22)),
#'                                        coding_only = FALSE,
#'                                        from_flatfile = TRUE,
#'                                        use_augmented_maf = TRUE)
#'
fancy_multisamp_ideogram = function(these_sample_ids,
                                    plot_title = "CN Segments Ideogram",
                                    plot_sub = "grch37",
                                    chr_anno_dist = 3,
                                    chr_select = paste0("chr", c(1:22)),
                                    include_cn2 = FALSE,
                                    kompare = FALSE,
                                    concordance = TRUE,
                                    coding_only = FALSE,
                                    from_flatfile = TRUE,
                                    use_augmented_maf = TRUE){

  #plot theme
  ideogram_theme = function(){
    theme(legend.position = "bottom", axis.title.x = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank())}

  #transform chr_anno_dist and set some other variables
  anno_dist = chr_anno_dist * -10000000
  if(length(these_sample_ids) == 2){
    seg_dist = 0.18
    seg_size = 4
    seg_size_cent = 5
  }else if(length(these_sample_ids) == 3){
    seg_dist = 0.28
    seg_size = 3.5
    seg_size_cent = 4
  }else if(length(these_sample_ids) == 4){
    seg_dist = 0.12}

  #chr segment coordinates
  grch37_end = GAMBLR::chromosome_arms_grch37[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44),3]
  grch37_cent_start = GAMBLR::chromosome_arms_grch37[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43),3]
  grch37_cent_end = GAMBLR::chromosome_arms_grch37[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44),2]

  #build chr table for segment plotting
  chr = paste0("chr", c(1:22))
  chr_start = c(0)
  chr_end = grch37_end
  cent_start = grch37_cent_start
  cent_end =  grch37_cent_end
  y = c(1:22)
  yend = c(1:22)

  #transform to data frame
  segment_data = data.frame(chr, chr_start, chr_end, cent_start, cent_end, y, yend)
  segment_data$chr = as.factor(segment_data$chr)

  #sub-setting maf based on user-defined parameters
  segment_data = segment_data[segment_data$chr %in% chr_select, ]
  segment_data = droplevels(segment_data)

  if(kompare){
    #call cnvKompare to retreive CN segments shared (or not) shared between selected samples.
    cnv_komp = cnvKompare(sample_ids = these_sample_ids)

    #select concordant or discordant CN segments for plotting.
    if(concordance){
      cnv_cord = cnv_komp$concordant_cytobands
    }else{
      cnv_cord = cnv_komp$discordant_cytobands
    }

    #transform log.ratio to CN states
    cnv_cord$CN_tmp = 2*2^cnv_cord$log.ratio
    cnv_cord$CN = round(cnv_cord$CN_tmp) %>%
      as.factor()

    colnames(cnv_cord)[2] = "chrom"
    colnames(cnv_cord)[3] = "start"
    colnames(cnv_cord)[4] = "end"

    cn_states = dplyr::select(cnv_cord, ID, chrom, start, end, CN)
  }else{
    #load CN data
    cn_states = get_sample_cn_segments(multiple_samples = TRUE, sample_list = these_sample_ids, streamlined = FALSE)
  }

  #convert chr into y coordinates
  cn_states$ycoord = cn_states$chrom

  #paste chr in chromosomecolumn, if not there
  if(!str_detect(cn_states$chrom, "chr")){
    cn_states = mutate(cn_states, chrom = paste0("chr", chrom))}

  #transform data types
  cols.int = c("start", "end", "ycoord")
  cn_states[cols.int] = sapply(cn_states[cols.int], as.integer)
  cn_states$chrom = as.factor(cn_states$chrom)
  cn_states$ID = as.factor(cn_states$ID)

  #sub-setting maf based on user-defined parameters
  cn_states = cn_states[cn_states$chrom %in% chr_select, ]
  cn_states = droplevels(cn_states)

  #retrieve sample names
  samples = levels(cn_states$ID)

  #first sample
  sample1 = samples[1]
  sample1_cn = dplyr::filter(cn_states, ID == sample1)
  subset_cnstates(cn_segments = sample1_cn, samplen = 1, include_2 = include_cn2)
  sample1_cn$CN = as.factor(sample1_cn$CN)
  sample1_cn = droplevels(sample1_cn)

  #second sample
  sample2 = samples[2]
  sample2_cn = dplyr::filter(cn_states, ID == sample2)
  subset_cnstates(cn_segments = sample2_cn, samplen = 2, include_2 = include_cn2)
  sample2_cn$CN = as.factor(sample2_cn$CN)
  sample2_cn = droplevels(sample2_cn)

  #third sample (if provided...)
  if(length(these_sample_ids) > 2){
    sample3 = samples[3]
    sample3_cn = dplyr::filter(cn_states, ID == sample3)
    subset_cnstates(cn_segments = sample3_cn, samplen = 3, include_2 = include_cn2)
    sample3_cn$CN = as.factor(sample3_cn$CN)
    sample3_cn = droplevels(sample3_cn)}

  #fourth sample (if provided...)
  if(length(these_sample_ids) > 3){
    sample4 = samples[4]
    sample4_cn = dplyr::filter(cn_states, ID == sample4)
    subset_cnstates(cn_segments = sample4_cn, samplen = 4, include_2 = include_cn2)
    sample4_cn$CN = as.factor(sample4_cn$CN)
    sample4_cn = droplevels(sample4_cn)}

  #get colours and combine palette for indels and cn states
  ideogram_palette = get_gambl_colours("copy_number")
  if(include_cn2){
    selected_colours = ideogram_palette[c(17:11)]
    names(selected_colours)[c(1:6)] = c("CN0", "CN1", "CN2", "CN3", "CN4", "CN5", "CN6+")
  }else{
    selected_colours = ideogram_palette[c(17,16,14:11)]
    names(selected_colours)[c(1:6)] = c("CN0", "CN1", "CN3", "CN4", "CN5", "CN6+")
  }

  #plot
  if(length(these_sample_ids) >= 2 & length(these_sample_ids) < 4){
    p = ggplot() + geom_segment(data = segment_data, aes(x = chr_start, xend = chr_end, y = y + seg_dist, yend = yend + seg_dist, label = chr), color = "#99A1A6", lineend = "butt", size = seg_size, stat = "identity", position = position_dodge()) + #chr contigs
      geom_segment(data = segment_data, aes(x = cent_start, xend = cent_end, y = y + seg_dist, yend = yend + seg_dist), color = "white", size = seg_size_cent, stat = "identity", position = position_dodge()) + #centromeres
      annotate(geom = "text", x = -2000000, y = segment_data$yend + seg_dist, label = sample1, color = "black", size = 3, hjust = 1) + #sample name annotations
      {if("0" %in% levels(sample1_cn$CN)) geom_segment(data = cn_0_sample1, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN0"), size = seg_size, stat = "identity", position = position_dodge())} +  #cn0
      {if("1" %in% levels(sample1_cn$CN)) geom_segment(data = cn_1_sample1, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN1"), size = seg_size, stat = "identity", position = position_dodge())} +  #cn1
      {if("2" %in% levels(sample1_cn$CN)) geom_segment(data = cn_2_sample1, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN2"), size = seg_size, stat = "identity", position = position_dodge())} +  #cn2
      {if("3" %in% levels(sample1_cn$CN)) geom_segment(data = cn_3_sample1, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN3"), size = seg_size, stat = "identity", position = position_dodge())} + #cn3
      {if("4" %in% levels(sample1_cn$CN)) geom_segment(data = cn_4_sample1, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN4"), size = seg_size, stat = "identity", position = position_dodge())} + #cn4
      {if("5" %in% levels(sample1_cn$CN)) geom_segment(data = cn_5_sample1, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN5"), size = seg_size, stat = "identity", position = position_dodge())} + #cn5
      {if("6+" %in% levels(sample1_cn$CN)) geom_segment(data = cn_6_sample1, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN6+"), size = seg_size, stat = "identity", position = position_dodge())} + #second sample
      geom_segment(data = segment_data, aes(x = chr_start, xend = chr_end, y = y - seg_dist, yend = yend - seg_dist, label = chr), color = "#99A1A6", lineend = "butt", size = seg_size, stat = "identity", position = position_dodge()) + #chr contigs
      geom_segment(data = segment_data, aes(x = cent_start, xend = cent_end, y = y - seg_dist, yend = yend - seg_dist), color = "white", size = seg_size_cent, stat = "identity", position = position_dodge()) + #centromeres
      annotate(geom = "text", x = -2000000, y = segment_data$yend - seg_dist, label = sample2, color = "black", size = 3, hjust = 1) + #sample name annotations
      {if("0" %in% levels(sample2_cn$CN)) geom_segment(data = cn_0_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN0"), size = seg_size, stat = "identity", position = position_dodge())} +  #cn0
      {if("1" %in% levels(sample2_cn$CN)) geom_segment(data = cn_1_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN1"), size = seg_size, stat = "identity", position = position_dodge())} +  #cn1
      {if("2" %in% levels(sample2_cn$CN)) geom_segment(data = cn_2_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN2"), size = seg_size, stat = "identity", position = position_dodge())} +  #cn2
      {if("3" %in% levels(sample2_cn$CN)) geom_segment(data = cn_3_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN3"), size = seg_size, stat = "identity", position = position_dodge())} + #cn3
      {if("4" %in% levels(sample2_cn$CN)) geom_segment(data = cn_4_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN4"), size = seg_size, stat = "identity", position = position_dodge())} + #cn4
      {if("5" %in% levels(sample2_cn$CN)) geom_segment(data = cn_5_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN5"), size = seg_size, stat = "identity", position = position_dodge())} + #cn5
      {if("6+" %in% levels(sample2_cn$CN)) geom_segment(data = cn_6_sample2, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN6+"), size = seg_size, stat = "identity", position = position_dodge())} #third sample
    if(length(these_sample_ids) > 2){
      p = p + geom_segment(data = segment_data, aes(x = chr_start, xend = chr_end, y = y + seg_dist - 0.28, yend = yend + seg_dist - 0.28, label = chr), color = "#99A1A6", lineend = "butt", size = seg_size, stat = "identity", position = position_dodge()) + #chr contigs
        geom_segment(data = segment_data, aes(x = cent_start, xend = cent_end, y = y + seg_dist - 0.28, yend = yend + seg_dist - 0.28), color = "white", size = seg_size_cent, stat = "identity", position = position_dodge()) + #centromeres
        annotate(geom = "text", x = -2000000, y = segment_data$yend + seg_dist - 0.28, label = sample3, color = "black", size = 3, hjust = 1) + #sample name annotations
        {if("0" %in% levels(sample3_cn$CN)) geom_segment(data = cn_0_sample3, aes(x = start, xend = end, y = ycoord + seg_dist - 0.28, yend = ycoord + seg_dist - 0.28, color = "CN0"), size = seg_size, stat = "identity", position = position_dodge())} +  #cn0
        {if("1" %in% levels(sample3_cn$CN)) geom_segment(data = cn_1_sample3, aes(x = start, xend = end, y = ycoord + seg_dist - 0.28, yend = ycoord + seg_dist - 0.28, color = "CN1"), size = seg_size, stat = "identity", position = position_dodge())} +  #cn1
        {if("2" %in% levels(sample3_cn$CN)) geom_segment(data = cn_2_sample3, aes(x = start, xend = end, y = ycoord + seg_dist - 0.28, yend = ycoord + seg_dist - 0.28, color = "CN2"), size = seg_size, stat = "identity", position = position_dodge())} +  #cn2
        {if("3" %in% levels(sample3_cn$CN)) geom_segment(data = cn_3_sample3, aes(x = start, xend = end, y = ycoord + seg_dist - 0.28, yend = ycoord + seg_dist - 0.28, color = "CN3"), size = seg_size, stat = "identity", position = position_dodge())} + #cn3
        {if("4" %in% levels(sample3_cn$CN)) geom_segment(data = cn_4_sample3, aes(x = start, xend = end, y = ycoord + seg_dist - 0.28, yend = ycoord + seg_dist - 0.28, color = "CN4"), size = seg_size, stat = "identity", position = position_dodge())} + #cn4
        {if("5" %in% levels(sample3_cn$CN)) geom_segment(data = cn_5_sample3, aes(x = start, xend = end, y = ycoord + seg_dist - 0.28, yend = ycoord + seg_dist - 0.28, color = "CN5"), size = seg_size, stat = "identity", position = position_dodge())} + #cn5
        {if("6+" %in% levels(sample3_cn$CN)) geom_segment(data = cn_6_sample3, aes(x = start, xend = end, y = ycoord + seg_dist - 0.28, yend = ycoord + seg_dist - 0.28, color = "CN6+"), size = seg_size, stat = "identity", position = position_dodge())}} #cn6+
    p = p + geom_segment(data = segment_data, aes(x = chr_start, xend = chr_end, y = y - 0.5, yend = yend - 0.5, label = chr), color = "white", lineend = "butt", size = 2, stat = "identity", position = position_dodge()) + #white space between chromosome groups (upper)
      geom_segment(data = segment_data, aes(x = chr_start, xend = chr_end, y = y + 0.5, yend = yend + 0.5, label = chr), color = "white", lineend = "butt", size = 2, stat = "identity", position = position_dodge()) + #white space between chromosome groups (bottom)
      annotate(geom = "text", x = anno_dist, y = segment_data$yend, label = segment_data$chr, color = "black", size = 5, hjust = 1) + #chr labels
      labs(title = plot_title, subtitle = plot_sub) + #plot titles
      scale_colour_manual(name = "", values = selected_colours) + #colours and legend
      scale_x_continuous(breaks = seq(0, max(segment_data$chr_end), by = 30000000)) + #x-axis boundaries
      scale_y_reverse() + #flip ideogram
      theme_cowplot() +  #theme
      ideogram_theme() #more theme

    #plotting has its own plotting chunk, due to re-scaling of geom_segments and segment widths etc.
  }else if(length(these_sample_ids) == 4){
    p = ggplot() + geom_segment(data = segment_data, aes(x = chr_start, xend = chr_end, y = y - (seg_dist - 0.48), yend = yend - (seg_dist - 0.48), label = chr), color = "#99A1A6", lineend = "butt", size = 3, stat = "identity", position = position_dodge()) + #chr contigs
      geom_segment(data = segment_data, aes(x = cent_start, xend = cent_end, y = y - (seg_dist - 0.48), yend = yend - (seg_dist - 0.48)), color = "white", size = 3, stat = "identity", position = position_dodge()) + #centromeres
      annotate(geom = "text", x = -2000000, y = segment_data$yend - (seg_dist - 0.48), label = sample1, color = "black", size = 3, hjust = 1) + #sample name annotations
      {if("0" %in% levels(sample1_cn$CN)) geom_segment(data = cn_0_sample1, aes(x = start, xend = end, y = ycoord - (seg_dist - 0.48), yend = ycoord - (seg_dist - 0.48), color = "CN0"), size = 3, stat = "identity", position = position_dodge())} +  #cn0
      {if("1" %in% levels(sample1_cn$CN)) geom_segment(data = cn_1_sample1, aes(x = start, xend = end, y = ycoord - (seg_dist - 0.48), yend = ycoord - (seg_dist - 0.48), color = "CN1"), size = 3, stat = "identity", position = position_dodge())} +  #cn1
      {if("2" %in% levels(sample1_cn$CN)) geom_segment(data = cn_2_sample1, aes(x = start, xend = end, y = ycoord - (seg_dist - 0.48), yend = ycoord - (seg_dist - 0.48), color = "CN2"), size = 3, stat = "identity", position = position_dodge())} +  #cn2
      {if("3" %in% levels(sample1_cn$CN)) geom_segment(data = cn_3_sample1, aes(x = start, xend = end, y = ycoord - (seg_dist - 0.48), yend = ycoord - (seg_dist - 0.48), color = "CN3"), size = 3, stat = "identity", position = position_dodge())} + #cn3
      {if("4" %in% levels(sample1_cn$CN)) geom_segment(data = cn_4_sample1, aes(x = start, xend = end, y = ycoord - (seg_dist - 0.48), yend = ycoord - (seg_dist - 0.48), color = "CN4"), size = 3, stat = "identity", position = position_dodge())} + #cn4
      {if("5" %in% levels(sample1_cn$CN)) geom_segment(data = cn_5_sample1, aes(x = start, xend = end, y = ycoord - (seg_dist - 0.48), yend = ycoord - (seg_dist - 0.48), color = "CN5"), size = 3, stat = "identity", position = position_dodge())} + #cn5
      {if("6+" %in% levels(sample1_cn$CN)) geom_segment(data = cn_6_sample1, aes(x = start, xend = end, y = ycoord - (seg_dist - 0.48), yend = ycoord - (seg_dist - 0.48), color = "CN6+"), size = 3, stat = "identity", position = position_dodge())} + #cn6+
      geom_segment(data = segment_data, aes(x = chr_start, xend = chr_end, y = y - seg_dist, yend = yend - seg_dist, label = chr), color = "#99A1A6", lineend = "butt", size = 3, stat = "identity", position = position_dodge()) + #chr contigs
      geom_segment(data = segment_data, aes(x = cent_start, xend = cent_end, y = y - seg_dist, yend = yend - seg_dist), color = "white", size = 3, stat = "identity", position = position_dodge()) + #centromeres
      annotate(geom = "text", x = -2000000, y = segment_data$yend - seg_dist, label = sample2, color = "black", size = 3, hjust = 1) + #sample name annotations
      {if("0" %in% levels(sample2_cn$CN)) geom_segment(data = cn_0_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN0"), size = 3, stat = "identity", position = position_dodge())} +  #cn0
      {if("1" %in% levels(sample2_cn$CN)) geom_segment(data = cn_1_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN1"), size = 3, stat = "identity", position = position_dodge())} +  #cn1
      {if("2" %in% levels(sample2_cn$CN)) geom_segment(data = cn_2_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN2"), size = 3, stat = "identity", position = position_dodge())} +  #cn2
      {if("3" %in% levels(sample2_cn$CN)) geom_segment(data = cn_3_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN3"), size = 3, stat = "identity", position = position_dodge())} + #cn3
      {if("4" %in% levels(sample2_cn$CN)) geom_segment(data = cn_4_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN4"), size = 3, stat = "identity", position = position_dodge())} + #cn4
      {if("5" %in% levels(sample2_cn$CN)) geom_segment(data = cn_5_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN5"), size = 3, stat = "identity", position = position_dodge())} + #cn5
      {if("6+" %in% levels(sample2_cn$CN)) geom_segment(data = cn_6_sample2, aes(x = start, xend = end, y = ycoord - seg_dist, yend = ycoord - seg_dist, color = "CN6+"), size = 3, stat = "identity", position = position_dodge())} + #cn6+
      geom_segment(data = segment_data, aes(x = chr_start, xend = chr_end, y = y + seg_dist, yend = yend + seg_dist, label = chr), color = "#99A1A6", lineend = "butt", size = 3, stat = "identity", position = position_dodge()) + #chr contigs
      geom_segment(data = segment_data, aes(x = cent_start, xend = cent_end, y = y + seg_dist, yend = yend + seg_dist), color = "white", size = 3, stat = "identity", position = position_dodge()) + #centromeres
      annotate(geom = "text", x = -2000000, y = segment_data$yend + seg_dist, label = sample3, color = "black", size = 3, hjust = 1) + #sample name annotations
      {if("0" %in% levels(sample3_cn$CN)) geom_segment(data = cn_0_sample3, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN0"), size = 3, stat = "identity", position = position_dodge())} +  #cn0
      {if("1" %in% levels(sample3_cn$CN)) geom_segment(data = cn_1_sample3, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN1"), size = 3, stat = "identity", position = position_dodge())} +  #cn1
      {if("2" %in% levels(sample3_cn$CN)) geom_segment(data = cn_2_sample3, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN2"), size = 3, stat = "identity", position = position_dodge())} +  #cn2
      {if("3" %in% levels(sample3_cn$CN)) geom_segment(data = cn_3_sample3, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN3"), size = 3, stat = "identity", position = position_dodge())} + #cn3
      {if("4" %in% levels(sample3_cn$CN)) geom_segment(data = cn_4_sample3, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN4"), size = 3, stat = "identity", position = position_dodge())} + #cn4
      {if("5" %in% levels(sample3_cn$CN)) geom_segment(data = cn_5_sample3, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN5"), size = 3, stat = "identity", position = position_dodge())} + #cn5
      {if("6+" %in% levels(sample3_cn$CN)) geom_segment(data = cn_6_sample3, aes(x = start, xend = end, y = ycoord + seg_dist, yend = ycoord + seg_dist, color = "CN6+"), size = 3, stat = "identity", position = position_dodge())} + #cn6+
      geom_segment(data = segment_data, aes(x = chr_start, xend = chr_end, y = y + (seg_dist - 0.48), yend = yend + (seg_dist - 0.48), label = chr), color = "#99A1A6", lineend = "butt", size = 3, stat = "identity", position = position_dodge()) + #chr contigs
      geom_segment(data = segment_data, aes(x = cent_start, xend = cent_end, y = y + (seg_dist - 0.48), yend = yend  + (seg_dist - 0.48)), color = "white", size = 3, stat = "identity", position = position_dodge()) + #centromeres
      annotate(geom = "text", x = -2000000, y = segment_data$yend  + (seg_dist - 0.48), label = sample4, color = "black", size = 3, hjust = 1) + #sample name annotations
      {if("0" %in% levels(sample4_cn$CN)) geom_segment(data = cn_0_sample4, aes(x = start, xend = end, y = ycoord + (seg_dist - 0.48), yend = ycoord + (seg_dist - 0.48), color = "CN0"), size = 3, stat = "identity", position = position_dodge())} +  #cn0
      {if("1" %in% levels(sample4_cn$CN)) geom_segment(data = cn_1_sample4, aes(x = start, xend = end, y = ycoord + (seg_dist - 0.48), yend = ycoord + (seg_dist - 0.48), color = "CN1"), size = 3, stat = "identity", position = position_dodge())} +  #cn1
      {if("2" %in% levels(sample4_cn$CN)) geom_segment(data = cn_2_sample4, aes(x = start, xend = end, y = ycoord + (seg_dist - 0.48), yend = ycoord + (seg_dist - 0.48), color = "CN2"), size = 3, stat = "identity", position = position_dodge())} +  #cn2
      {if("3" %in% levels(sample4_cn$CN)) geom_segment(data = cn_3_sample4, aes(x = start, xend = end, y = ycoord + (seg_dist - 0.48), yend = ycoord + (seg_dist - 0.48), color = "CN3"), size = 3, stat = "identity", position = position_dodge())} + #cn3
      {if("4" %in% levels(sample4_cn$CN)) geom_segment(data = cn_4_sample4, aes(x = start, xend = end, y = ycoord + (seg_dist - 0.48), yend = ycoord + (seg_dist - 0.48), color = "CN4"), size = 3, stat = "identity", position = position_dodge())} + #cn4
      {if("5" %in% levels(sample4_cn$CN)) geom_segment(data = cn_5_sample4, aes(x = start, xend = end, y = ycoord + (seg_dist - 0.48), yend = ycoord + (seg_dist - 0.48), color = "CN5"), size = 3, stat = "identity", position = position_dodge())} + #cn5
      {if("6+" %in% levels(sample4_cn$CN)) geom_segment(data = cn_6_sample4, aes(x = start, xend = end, y = ycoord + (seg_dist - 0.48), yend = ycoord + (seg_dist - 0.48), color = "CN6+"), size = 3, stat = "identity", position = position_dodge())} + #cn6+
      geom_segment(data = segment_data, aes(x = chr_start, xend = chr_end, y = y - 0.5, yend = yend - 0.5, label = chr), color = "white", lineend = "butt", size = 2, stat = "identity", position = position_dodge()) + #white space between chromosome groups (upper)
      geom_segment(data = segment_data, aes(x = chr_start, xend = chr_end, y = y + 0.5, yend = yend + 0.5, label = chr), color = "white", lineend = "butt", size = 2, stat = "identity", position = position_dodge()) + #white space between chromosome groups (bottom)
      annotate(geom = "text", x = anno_dist, y = segment_data$yend, label = segment_data$chr, color = "black", size = 5, hjust = 1) + #chr labels
      labs(title = plot_title, subtitle = plot_sub) + #plot titles
      scale_colour_manual(name = "", values = selected_colours) + #colours and legend
      scale_x_continuous(breaks = seq(0, max(segment_data$chr_end), by = 30000000)) + #x-axis boundaries
      scale_y_reverse() + #flip ideogram
      theme_cowplot() +  #theme
      ideogram_theme() #more theme
  }
  return(p)
}


#' Construct pdf with sample-level plots, using minimum of arguments
#'
#' @param this_sample Sample ID to be plotted in report.
#' @param export_individual_plots Boolean parameter, set to TRUE to export individual plots.
#' @param out Path to output folder.
#' @param seq_data Optional parameter with copy number df already loaded into R.
#' @param seq_path Optional parameter with path to external cn file.
#' @param maf_data Optional parameter with maf like df already loaded into R.
#' @param maf_path Optional parameter with path to external maf like file.
#'
#' @return Nothing.
#' @export
#'
#' @examples
#' comp_report(this_sample = "HTMCP-01-06-00422-01A-01D", out = "reports/", export_individual_plots = TRUE)
#'
comp_report = function(this_sample,
                       export_individual_plots = FALSE,
                       out,
                       seq_data,
                       seq_path = NULL,
                       maf_data,
                       maf_path = NULL){

  if(!missing(maf_data)){
    maf = maf_data
    maf = as.data.frame(maf)
    colnames(maf)[variant_type_col_maf] = "Variant_Type"
    colnames(maf)[chromosome_col_maf] = "Chromosome"
    colnames(maf)[start_col_maf] = "Start_Position"
    colnames(maf)[end_col_maf] = "End_Position"

  }else if (!is.null(maf_path)){
    maf = fread_maf(maf_path)
    maf = as.data.frame(maf)
    colnames(maf)[variant_type_col_maf] = "Variant_Type"
    colnames(maf)[chromosome_col_maf] = "Chromosome"
    colnames(maf)[start_col_maf] = "Start_Position"
    colnames(maf)[end_col_maf] = "End_Position"
  }

  if(!missing(seq_data)){
    seq = seq_data
    seq = as.data.frame(seq)
    colnames(seq)[chrom_col_seq] = "chrom"
    colnames(seq)[start_col_seq] = "start"
    colnames(seq)[end_col_seq] = "end"
    colnames(seq)[cn_col_seq] = "CN"

  }else if (!is.null(seq_path)){
    seq = read.table(seq_path, sep = "\t", header = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
    seq = as.data.frame(seq)
    colnames(seq)[chrom_col_seq] = "chrom"
    colnames(seq)[start_col_seq] = "start"
    colnames(seq)[end_col_seq] = "end"
    colnames(seq)[cn_col_seq] = "CN"
  }

  #read maf and seq data into r (avoid calling assign_cn_to_ssm and get_cn_segments for every plotting function)
  if(missing(maf_data) && is.null(maf_path)){
    maf = assign_cn_to_ssm(this_sample = this_sample, coding_only = FALSE, from_flatfile = TRUE, use_augmented_maf = TRUE)$maf
  }

  if(missing(seq_data) && is.null(seq_path)){
    seq = get_sample_cn_segments(this_sample = this_sample, multiple_samples = FALSE, streamlined = FALSE, from_flatfile = TRUE)
  }

  #execute a collection of sample-level plots with default parameters
  #page 1
  ssm_chr = fancy_v_chrcount(this_sample = this_sample, maf_data = maf, plot_title = "", plot_subtitle = "A. SSM Distribution Per Chromosome.", hide_legend = TRUE)
  sv_chr = fancy_v_chrcount(this_sample = this_sample, plot_title = "", plot_subtitle = "B. SV Distribution Per Chromosome.", ssm = FALSE, hide_legend = TRUE)
  ssm_count = fancy_v_count(this_sample = this_sample,  maf_data = maf, plot_title = "", plot_subtitle = "C. SSM Counts.", hide_legend = TRUE)
  violine_plot = fancy_v_sizedis(this_sample = this_sample,  maf_data = maf, plot_title = "", plot_subtitle = "D. SSM Size Distributions.")
  sv_count = fancy_v_count(this_sample = this_sample, plot_title = "", plot_subtitle = "E. SV Counts.", ssm = FALSE, variant_select = c("DEL", "DUP"), hide_legend = TRUE)
  sv_size = fancy_sv_sizedens(this_sample = this_sample, plot_title = "", plot_subtitle = "F. SV Size Density.", hide_legend = TRUE)
  snv_plot = fancy_snv_chrdistplot(this_sample = this_sample,  maf_data = maf, plot_title = "", plot_subtitle = "G. SNV Distribution Per Chromosome.")
  cns = fancy_cnbar(this_sample = this_sample, seq_data = seq, plot_title = "", plot_subtitle = "H. CN states.")

  #page 2 ideogram
  cnv_ideogram = fancy_ideogram(this_sample = this_sample, seq_data = seq, maf_data = maf, plot_title = "", plot_subtitle = "F. Ideogram.")

  #build pdf report
  pdf(paste0(out, this_sample, "_report.pdf"), width = 17, height = 12)
  page1 = grid.arrange(ssm_chr, sv_chr, ssm_count, violine_plot, sv_count, sv_size, snv_plot, cns, nrow = 3, ncol = 6, name = "Report", top = textGrob(paste0(this_sample, " - Report"), gp = gpar(fontsize = 15, fontface = "bold")), bottom = "Page 1", layout_matrix = rbind(c(1,1,1,2,2,2), c(3,3,4,4,5,5), c(6,6,7,7,8,8)))
  page2 = grid.arrange(cnv_ideogram,  nrow = 4, ncol = 4, name = "Report", bottom = "Page 2", layout_matrix = rbind(c(1,1,1,1), c(1,1,1,1), c(1,1,1,1), c(1,1,1,1)))
  dev.off()

  #export individual plots
  if(export_individual_plots){
    ggsave(ssm_chr, file = paste0(out, this_sample, "_ssm_dist_chr.pdf"), limitsize = FALSE, width = 17, height = 12, units = c("in"), dpi = 300)
    ggsave(sv_chr, file = paste0(out, this_sample, "_sv_dist_chr.pdf"), limitsize = FALSE, width = 17, height = 12, units = c("in"), dpi = 300)
    ggsave(snv_plot, file = paste0(out, this_sample, "_snv_dist_chr.pdf"), limitsize = FALSE, width = 17, height = 12, units = c("in"), dpi = 300)
    ggsave(ssm_count, file = paste0(out, this_sample, "_ssm_counts.pdf"), limitsize = FALSE, width = 12, height = 12, units = c("in"), dpi = 300)
    ggsave(sv_count, file = paste0(out, this_sample, "_sv_counts.pdf"), limitsize = FALSE, width = 12, height = 12, units = c("in"), dpi = 300)
    ggsave(sv_size, file = paste0(out, this_sample, "_sv_size_dens.pdf"), limitsize = FALSE, width = 12, height = 12, units = c("in"), dpi = 300)
    ggsave(cns, file = paste0(out, this_sample, "_cn_states.pdf"), limitsize = FALSE, width = 12, height = 12, units = c("in"), dpi = 300)
    ggsave(violine_plot, file = paste0(out, this_sample, "_sv_size_dist.pdf"), limitsize = FALSE, width = 12, height = 12, units = c("in"), dpi = 300)
    ggsave(cnv_ideogram, file = paste0(out, this_sample, "_cnv_ideo.pdf"), limitsize = FALSE, width = 17, height = 12, units = c("in"), dpi = 300)
  }
  return()
}


#' Create a circos plots visualizing SVS and SSM with optional gene annotations.
#'
#' @param this_sample Sample to be plotted.
#' @param gene_list Optional parameter to annotate genes on the circos plot from a list of genes (df). Is compatible with gene_to_region (return_as = "bed") output format. See examples.
#' @param ssm_calls Boolean parameter for plotting ssm. Default is TRUE.
#' @param sv_calls Boolean parameter for plotting SVs, default is TRUE.
#' @param chr_select Optional argument for subset on selected chromosomes, default is all autosomes.
#' @param vaf_cutoff Threshold for filtering variants on VAF (events with a VAF > cutoff will be retained).
#' @param coding_only Optional. Set to TRUE to restrict to plotting only coding mutations.
#' @param from_flatfile If set to TRUE the function will use flat files instead of the database.
#' @param use_augmented_maf Boolean statement if to use augmented maf, default is TRUE.
#' @param projection Genomic projection for variants and circos plot. Accepted values are grch37 and hg38, default is grch37.
#' @param out Path to output folder, to where the plot will be exported.
#' @param plot_title Optional parameter for naming your plot, default is this_sample.
#' @param pdf Set to FALSE for png, default is TRUE (pdf).
#' @param file_name Optional parameter for specifying the file name of generated circos plot, default is "{this_sample}_circos.pdf". If pdf is set to FALSE, a png will be generated, thus the .png extension needs to be attached to the file_name.
#'
#' @return Nothing.
#' @import tidyverse RCircos
#' @export
#'
#' @examples
#' #retrieve gene names for FL genes
#' fl_genes = dplyr::filter(lymphoma_genes, FL == TRUE) %>%
#'  dplyr::select(Gene) %>%
#'  pull(Gene)
#'
#' # get regions for selected genes
#' fl_genes_list = gene_to_region(gene_symbol = fl_genes, return_as = "bed")
#'
#' fancy_circos_plot_new(this_sample = "DOHH-2",
#'                       ssm_calls = FALSE,
#'                       sv_calls = TRUE,
#'                       gene_list = fl_genes_list,
#'                       chr_select = c("chr8", "chr14", "chr18"),
#'                       coding_only = FALSE,
#'                       projection = "grch37",
#'                       out = "../../plots/",
#'                       plot_title = "DOHH-2 (SVs) Example Plot",
#'                       pdf = FALSE,
#'                       file_name = "dohh2_example.png")
#'
fancy_circos_plot = function(this_sample,
                             gene_list,
                             ssm_calls = TRUE,
                             sv_calls = TRUE,
                             chr_select = paste0("chr", c(1:22)),
                             vaf_cutoff = 0,
                             coding_only = FALSE,
                             from_flatfile = TRUE,
                             use_augmented_maf = TRUE,
                             projection = "grch37",
                             plot_title = paste0(this_sample),
                             out,
                             pdf = TRUE,
                             file_name = paste0(this_sample, "_circos.pdf")){

  #set track properties based on selected plotting data
  if(ssm_calls && sv_calls && !missing(gene_list)){
    gene_con_track = 1
    gene_name_track = 2
    sv_del_track = 4
    sv_dup_track = 5
    ssm_snp_track = 6
    ssm_dnp_track = 7
    ssm_del_track = 8
    ssm_ins_track = 9
    trans_track = 10
  }

  if(ssm_calls && sv_calls && missing(gene_list)){
    sv_del_track = 1
    sv_dup_track = 2
    ssm_snp_track = 3
    ssm_dnp_track = 4
    ssm_del_track = 5
    ssm_ins_track = 6
    trans_track = 7
  }

  if(ssm_calls && !sv_calls && !missing(gene_list)){
    gene_con_track = 1
    gene_name_track = 2
    ssm_snp_track = 4
    ssm_dnp_track = 5
    ssm_del_track = 6
    ssm_ins_track = 7
  }

  if(ssm_calls && !sv_calls && missing(gene_list)){
    ssm_snp_track = 1
    ssm_dnp_track = 2
    ssm_del_track = 3
    ssm_ins_track = 4
  }

  if(!ssm_calls && sv_calls && !missing(gene_list)){
    gene_con_track = 1
    gene_name_track = 2
    sv_del_track = 4
    sv_dup_track = 5
    trans_track = 6
  }

  if(!ssm_calls && sv_calls && missing(gene_list)){
    sv_del_track = 1
    sv_dup_track = 2
    trans_track = 3
  }

  #sanity checking incoming gene list and renaming columns if needed
  if(!missing(gene_list)){
    #check type of incoming gene_list
    if(!is.list(gene_list)){
      message("Please ensure that incoming gene list is in fact a data frame with the following columns: chr:start:end:gene")
    }

    #rename columns, if needed (Rcircos is expecting the column names to always be the same...)
    if(!"Chromosome" %in% colnames(gene_list)[1]){
      colnames(gene_list)[1] = "Chromosome"
    }
    if(!"chromStart" %in% colnames(gene_list)[2]){
      colnames(gene_list)[2] = "chromStart"
    }
    if(!"chromEnd" %in% colnames(gene_list)[3]){
      colnames(gene_list)[3] = "chromEnd"
    }
    if(!"Gene" %in% colnames(gene_list)[4]){
      colnames(gene_list)[4] = "Gene"
    }

    #add "chr" prefix, if needed
    if(!str_detect(gene_list$Chromosome, "chr")){
      gene_list = mutate(gene_list, Chromosome = paste0("chr", Chromosome))
    }

    #filter gene list on selected chromosomes
    gene_list = gene_list[gene_list$Chromosome %in% chr_select, ]
  }

  #get SSM data
  if(ssm_calls){
    maf = assign_cn_to_ssm(this_sample = this_sample, coding_only = coding_only, from_flatfile = from_flatfile, use_augmented_maf = use_augmented_maf)$maf #get maf data
    maf_tmp = dplyr::select(maf, Chromosome, Start_Position, End_Position, Variant_Type) #select appropriate columns
    maf_tmp$Variant_Size = maf_tmp$End_Position - maf_tmp$Start_Position # calcualte variant size
    maf_tmp$Variant_Type = as.factor(maf_tmp$Variant_Type) #transform Variant_Type to factor
    maf_tmp[maf_tmp==0] <- 1 #transform all lenght coordinates == 0 to 1

    if(!str_detect(maf_tmp$Chromosome, "chr")){ #add chr prefic, if missing...
      maf_tmp = mutate(maf_tmp, Chromosome = paste0("chr", Chromosome))
    }

    maf_tmp = maf_tmp[maf_tmp$Chromosome %in% chr_select, ] #filter incoming maf on selected chromosomes
    ssm_del = dplyr::filter(maf_tmp, Variant_Type == "DEL") #subset on deletions
    ssm_ins = dplyr::filter(maf_tmp, Variant_Type == "INS") #subset on insertions
    ssm_snp = dplyr::filter(maf_tmp, Variant_Type == "SNP") #subset on single nucleotide polymorphism
    ssm_dnp = dplyr::filter(maf_tmp, Variant_Type == "DNP") #subset on dinucleotide polymorphism
    message(paste0(nrow(ssm_del) + nrow(ssm_dnp) + nrow(ssm_ins) + nrow(ssm_snp)), " SSMs found for ", this_sample)
  }

  #get SVs
  if(sv_calls){
    svs = get_combined_sv(sample_ids = this_sample, projection = projection)

    #filter on vaf
    svs = dplyr::filter(svs, VAF_tumour > vaf_cutoff)

    #subset on relevant variables
    svs_df = dplyr::select(svs, CHROM_A, START_A, END_A, CHROM_B, START_B, END_B, manta_name)

    #split manta_name variable
    svs_df = data.frame(svs_df$CHROM_A, svs_df$START_A, svs_df$END_A, svs_df$CHROM_B, svs_df$START_B, svs_df$END_B, do.call(rbind, strsplit(svs_df$manta_name, split = ":", fixed = TRUE)))

    #rename variables
    colnames(svs_df)[1:7] = c("CHROM_A", "START_A", "END_A", "CHROM_B", "START_B", "END_B", "TYPE")

      #add chr prefix, if missing
    if(!str_detect(svs_df$CHROM_A, "chr")[1]){
      svs_df = mutate(svs_df, CHROM_A = paste0("chr", CHROM_A))
      }

    if(!str_detect(svs_df$CHROM_B, "chr")[4]){
      svs_df = mutate(svs_df, CHROM_B = paste0("chr", CHROM_B))
      }

    #subset on selected chromosomes
    svs_df = svs_df[svs_df$CHROM_A %in% chr_select, ]
    svs_df = svs_df[svs_df$CHROM_B %in% chr_select, ]

    #subset df on SV type
    sv_trans = dplyr::filter(svs_df, TYPE == "MantaBND") %>%
      dplyr::select(CHROM_A, START_A, END_A, CHROM_B, START_B, END_B)

    sv_del = dplyr::filter(svs_df, TYPE == "MantaDEL") %>%
      dplyr::select(CHROM_A, START_A, END_A)

    sv_dup = dplyr::filter(svs_df, TYPE == "MantaDUP") %>%
      dplyr::select(CHROM_A, START_A, END_A)

    #calculate sizes
    sv_del$SIZE = sv_del$END_A - sv_del$START_A
    sv_dup$SIZE = sv_dup$END_A - sv_dup$START_A

    message(paste0(nrow(sv_trans) + nrow(sv_del) + nrow(sv_dup)), " SVs found for ", this_sample)
  }

  #plotting
  #define reference build
  if(projection == "grch37"){
    data(UCSC.HG19.Human.CytoBandIdeogram)
    cytobands = UCSC.HG19.Human.CytoBandIdeogram
  }else if(projection == "hg38"){
    data(UCSC.HG38.Human.CytoBandIdeogram)
    cytobands = UCSC.HG38.Human.CytoBandIdeogram
  }

  #get chr excluded (reversed of chr included, since RCircos only accept this)
  chr_all = paste0("chr", c(1:22, "X", "Y"))
  chr_exclude = setdiff(chr_all, chr_select)

  #set core components
  suppressMessages(RCircos::RCircos.Set.Core.Components(cyto.info = cytobands, chr.exclude = chr_exclude, tracks.inside = 10))

  #set plot parameters
  RCircos.params = RCircos.Get.Plot.Parameters()

  #define plotting parameters
  out.file = paste0(out, file_name)

  if(pdf){
    pdf(out.file, height = 7, width = 7)
  }else{
    png(out.file, height = 7, width = 7, units = "in", res = 300)
  }

  RCircos.Set.Plot.Area(margins = 0);

  #create empty plot
  RCircos.Chromosome.Ideogram.Plot()

  #add tracks to plot
  #add gene names
  if(!missing(gene_list)){
    RCircos.Gene.Connector.Plot(gene_list, gene_con_track, "in")
    RCircos.Gene.Name.Plot(gene_list, 4, gene_name_track, "in")
  }

  #aadd tracks
  if(ssm_calls){
    #ssm deletions
    RCircos.params$track.background = "steelblue2"
    RCircos.params$max.layers = 1
    RCircos.Reset.Plot.Parameters(RCircos.params)
    RCircos.Tile.Plot(ssm_del, track.num = ssm_del_track, side = "in")

    #ssm insertions
    RCircos.params$track.background = "sienna2"
    RCircos.params$max.layers = 1
    RCircos.Reset.Plot.Parameters(RCircos.params)
    RCircos.Tile.Plot(ssm_ins, track.num = ssm_ins_track, side = "in")

    #ssm snp
    RCircos.params$track.background = "seagreen"
    RCircos.params$max.layers = 1
    RCircos.Reset.Plot.Parameters(RCircos.params)
    RCircos.Tile.Plot(ssm_snp, track.num = ssm_snp_track, side = "in")

    #ssm dnp
    RCircos.params$track.background = "tomato4"
    RCircos.params$max.layers = 1
    RCircos.Reset.Plot.Parameters(RCircos.params)
    RCircos.Tile.Plot(ssm_dnp, track.num = ssm_dnp_track, side = "in")
  }

  if(sv_calls){
    #translocations
    RCircos.Link.Plot(sv_trans, track.num = trans_track, by.chromosome = FALSE)

    #duplications
    RCircos.params$track.background = "sienna2"
    RCircos.params$max.layers = 1
    RCircos.Reset.Plot.Parameters(RCircos.params)
    RCircos.Tile.Plot(sv_dup, track.num = sv_dup_track, side = "in")

    #deletions
    RCircos.params$track.background = "steelblue2"
    RCircos.params$max.layers = 1
    RCircos.Reset.Plot.Parameters(RCircos.params)
    RCircos.Tile.Plot(sv_del, track.num = sv_del_track, side = "in")
  }

  #add plot title and legends
  if(sv_calls && !ssm_calls){
    text(x = 0, y = 2.5, plot_title, font = 2, cex = 1.2)
    legend(x = 2, y = -1.2, legend = c("Del", "Dup"), bty = "n", col = c("steelblue2", "sienna2"), pch = 19, text.col = "black", horiz = FALSE, pt.cex = 2, title = "SVs", inset = c(0.1, 0.1))
  }

  if(!sv_calls && ssm_calls){
    text(x = 0, y = 2.5, plot_title, font = 2, cex = 1.2)
    legend(x = 2, y = -1.8, legend = c("SNP", "DNP", "Del", "Ins"), bty = "n", col = c("seagreen", "tomato4", "steelblue2", "sienna2"), pch = 19, text.col = "black", horiz = FALSE, pt.cex = 2, title = "SSM", inset = c(0.1, 0.1))
  }

  if(sv_calls && ssm_calls){
    text(x = 0, y = 2.5, plot_title, font = 2, cex = 1.2)
    legend(x = 2, y = -1.2, legend = c("Del", "Dup"), bty = "n", col = c("steelblue2", "sienna2"), pch = 19, text.col = "black", horiz = FALSE, pt.cex = 2, title = "SVs", inset = c(0.1, 0.1))
    legend(x = 2, y = -1.8, legend = c("SNP", "DNP", "Del", "Ins"), bty = "n", col = c("seagreen", "tomato4", "steelblue2", "sienna2"), pch = 19, text.col = "black", horiz = FALSE, pt.cex = 2, title = "SSM", inset = c(0.1, 0.1))
  }

  invisible(dev.off())
}


#' Generate plot visualizing SV sizes. Subset on variant type, filter on VAF, size etc.
#'
#' @param this_sample Sample to be plotted.
#' @param maf_data Optional parameter with copy number df already loaded into R.
#' @param maf_path Optional parameter with path to external cn file.
#' @param chrom_a_col Index of column holding chromosome (to be used with either maf_data or maf_path).
#' @param start_a_col Index of column holding start coordiantes (to be used with either maf_data or maf_path).
#' @param end_a_col Index of column holding end coordinates (to be used with either maf_data or maf_path).
#' @param variant_type_col Index of column holding variant type information (to be used with either maf_data or maf_path).
#' @param vaf_cutoff Threshold for filtering variants on VAF (events with a VAF > cutoff will be retained).
#' @param size_cutoff Threshold for filtering variants on size, default is 50bp.
#' @param adjust_value A multiplicate bandwidth adjustment. This makes it possible to adjust the bandwidth while still using the a bandwidth estimator. For example, adjust = 1/2 means use half of the default bandwidth.
#' @param trim If FALSE, the default, each density is computed on the full range of the data.
#' @param chr_select Optional argument for subsetting on selected chromosomes, default is all autosomes.
#' @param hide_legend Set to True to remove legend from plot, default is FALSE.
#' @param plot_title Title of plot (default to sample ID).
#' @param plot_subtitle Subtitle for created plot.
#' @param projection Genomic projection for SVs and circos plot. Accepted values are grch37 and hg38.
#'
#' @return Nothing.
#' @import tidyverse
#' @export
#'
#' @examples
#' myplot = fancy_sv_sizedens(this_sample = "HTMCP-01-06-00422-01A-01D")
#' myplot2 = fancy_sv_sizedens(this_sample = "HTMCP-01-06-00422-01A-01D", size_cutoff = 0, chr_select = c("chr1", "chr2"))
#'
fancy_sv_sizedens = function(this_sample,
                            maf_data,
                            maf_path = NULL,
                            chrom_a_col = 3,
                            start_a_col = 4,
                            end_a_col = 5,
                            variant_type_col = 9,
                            vaf_cutoff = 0,
                            size_cutoff = 50,
                            adjust_value = 1,
                            trim = FALSE,
                            hide_legend = FALSE,
                            chr_select = paste0("chr", c(1:22)),
                            plot_title = paste0(this_sample),
                            plot_subtitle = paste0("SV sizes for Manta calls. Dashed line annotates mean variant size.\nVAF cut off: ", vaf_cutoff,", SV size cut off: ", size_cutoff),
                            projection = "grch37"){
  if(!missing(maf_data)){
    svs = maf_data
    svs = as.data.frame(svs)
    colnames(svs)[chrom_a_col] = "CHROM_A"
    colnames(svs)[start_a_col] = "START_A"
    colnames(svs)[end_a_col] = "END_A"
    colnames(svs)[variant_type_col] = "manta_name"

  }else if (!is.null(maf_path)){
    svs = maf_data
    svs = as.data.frame(svs)
    colnames(svs)[chrom_a_col] = "CHROM_A"
    colnames(svs)[start_a_col] = "START_A"
    colnames(svs)[end_a_col] = "END_A"
    colnames(svs)[variant_type_col] = "manta_name"
  }

  #get variants, filter and subset
  if(missing(maf_data) && is.null(maf_path)){
    svs = get_combined_sv(sample_ids = this_sample, projection = projection) %>%
      dplyr::filter(VAF_tumour > vaf_cutoff) %>%
      dplyr::select(CHROM_A, START_A, END_A, manta_name)
  }

  #split manta_name variable
  svs_df = data.frame(svs$CHROM_A, svs$START_A, svs$END_A, do.call(rbind, strsplit(svs$manta_name, split = ":", fixed = TRUE)))

  #rename variables
  names(svs_df)[1] = "chrom"
  names(svs_df)[2] = "start"
  names(svs_df)[3] = "end"
  names(svs_df)[4] = "type"

  #subset df on SV type
  manta_sv = dplyr::filter(svs_df, type %in% c("MantaDEL", "MantaDUP")) %>%
    dplyr::select(chrom, start, end, type)

  #calculate sizes
  manta_sv$size = manta_sv$end - manta_sv$start
  manta_sv$type = as.factor(manta_sv$type)

  #add chr prefix, if missing
  if(!str_detect(manta_sv$chrom, "chr")[1]){
    manta_sv = mutate(manta_sv, chrom = paste0("chr", chrom))
  }

  #subset on selected chromosomes
  manta_sv = manta_sv[manta_sv$chrom %in% chr_select, ]

  #filter out varaints < 50 bp
  manta_sv = dplyr::filter(manta_sv, size >= size_cutoff)
  manta_sv$row_num = seq.int(nrow(manta_sv))

  del_col = get_gambl_colours("indels")[[1]]
  dup_col = get_gambl_colours("indels")[[2]]

  #plotting
  p = ggplot(manta_sv, aes(x = size, fill = type)) +
        geom_density(alpha = 0.7, color = NA, adjust = adjust_value, trim = trim) +
        labs(title = plot_title, subtitle = plot_subtitle, x = "Size (bp)", y = "") +
        scale_fill_manual(values = c(del_col, dup_col)) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        {if(hide_legend)theme(legend.position = "none")} +
        theme_cowplot()

  return(p)
}


#' Visualize (stacked barplot) genomic read-subsets across a selection of samples.
#'
#' @param these_samples Data frame with sample IDs (to be plotted) in the first column.
#' @param metadata Optional argument, used to derive sample IDs if sample_table is Null.
#' @param these_samples_metadata GAMBL metadata subset to the cases you want to process.
#' @param comparison_group Optional argument for plotting mean alignment metrics. Default is plotting the mean for samples provided. This parameter takes a list of sample IDs.
#' @param seq_type Subset qc metrics to a specific seq_type, default is genome.
#' @param add_mean Set to TRUE to superimpose mean values of plotted variables. Default is TRUE.
#' @param add_corrected_coverage Set to TRUE to add corrected coverage for selected samples.
#' @param keep_cohort If no df with sample IDs is supplied (these_samples = NULL) the function calls get_gambl_metadata and subsets on selected cohort.
#' @param keep_pathology If no df with sample IDs is supplied (these_samples = NULL) the function calls get_gambl_metadata and subsets on selected pathology.
#' @param this_color_palette Optional parameter that holds the selected colours for the plotted bars.
#' @param plot_sub Optional parameter, add a subtitle to alignment metric plot.
#'
#' @return plot as ggplot object.
#' @import tidyverse cowplot
#' @export
#'
#' @examples
#' #Example 1 - using these_samples parameter
#' #subset on FL cases with QC metrics available and plot
#' kridel_fl = get_gambl_metadata() %>%
#'  dplyr::filter(pathology == "FL", cohort == "FL_Kridel") %>%
#'  dplyr::select(sample_id) %>%
#'  pull(sample_id)
#'
#' my_plot_1 = fancy_alignment_plot(these_samples = kridel_fl, seq_type = "genome")
#'
#' #Example 2 - using already filtered metadata (these_samples_metadata)
#' fl_metadata = get_gambl_metadata() %>%
#'  dplyr::filter(pathology == "FL", cohort == "FL_Kridel")
#'
#' my_plot_2 = fancy_alignment_plot(these_samples_metadata = fl_metadata, seq_type = "genome")
#'
#' #Example 3 - using in-house metadata fitlering options
#' my_plot_3 = fancy_alignment_plot(keep_cohort = "FL_Kridel", keep_pathology = "FL", seq_type = "genome")
#'
fancy_alignment_plot = function(these_samples,
                                metadata,
                                these_samples_metadata,
                                comparison_group,
                                seq_type = "genome",
                                add_mean = TRUE,
                                add_corrected_coverage = TRUE,
                                keep_cohort,
                                keep_pathology,
                                this_color_palette = c("TotalReads" = "#3D405B",
                                                       "TotalUniquelyMapped" = "#81B29A",
                                                       "TotalDuplicatedreads" = "#E07A5F"),
                                plot_sub = ""){

  #get gambl metadata (if not supplied)
  if(missing(metadata)){
    this_meta = get_gambl_metadata(seq_type_filter = seq_type)
  }else{
    this_meta = metadata
  }

  if(!missing(these_samples_metadata)){
    these_samples = dplyr::select(these_samples_metadata, sample_id) %>%
      as.data.frame(strings.as.factors = FALSE) %>%
      pull(sample_id)
  }

  #filter metadata on selected cohort/pathology
  if(missing(these_samples)){
    if(!missing(keep_cohort) && missing(keep_pathology)){
      these_samples = dplyr::filter(this_meta, cohort == keep_cohort) %>%
        pull(sample_id)
    }

    if(!missing(keep_pathology) && missing(keep_cohort)){
      these_samples = dplyr::filter(this_meta, pathology == keep_pathology) %>%
        pull(sample_id)
    }

    if(!missing(keep_cohort) && !missing(keep_pathology)){
      these_samples = dplyr::filter(this_meta, pathology == keep_pathology, cohort == keep_cohort) %>%
        pull(sample_id)
    }

    if(missing(keep_cohort) && missing(keep_pathology)){
      these_samples = dplyr::select(this_meta, sample_id) %>%
        pull(sample_id)
    }
  }

  #get qc data for selected samples
  qc_metrics = collate_results(sample_table = these_samples, seq_type_filter = seq_type)

  message(paste0("QC Metric successfully retreived for ", nrow(qc_metrics),
                 " samples out of a total of ", nrow(these_samples), " samples in input sample table."))

  #subset alignment metrics
  melt_align = dplyr::select(qc_metrics, c(sample_id, TotalReads, TotalUniquelyMapped, TotalDuplicatedreads)) %>%
    melt(id.var = "sample_id") %>%
    arrange(sample_id)

  mean_cov_df = data.frame(Metric = c("TotalReads", "TotalUniquelyMapped", "TotalDuplicatedreads"),
                           Value = c (mean(melt_align$value[melt_align$variable == "TotalReads"]),
                                      mean(melt_align$value[melt_align$variable == "TotalUniquelyMapped"]),
                                      mean(melt_align$value[melt_align$variable == "TotalDuplicatedreads"])))

  if(!missing(comparison_group)){
    comp_data = collate_results(sample_table = comparison_group, seq_type_filter = seq_type) %>%
      dplyr::select(sample_id, TotalReads, TotalUniquelyMapped, TotalDuplicatedreads) %>%
      melt(id.var = "sample_id") %>%
      arrange(sample_id)

    mean_cov_df_comp = data.frame(Metric = c("TotalReads", "TotalUniquelyMapped", "TotalDuplicatedreads"),
                                  Value = c (mean(comp_data$value[comp_data$variable == "TotalReads"]),
                                             mean(comp_data$value[comp_data$variable == "TotalUniquelyMapped"]),
                                             mean(comp_data$value[comp_data$variable == "TotalDuplicatedreads"])))
  }

  #corrected mean coverage
  melt_cov = dplyr::select(qc_metrics, c(sample_id, MeanCorrectedCoverage)) %>%
    melt(id.var = "sample_id") %>%
    arrange(sample_id)

  #plot alignment data
  p = ggplot() +
    geom_bar(melt_align, mapping = aes(x = sample_id, y = value, fill = variable), position = "dodge", stat = "identity") +
    {if(add_mean)geom_hline(mean_cov_df, mapping = aes(yintercept = Value, color = Metric))} +
    {if(!missing(comparison_group)) geom_hline(mean_cov_df_comp, mapping = aes(yintercept = Value, linetype = Metric), color = "#54555E")} +
    {if(add_corrected_coverage)geom_point(melt_cov, mapping = aes(x = sample_id, y = value * 25000000, shape = variable), fill = "#A892B3", color = "#5C4966", size = 3, position = "dodge", stat = "identity")} +
    {if(add_corrected_coverage)scale_y_continuous(expand = c(0, 0), breaks = seq(0, max(melt_align$value), by = 3e+08), sec.axis = sec_axis(~ . / 2500000000, name = "", labels = function(b){paste0(round(b * 100, 0), "X")}))} +
    {if(!add_corrected_coverage)scale_y_continuous(expand = c(0, 0), breaks = seq(0, max(melt_align$value), by = 3e+08))} +
    labs(title = "Alignment Summary", subtitle = plot_sub, x = "", y = "Reads (n)") +
    scale_fill_manual(values = this_color_palette) +
    scale_linetype_manual(values = c("TotalReads" = "solid", "TotalUniquelyMapped" = "dashed", "TotalDuplicatedreads" = "dotted")) +
    scale_shape_manual(values = c("MeanCorrectedCoverage" = 21)) +
    scale_color_manual(values = c(this_color_palette)) +
    theme_cowplot() +
    labs(linetype = "Comparison Group", shape = "Corrected Coverage (right y-axis)", fill = "Alignment Metrics", color = "Alignment Metrics (Mean)") +
    theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0.5), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank())

  return(p)
}



#' Plot for visualizing QC metrics and allowing for grouping by different metadata columns.
#'
#' @param these_samples Data frame with sample IDs (to be plotted) in the first column (has to be named sample_id).
#' @param keep_cohort Optional parameter to be used when these_sample is NULL. Calls get_gambl_metadata() and filters on the cohort supplied in this parameter.
#' @param keep_pathology Optional parameter to be used when these_sample is NULL. Calls get_gambl_metadata() and filters on the pathology supplied in this parameter.
#' @param seq_type Selected seq type for incoming QC metrics.
#' @param metadata Optional, user can provide a metadata df to subset sample IDs from.
#' @param these_samples_metadata GAMBL metadata subset to the cases you want to process.
#' @param sort_by Plotting parameter, set sorting column for bar plots.
#' @param plot_data Plotting parameter, define the data type to be plotted.
#' @param fill_by Parameter for specifying fill variable for grouped bar plot. Can be any factor from incoming metadata, e.g pathology, cohort, etc.
#' @param labels If HTML plot version is rendered, you can specify what labels should be visible when hovering over the dots. Default is sample id and cohort. This parameter expects a vector of charachters.
#' @param interactive Boolean parameter for generating interactive plot (HTML). Default is FALSE.
#' @param comparison_samples Optional parameter, give the function a list of sample IDs to be compared against the main plotting group. Pathology is default.
#' @param plot_title Plotting parameter, plot title.
#' @param y_axis_lab Plotting parameter, label of y-axis.
#' @param return_plotdata Optional parameter, if set to TRUE a list of acceptable data types for plotting will be returned, and nothing else.
#'
#' @return plot as ggplot object.
#' @import tidyverse cowplot ggbeeswarm plotly
#' @export
#'
#' @examples
#' #Example 1 - using these_samples parameter
#' #subset on FL cases with QC metrics available and plot
#' kridel_fl = get_gambl_metadata() %>%
#'  dplyr::filter(pathology == "FL", cohort == "FL_Kridel") %>%
#'  dplyr::select(sample_id) %>%
#'
#' my_plot_1 = fancy_qc_plot(these_samples = kridel_fl,
#'                           seq_type = "genome",
#'                           interactive = FALSE,
#'                           plot_data = "AverageBaseQuality",
#'                           y_axis_lab = "Average Base Quality",
#'                           plot_title = "Average Base Quality For FL_Kridel")
#'
#' #Example 2 - using already filtered metadata (these_samples_metadata)
#' fl_metadata = get_gambl_metadata() %>%
#'  dplyr::filter(pathology == "FL", cohort == "FL_Kridel")
#'
#' my_plot_2 = fancy_qc_plot(these_samples_metadata = fl_metadata,
#'                           seq_type = "genome",
#'                           interactive = TRUE,
#'                           labels = c("cohort", "pathology")
#'                           plot_data = "AverageBaseQuality",
#'                           y_axis_lab = "Average Base Quality",
#'                           plot_title = "Average Base Quality For FL_Kridel")
#'
#' #Example 3 - using in-house metadata filtering options
#' my_plot_3 = fancy_qc_plot(keep_cohort = "FL_Kridel",
#'                           keep_pathology = "FL",
#'                           seq_type = "genome",
#'                           plot_data = "AverageBaseQuality",
#'                           y_axis_lab = "Average Base Quality",
#'                           plot_title = "Average Base Quality For FL_Kridel")
#'
fancy_qc_plot = function(these_samples,
                         keep_cohort,
                         keep_pathology,
                         seq_type = "genome",
                         metadata,
                         these_samples_metadata,
                         plot_data,
                         fill_by = "pathology",
                         labels = c("sample_id", "cohort"),
                         interactive = FALSE,
                         comparison_samples,
                         plot_title = "",
                         y_axis_lab = "",
                         return_plotdata = FALSE){

  #return a list of acceptable data types for plotting
  if(return_plotdata){
    plotting_variables = c("AverageBaseQuality", "AverageInsertSize", "AverageReadLength",
                           "PairsOnDiffCHR", "TotalReads", "TotalUniquelyMapped",
                           "TotalUnmappedreads", "TotalDuplicatedreads", "ProportionReadsDuplicated",
                           "ProportionReadsMapped", "MeanCorrectedCoverage", "ProportionCoverage10x", "ProportionCoverage30x")

    return(plotting_variables)
  }

  #get gambl metadata (if not supplied)
  if(missing(metadata)){
    this_meta = get_gambl_metadata(seq_type_filter = seq_type)
  }else{
    this_meta = metadata
  }

  if(!missing(these_samples_metadata)){
    these_samples = dplyr::select(these_samples_metadata, sample_id) %>%
      as.data.frame(strings.as.factors = FALSE) %>%
      pull(sample_id)
  }

  #filter metadata on selected cohort/pathology
  if(missing(these_samples) && missing(these_samples_metadata)){
    if(!missing(keep_cohort) && missing(keep_pathology)){
      these_samples = dplyr::filter(this_meta, cohort == keep_cohort) %>%
        pull(sample_id)
    }

    if(!missing(keep_pathology) && missing(keep_cohort)){
      these_samples = dplyr::filter(this_meta, pathology == keep_pathology) %>%
        pull(sample_id)
    }

    if(!missing(keep_cohort) && !missing(keep_pathology)){
      these_samples = dplyr::filter(this_meta, pathology == keep_pathology, cohort == keep_cohort) %>%
        pull(sample_id)
    }

    if(missing(keep_cohort) && missing(keep_pathology)){
      these_samples = dplyr::select(this_meta, sample_id) %>%
        pull(sample_id)
    }
  }

  #get QC data for selected samples
  qc_metrics = collate_results(sample_table = these_samples, seq_type_filter = seq_type)
  message(paste0("QC Metric successfully retreived for ", nrow(qc_metrics), " samples out of a total of ", nrow(these_samples), " samples in input sample table."))

  #aggregate sample list with metadata columns
  qc_meta = qc_metrics %>%
    inner_join(this_meta) %>%
    mutate_if(is.integer, as.factor) %>%
    mutate_if(is.character, as.factor)

  qc_meta$group = "main_sample"

  #Retrieve QC metrics for comparison samples, if provided.
  if(!missing(comparison_samples)){
    comp_data = collate_results(sample_table = comparison_samples, seq_type_filter = seq_type)

    #aggregate sample list with metadata columns
    comp_meta = comp_data %>% inner_join(this_meta)
    comp_meta = mutate_if(comp_meta, is.integer, as.factor)
    comp_meta = mutate_if(comp_meta, is.character, as.factor)
    comp_meta$group = "comparison_sample"
    qc_meta = rbind(qc_meta, comp_meta)
  }

  #get gambl colours for selected fill and subset to levels in selected factor
  col_gambl = get_gambl_colours(fill_by) %>%
    as.data.frame()

  col_gambl$factors = rownames(col_gambl)
  colnames(col_gambl)[1] = "hex"
  row.names(col_gambl) <- NULL

  levels_fill = levels(qc_meta[[fill_by]]) %>%
    as.data.frame()

  colnames(levels_fill)[1] = "factors"
  sub_cols = dplyr::left_join(levels_fill, col_gambl, by = "factors")
  list_col = as.list(sub_cols$hex)

  #plotting
  p = ggplot(qc_meta) +
    {if(interactive)aes_string(x = paste0("group"), y = plot_data, fill = fill_by, label1 = labels[1], label2 = labels[2])} +
    {if(!interactive)aes_string(x = paste0("group"), y = plot_data, fill = fill_by)} +
    geom_boxplot(mapping = aes(x = group), outlier.shape = NA) +
    geom_quasirandom() +
    labs(title = plot_title, x = "", y = y_axis_lab) +
    theme_cowplot() +
    scale_fill_manual(values = c(list_col)) +
    theme(legend.position = "right", axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank())

  #make plot interactive (html) with plotly.
  if(interactive){
    p = ggplotly(p)
  }
  return(p)
}


#' Visualize proportional coverage (10X and 30X) for selected samples and add comparison group (optional).
#'
#' @param these_samples Data frame with sample IDs (to be plotted) in the first column.
#' @param metadata Optional, user can provide a metadata df to subset sample IDs from.
#' @param these_samples_metadata GAMBL metadata subset to the cases you want to process.
#' @param keep_cohort Optional parameter to be used when these_sample is NULL. Calls get_gambl_metadata() and filters on the cohort supplied in this parameter.
#' @param keep_pathology Optional parameter to be used when these_sample is NULL. Calls get_gambl_metadata() and filters on the pathology supplied in this parameter.
#' @param comparison_samples Optional parameter, give the function a list of sample IDs to be compared against the main plotting group.
#' @param seq_type Selected seq type for incoming QC metrics.
#' @param plot_subtitle Plotting parameter, subtitle of generated plot.
#'
#' @return plot as ggplot object.
#' @import tidyverse cowplot
#' @export
#'
#' @examples
#' #Example 1 - using these_samples parameter
#' #subset on FL cases with QC metrics available and plot
#' kridel_fl = get_gambl_metadata() %>%
#'  dplyr::filter(pathology == "FL", cohort == "FL_Kridel") %>%#
#'  dplyr::select(sample_id) %>%
#'  pull(sample_id)
#'
#' my_plot_1 = fancy_propcov_plot(these_samples = kridel_fl, seq_type = "genome")
#'
#' #Example 2 - using already filtered metadata (these_samples_metadata)
#' fl_metadata = get_gambl_metadata() %>%
#'  dplyr::filter(pathology == "FL", cohort == "FL_Kridel")
#'
#' my_plot_2 = fancy_propcov_plot(these_samples_metadata = fl_metadata, seq_type = "genome")
#'
#' #Example 3 - using in-house metadata fitlering options
#' my_plot_3 = fancy_propcov_plot(keep_cohort = "FL_Kridel", keep_pathology = "FL", seq_type = "genome")
#'
fancy_propcov_plot = function(these_samples,
                              metadata,
                              these_samples_metadata,
                              keep_cohort,
                              keep_pathology,
                              comparison_samples,
                              seq_type = "genome",
                              plot_subtitle = ""){

  #get gambl metadata (if not supplied)
  if(missing(metadata)){
    this_meta = get_gambl_metadata(seq_type_filter = seq_type)
  }else{
    this_meta = metadata
  }

  if(!missing(these_samples_metadata)){
    these_samples = dplyr::select(these_samples_metadata, sample_id) %>%
      as.data.frame(strings.as.factors = FALSE) %>%
      pull(sample_id)
  }

  #filter metadata on selected cohort/pathology
  if(missing(these_samples) && missing(these_samples_metadata)){
    if(!missing(keep_cohort) && missing(keep_pathology)){
      these_samples = dplyr::filter(this_meta, cohort == keep_cohort) %>%
        pull(sample_id)
    }

    if(!missing(keep_pathology) && missing(keep_cohort)){
      these_samples = dplyr::filter(this_meta, pathology == keep_pathology) %>%
        pull(sample_id)
    }

    if(!missing(keep_cohort) && !missing(keep_pathology)){
      these_samples = dplyr::filter(this_meta, pathology == keep_pathology, cohort == keep_cohort) %>%
        pull(sample_id)
    }

    if(missing(keep_cohort) && missing(keep_pathology)){
      these_samples = dplyr::select(this_meta, sample_id) %>%
        pull(sample_id)
    }
  }

  #retrieve data for comparison, provided as a df with sample IDs in the first column (subset from gambl metadata)
  if(!missing(comparison_samples)){
    comp_data = collate_qc_results(sample_table = comparison_samples, seq_type_filter = seq_type) %>%
      dplyr::select(ProportionCoverage10x, ProportionCoverage30x) %>%
      gather(Type, Value)

    comp_data$Type = as.factor(comp_data$Type)

    levels(comp_data$Type)[levels(comp_data$Type)=="ProportionCoverage10x"] = "comparison_group_10X"
    levels(comp_data$Type)[levels(comp_data$Type)=="ProportionCoverage30x"] = "comparison_group_30X"
  }

  #get QC data for selected samples
  qc_metrics = collate_results(sample_table = these_samples, seq_type_filter = seq_type)
  message(paste0("QC Metric successfully retreived for ", nrow(qc_metrics), " samples out of a total of ", nrow(these_samples), " samples in input sample table."))

  #data wrangling steps
  sub_metrics = dplyr::select(qc_metrics, ProportionCoverage10x, ProportionCoverage30x) %>%
    gather(Type, Value)

  sub_metrics$Type = as.factor(sub_metrics$Type)

  levels(sub_metrics$Type)[levels(sub_metrics$Type)=="ProportionCoverage10x"] = "selected_samples_10X"
  levels(sub_metrics$Type)[levels(sub_metrics$Type)=="ProportionCoverage30x"] = "selected_samples_30X"

  #combine comparison data with sample data
  if(!missing(comparison_samples)){
    sub_metrics = rbind(sub_metrics, comp_data) %>%
      mutate(Type = factor(Type, levels = c("selected_samples_10X", "comparison_group_10X", "selected_samples_30X", "comparison_group_30X")))
  }

  #plotting
  p = ggplot(data = sub_metrics, aes(x = Type, y = Value, fill = Type)) +
    geom_violin(trim = FALSE, scale = "width") +
    stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2, fill = "white") +
    ylim(0, 1) +
    labs(title = "Proportion Coverage", subtitle = plot_subtitle, x = "", y = "Fraction") +
    theme_cowplot() +
    scale_fill_manual(values = c("#dda15e", "#606c38", "#433D6B", "#6B3254")) +
    theme(legend.position = "right", legend.title = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(), panel.background = element_blank())

  return(p)
}


#' Visualize proportional metrics for selected samples.
#'
#' @param these_samples Data frame with sample IDs (to be plotted) in the first column.
#' @param metadata Optional, user can provide a metadata df to subset sample IDs from.
#' @param these_samples_metadata GAMBL metadata subset to the cases you want to process.
#' @param keep_cohort Optional parameter to be used when these_sample is NULL. Calls get_gambl_metadata() and filters on the cohort supplied in this parameter.
#' @param keep_pathology Optional parameter to be used when these_sample is NULL. Calls get_gambl_metadata() and filters on the pathology supplied in this parameter.
#' @param seq_type Selected seq type for incoming QC metrics.
#' @param plot_subtitle Plotting parameter, subtitle of generated plot.
#'
#' @return plot as ggplot object.
#' @import tidyverse cowplot
#' @export
#'
#' @examples
#' #Example 1 - using these_samples parameter
#' #subset on FL cases with QC metrics available and plot
#' kridel_fl = get_gambl_metadata() %>%
#'  dplyr::filter(pathology == "FL", cohort == "FL_Kridel") %>%#
#'  dplyr::select(sample_id) %>%
#'  pull(sample_id)
#'
#' my_plot_1 = fancy_proportions_plot(these_samples = kridel_fl, seq_type = "genome")
#'
#' #Example 2 - using already filtered metadata (these_samples_metadata)
#' fl_metadata = get_gambl_metadata() %>%
#'  dplyr::filter(pathology == "FL", cohort == "FL_Kridel")
#'
#' my_plot_2 = fancy_proportions_plot(these_samples_metadata = fl_metadata, seq_type = "genome")
#'
#' #Example 3 - using in-house metadata fitlering options
#' my_plot_3 = fancy_proportions_plot(keep_cohort = "FL_Kridel", keep_pathology = "FL", seq_type = "genome")
#'
fancy_proportions_plot = function(these_samples,
                                  metadata,
                                  these_samples_metadata,
                                  keep_cohort,
                                  keep_pathology,
                                  seq_type = "genome",
                                  plot_subtitle = ""){

  #get gambl metadata (if not supplied)
  if(missing(metadata)){
    this_meta = get_gambl_metadata(seq_type_filter = seq_type)
  }else{
    this_meta = metadata
  }

  if(!missing(these_samples_metadata)){
    these_samples = dplyr::select(these_samples_metadata, sample_id) %>%
      as.data.frame(strings.as.factors = FALSE) %>%
      pull(sample_id)
  }

  #filter metadata on selected cohort/pathology
  if(missing(these_samples) && missing(these_samples_metadata)){
    if(!missing(keep_cohort) && missing(keep_pathology)){
      these_samples = dplyr::filter(this_meta, cohort == keep_cohort) %>%
        pull(sample_id)
    }

    if(!missing(keep_pathology) && missing(keep_cohort)){
      these_samples = dplyr::filter(this_meta, pathology == keep_pathology) %>%
        pull(sample_id)
    }

    if(!missing(keep_cohort) && !missing(keep_pathology)){
      these_samples = dplyr::filter(this_meta, pathology == keep_pathology, cohort == keep_cohort) %>%
        pull(sample_id)
    }

    if(missing(keep_cohort) && missing(keep_pathology)){
      these_samples = dplyr::select(this_meta, sample_id) %>%
        pull(sample_id)
    }
  }

  #get QC data for selected samples
  qc_metrics = collate_results(sample_table = these_samples, seq_type_filter = seq_type)
  message(paste0("QC Metric successfully retreived for ", nrow(qc_metrics), " samples out of a total of ", nrow(these_samples), " samples in input sample table."))

  #data wrangling
  qc_sub = dplyr::select(qc_metrics, sample_id, ProportionReadsDuplicated, ProportionReadsMapped, ProportionCoverage10x, ProportionCoverage30x) %>%
    gather(Type, Value, -sample_id)

  qc_sub$Type = as.factor(qc_sub$Type)

  levels(qc_sub$Type)[levels(qc_sub$Type)=="ProportionReadsDuplicated"] = "Duplicated Reads"
  levels(qc_sub$Type)[levels(qc_sub$Type)=="ProportionReadsMapped"] = "Mapped Reads"
  levels(qc_sub$Type)[levels(qc_sub$Type)=="ProportionCoverage10x"] = "10X"
  levels(qc_sub$Type)[levels(qc_sub$Type)=="ProportionCoverage30x"] = "30X"

  #get means for each metric
  mean_df = data.frame(Metric = c("Duplicated Reads", "Mapped Reads", "10X", "30X"),
                       Value = c(mean(qc_sub$Value[qc_sub$Type == "Duplicated Reads"]),
                                 mean(qc_sub$Value[qc_sub$Type == "Mapped Reads"]),
                                 mean(qc_sub$Value[qc_sub$Type == "10X"]),
                                 mean(qc_sub$Value[qc_sub$Type == "30X"])))

  #set colors
  these_colors = c("Duplicated Reads" = "#B8794D",
                   "Mapped Reads" = "#366B32",
                   "10X" = "#DDA15E",
                   "30X" = "#433D6B")

  #plotting
  p = ggplot(qc_sub, aes(x = sample_id, y = Value, group = Type)) +
    geom_bar(aes(fill = Type), position = "dodge", stat = "identity", width = 0.7) +
    geom_hline(mean_df, mapping = aes(yintercept = Value, color = Metric), size = 0.5, linetype = "dashed") +
    labs(title = "QC Metrics As Proportions (of All Reads)", subtitle = plot_subtitle, x = "", y = "Proportion") +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, by = 0.1)) +
    scale_fill_manual(values = these_colors) +
    scale_color_manual(values = these_colors) +
    labs(fill = "Metrics", color = "Mean") +
    theme(legend.position = "right", axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5),
          axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.background = element_blank(),
          axis.line.x = element_line(color = "black", size = 0.5, linetype = 1),
          axis.line.y = element_line(color = "black", size = 0.5, linetype = 1),
          axis.ticks.y = element_line(color = "black", size = 0.5, linetype = 1))

  return(p)
}
