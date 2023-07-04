#' Oncogenes in grch37 genome build.
#'
#' A data frame with the coordinates of lymphoma oncogenes relative to the grch37 genome build.
#'
#' @format ## `grch37_oncogene`
#' A data frame with 19 rows and 5 columns.
#' \describe{
#'   \item{chrom}{Chromosomes without chr-prefix, 1:22.}
#'   \item{start}{Start coordinate for the specified oncogene.}
#'   \item{end}{End coordinate for the specified oncogene.}
#'   \item{gene}{Lymphoma oncogene.}
#'   \item{entrez}{ENTREZ ID for the specified oncogene.}
#' }
"grch37_oncogene"

#' Oncogenes in hg38 genome build.
#'
#' A data frame with the coordinates of lymphoma oncogenes relative to the hg38 genome build.
#'
#' @format ## `hg38_oncogene`
#' A data frame with 19 rows and 5 columns.
#' \describe{
#'   \item{chrom}{Chromosomes without chr-prefix, 1:22.}
#'   \item{start}{Start coordinate for the specified oncogene.}
#'   \item{end}{End coordinate for the specified oncogene.}
#'   \item{gene}{Lymphoma oncogene.}
#'   \item{entrez}{ENTREZ ID for the specified oncogene.}
#' }
"hg38_oncogene"

#' Chromosome Arms grch37.
#'
#' A data frame with the chromosome arm coordinates in respect to grch37.
#'
#' @format ## `chromosome_arms_grch37`
#' A data frame with 48 rows and 4 columns.
#' \describe{
#'   \item{chromosome}{Chromosomes without chr-prefix, 1:22, X and Y.}
#'   \item{start}{Start coordinates for the specified chromosome arm.}
#'   \item{end}{End coordinates for the specified chromosome arm.}
#'   \item{arm}{Chromosome arm, either p or q.}
#' }
"chromosome_arms_grch37"


#' Chromosome Arms hg38.
#'
#' A data frame with the chromosome arm coordinates in respect to hg38.
#'
#' @format ## `chromosome_arms_hg38`
#' A data frame with 48 rows and 4 columns.
#' \describe{
#'   \item{chromosome}{Chromosomes with chr-prefix, 1:22, X and Y.}
#'   \item{start}{Start coordinates for the specified chromosome arm.}
#'   \item{end}{End coordinates for the specified chromosome arm.}
#'   \item{arm}{Chromosome arm, either p or q.}
#' }
"chromosome_arms_hg38"


#' Double Hit Signature Genes With Weights.
#'
#' A data frame with double hit signature genes (both as ensembl IDs and Hugo symbols) and importance scores.
#'
#' @format ## `dhitsig_genes_with_weights`
#' A data frame with 104 rows and 3 columns.
#' \describe{
#'   \item{Ensembl_ID}{Ensembl IDs as factors, 104 levels.}
#'   \item{ImportanceScore}{Numeric column with importance scores.}
#'   \item{Hugo_Symbol}{Gene symbols in Hugo format as a factor with 104 levels.}
#' }
"dhitsig_genes_with_weights"


#' Genes Blacklist.
#'
#' A tibble with gene symbols (Hugo) that falls within blacklisted regions of the genome.
#'
#' @format ## `gene_blacklist`
#' A tibble with 291 rows.
#' \describe{
#'   \item{Gene}{Genes symbols in Hugo format.}
#' }
"gene_blacklist"


#' Grande et al. MAF.
#'
#' A MAF (data frame) drawn from the Grande et al. dataset.
#'
#' @format ## `grande_maf`
#' A MAF in data frame format. 12251 rows and 125 columns.
#' \describe{
#'   \item{Hugo_Symbol}{HUGO symbol for the gene (HUGO symbols are always in all caps). "Unknown" is used for regions that do not correspond to a gene}
#'   \item{Entrez_Gene_ID}{Entrez gene ID (an integer). "0" is used for regions that do not correspond to a gene region or Ensembl ID}
#'   \item{Centre}{One or more genome sequencing center reporting the variant}
#'   \item{NCBI_Build}{	The reference genome used for the alignment}
#'   \item{Chromosome}{The affected chromosome}
#'   \item{Start_Position}{Lowest numeric position of the reported variant on the genomic reference sequence. Mutation start coordinate}
#'   \item{End_Position}{Highest numeric genomic position of the reported variant on the genomic reference sequence. Mutation end coordinate}
#'   \item{Strand}{Genomic strand of the reported allele. Currently, all variants will report the positive strand: '+'}
#'   \item{Variant_Classification}{Translational effect of variant allele}
#'   \item{Variant_Type}{Type of mutation. TNP (tri-nucleotide polymorphism) is analogous to DNP (di-nucleotide polymorphism) but for three consecutive nucleotides. ONP (oligo-nucleotide polymorphism) is analogous to TNP but for consecutive runs of four or more (SNP, DNP, TNP, ONP, INS, DEL, or Consolidated)}
#'   \item{Reference_Allele}{The plus strand reference allele at this position. Includes the deleted sequence for a deletion or "-" for an insertion}
#'   \item{Tumor_Seq_Allele1}{Primary data genotype for tumor sequencing (discovery) allele 1. A "-" symbol for a deletion represents a variant. A "-" symbol for an insertion represents wild-type allele. Novel inserted sequence for insertion does not include flanking reference bases}
#'   \item{Tumor_Seq_Allele2}{Tumor sequencing (discovery) allele 2}
#'   \item{dbSNP_RS}{The rs-IDs from the dbSNP database, "novel" if not found in any database used, or null if there is no dbSNP record, but it is found in other databases}
#'   \item{dbSNP_Val_Status}{The dbSNP validation status is reported as a semicolon-separated list of statuses. The union of all rs-IDs is taken when there are multiple}
#'   \item{Tumor_Sample_Barcode}{Aliquot barcode for the tumor sample}
#'   \item{Matched_Norm_Sample_Barcode}{Aliquot barcode for the matched normal sample}
#'   \item{Match_Norm_Seq_Allele1}{Primary data genotype. Matched normal sequencing allele 1. A "-" symbol for a deletion represents a variant. A "-" symbol for an insertion represents wild-type allele. Novel inserted sequence for insertion does not include flanking reference bases (cleared in somatic MAF)}
#'   \item{Match_Norm_Seq_Allele2}{Matched normal sequencing allele 2}
#'   \item{Tumor_Validation_Allele1}{Secondary data from orthogonal technology. Tumor genotyping (validation) for allele 1. A "-" symbol for a deletion represents a variant. A "-" symbol for an insertion represents wild-type allele. Novel inserted sequence for insertion does not include flanking reference bases}
#'   \item{Tumor_Validation_Allele2}{	Secondary data from orthogonal technology. Tumor genotyping (validation) for allele 2}
#'   \item{Verification_Status}{Second pass results from independent attempt using same methods as primary data source. Generally reserved for 3730 Sanger Sequencing}
#'   \item{Validation_Status}{Second pass results from orthogonal technology}
#'   \item{Mutation_Status}{An assessment of the mutation as somatic, germline, LOH, post transcriptional modification, unknown, or none. The values allowed in this field are constrained by the value in the Validation_Status field}
#'   \item{Sequencing_Phase}{TCGA sequencing phase (if applicable). Phase should change under any circumstance that the targets under consideration change}
#'   \item{Sequence_Source}{Molecular assay type used to produce the analytes used for sequencing. Allowed values are a subset of the SRA 1.5 library_strategy field values. This subset matches those used at CGHub}
#'   \item{Validation_Method}{The assay platforms used for the validation call}
#'   \item{Score}{Boolean variable}
#'   \item{BAM_File}{Boolean column stating if BAM file exists or not}
#'   \item{Sequencer}{Instrument used to produce primary sequence data}
#'   \item{Tumor_Sample_UUID}{GDC aliquot UUID for tumor sample}
#'   \item{Matched_Norm_Sample_UUID}{GDC aliquot UUID for matched normal sample}
#'   \item{HGVSc}{The coding sequence of the variant in HGVS recommended format}
#'   \item{HGVSp}{The protein sequence of the variant in HGVS recommended format. "p.=" signifies no change in the protein}
#'   \item{HGVS_Short}{Same as the HGVSp column, but using 1-letter amino-acid codes}
#'   \item{Transcript_ID}{Ensembl ID of the transcript affected by the varian}
#'   \item{Exon_Number}{The exon number (out of total number)}
#'   \item{t_depth}{Read depth across this locus in tumor BAM}
#'   \item{t_ref_count}{Read depth supporting the reference allele in tumor BAM}
#'   \item{t_alt_count}{Read depth supporting the variant allele in tumor BAM}
#'   \item{n_depth}{Read depth across this locus in normal BAM}
#'   \item{n_ref_count}{Read depth supporting the reference allele in normal BAM (cleared in somatic MAF)}
#'   \item{n_alt_count}{Read depth supporting the variant allele in normal BAM (cleared in somatic MAF)}
#'   \item{all_effects}{A semicolon delimited list of all possible variant effects, sorted by priority}
#'   \item{Allele}{The variant allele used to calculate the consequence}
#'   \item{Gene}{Stable Ensembl ID of affected gene}
#'   \item{Feature}{Stable Ensembl ID of feature (transcript, regulatory, motif)}
#'   \item{Feature_type}{Type of feature. Currently one of Transcript, RegulatoryFeature, MotifFeature (or blank)}
#'   \item{Consequence}{Consequence type of this variant; sequence ontology terms}
#'   \item{cDNA_position}{Relative position of base pair in the cDNA sequence as a fraction. A "-" symbol is displayed as the numerator if the variant does not appear in cDNA}
#'   \item{CDS_position}{Relative position of base pair in coding sequence. A "-" symbol is displayed as the numerator if the variant does not appear in coding sequence}
#'   \item{Protein_position}{Relative position of affected amino acid in protein. A "-" symbol is displayed as the numerator if the variant does not appear in coding sequence}
#'   \item{Amino_acids}{Only given if the variation affects the protein-coding sequence}
#'   \item{Codons}{The alternative codons with the variant base in upper case}
#'   \item{Existing_variation}{Known identifier of existing variation}
#'   \item{ALLELE_NUM}{Allele number from input; 0 is reference, 1 is first alternate etc.}
#'   \item{DISTANCE}{Shortest distance from the variant to transcript}
#'   \item{STRAND_VEP}{}
#'   \item{SYMBOL}{The gene symbol}
#'   \item{SYMBOL_SOURCE}{The source of the gene symbol}
#'   \item{HGNC_ID}{Gene identifier from the HUGO Gene Nomenclature Committee if applicable}
#'   \item{BIOTYPE}{Biotype of transcript}
#'   \item{CANONICAL}{A flag (YES) indicating that the VEP-based canonical transcript, the longest translation, was used for this gene. If not, the value is null}
#'   \item{CCDS}{The CCDS identifier for this transcript, where applicable}
#'   \item{ENSP}{The Ensembl protein identifier of the affected transcript}
#'   \item{SWISSPROT}{UniProtKB/Swiss-Prot accession}
#'   \item{TREMBL}{UniProtKB/TrEMBL identifier of protein product}
#'   \item{UNIPARC}{UniParc identifier of protein product}
#'   \item{RefSeq}{RefSeq identifier for this transcript}
#'   \item{SIFT}{The SIFT prediction and/or score, with both given as prediction (score)}
#'   \item{PolyPhen}{The PolyPhen prediction and/or score}
#'   \item{EXON}{The exon number (out of total number)}
#'   \item{INTRON}{The intron number (out of total number)}
#'   \item{DOMAINS}{The source and identifier of any overlapping protein domains}
#'   \item{GMAF}{Non-reference allele and frequency of existing variant in 1000 Genomes}
#'   \item{GMAF_Allele}{Non-reference allele and frequency of existing variant in 1000 Genomes}
#'   \item{GMAF_AF}{Non-reference allele and frequency of existing variant in 1000 Genomes}
#'   \item{AFR_MAF}{Non-reference allele and frequency of existing variant in 1000 Genomes combined African population}
#'   \item{AMR_MAF}{Non-reference allele and frequency of existing variant in 1000 Genomes combined American population}
#'   \item{ASN_MAF}{Non-reference allele and frequency of existing variant in 1000 Genomes combined Asian population}
#'   \item{EAS_MAF}{Non-reference allele and frequency of existing variant in 1000 Genomes combined East Asian population}
#'   \item{EUR_MAF}{Non-reference allele and frequency of existing variant in 1000 Genomes combined European population}
#'   \item{SAS_MAF}{Non-reference allele and frequency of existing variant in 1000 Genomes combined South Asian population}
#'   \item{AA_MAF}{Non-reference allele and frequency of existing variant in NHLBI-ESP African American population}
#'   \item{EA_MAF}{Non-reference allele and frequency of existing variant in NHLBI-ESP European American population}
#'   \item{CLIN_SIG}{Clinical significance of variant from dbSNP as annotated in ClinVar}
#'   \item{SOMATIC}{Somatic status of each ID reported under Existing_variation (0, 1, or null)}
#'   \item{PUBMED}{Pubmed ID(s) of publications that cite existing variant}
#'   \item{MOTIF_NAME}{The source and identifier of a transcription factor binding profile aligned at this position}
#'   \item{MOTIF_POS}{The relative position of the variation in the aligned TFBP}
#'   \item{HIGH_INF_POS}{A flag indicating if the variant falls in a high information position of a transcription factor binding profile (TFBP) (Y, N, or null)}
#'   \item{MOTIF_SCORE_CHANGE}{The difference in motif score of the reference and variant sequences for the TFBP}
#'   \item{IMAPCT}{The impact modifier for the consequence type}
#'   \item{PICK}{Indicates if this block of consequence data was picked by VEP's pick feature (1 or null)}
#'   \item{VARIANT_CLASS}{Sequence Ontology variant class}
#'   \item{TSL}{Transcript support level, which is based on independent RNA analyses}
#'   \item{HGVS_OFFSET}{Indicates by how many bases the HGVS notations for this variant have been shifted}
#'   \item{PHENO}{Indicates if existing variant is associated with a phenotype, disease or trait (0, 1, or null)}
#'   \item{MINIMISED}{Alleles in this variant have been converted to minimal representation before consequence calculation (1 or null)}
#'   \item{ExAC_AF}{Global Allele Frequency from ExAC}
#'   \item{ExAC_AF_AFR}{African/African American Allele Frequency from ExAC}
#'   \item{ExAC_AF_AMR}{American Allele Frequency from ExAC}
#'   \item{ExAC_AF_EAS}{East Asian Allele Frequency from ExAC}
#'   \item{ExAC_AF_FIN}{Finnish Allele Frequency from ExAC}
#'   \item{ExAC_AF_NFE}{Non-Finnish European Allele Frequency from ExAC}
#'   \item{ExAC_AF_OTH}{Other Allele Frequency from ExAC}
#'   \item{ExAC_AF_SAS}{South Asian Allele Frequency from ExAC}
#'   \item{GENE_PHENO}{Indicates if gene that the variant maps to is associated with a phenotype, disease or trait (0, 1, or null)}
#'   \item{FILTER}{Copied from input VCF. This includes filters implemented directly by the variant caller and other external software used in the DNA-Seq pipeline. See below for additional details.}
#'   \item{flanking_bps}{The flanking basepairs}
#'   \item{variant_id}{Variant ID}
#'   \item{variant_qual}{Variant quality}
#'   \item{ExAC_AF_Adj}{Adjusted Global Allele Frequency from ExAC}
#'   \item{ExAC_AC_AN_Adj}{Adjusted Global Allele Frequency from ExAC}
#'   \item{ExAC_AC_AN}{Global Allele Frequency from ExAC}
#'   \item{ExAC_AC_AN_AFR}{African/African American Allele Frequency from ExAC}
#'   \item{ExAC_AC_AN_AMR}{American Allele Frequency from ExAC}
#'   \item{ExAC_AC_AN_EAS}{East Asian Allele Frequency from ExAC}
#'   \item{ExAC_AC_AN_FIN}{Finnish Allele Frequency from ExAC}
#'   \item{ExAC_AC_AN_NFE}{Non-Finnish European Allele Frequency from ExAC}
#'   \item{ExAC_AC_AN_OTH}{Other Allele Frequency from ExAC}
#'   \item{ExAC_AC_AN_SAS}{South Asian Allele Frequency from ExAC}
#'   \item{ExAC_FILTER}{Filter}
#' }
"grande_maf"


#' Gene Coordinates for grch37 (all).
#'
#' All gene coordinates in respect to grch37.
#'
#' @format ## `grch37_all_gene_coordinates`
#' A data frame with 2602 rows and 6 columns.
#' \describe{
#'   \item{ensembl_gene_id}{Ensembl gene ID}
#'   \item{chromosome}{The chromosome that the gene is residing on}
#'   \item{start}{The start coordinates for the gene}
#'   \item{end}{The end coordinates for the gene}
#'   \item{gene_name}{The gene name}
#'   \item{hugo_symbol}{Gene symbol in Hugo format}
#' }
"grch37_all_gene_coordinates"


#' grch37 Gene Coordinates.
#'
#' All gene coordinates in respect to grch37.
#'
#' @format ## `grch37_gene_coordinates`
#' A data frame with 63,763 rows and 6 columns.
#' \describe{
#'   \item{ensembl_gene_id}{Ensembl gene ID}
#'   \item{chromosome}{The chromosome that the gene is residing on}
#'   \item{start}{The start coordinates for the gene}
#'   \item{end}{The end coordinates for the gene}
#'   \item{gene_name}{The gene name}
#'   \item{hugo_symbol}{Gene symbol in Hugo format}
#' }
"grch37_gene_coordinates"


#' Lymphoma Genes (grch37).
#'
#' Lymphoma associated genes in respect to grch37.
#'
#' @format ## `grch37_lymphoma_genes_bed`
#' A data frame with 195 rows and 4 columns.
#' \describe{
#'   \item{chromosome_name}{The chromosome for which the gene is residing on}
#'   \item{start_position}{The start coordinate for the gene}
#'   \item{end_position}{The end coordinate for the gene}
#'   \item{hgnc_symbol}{Gene symbol in Hugo format}
#' }
"grch37_lymphoma_genes_bed"


#' grch37 Partner Genes.
#'
#' Translocation partners for oncogenes in with coordinates in respect to grch37.
#'
#' @format ## `grch37_partners`
#' A data frame with 31 rows and 5 columns.
#' \describe{
#'   \item{chrom}{The chromosome for which the gene is residing on}
#'   \item{start}{The start coordinate for the gene}
#'   \item{end}{The end coordinate for the gene}
#'   \item{gene}{Gene symbol in Hugo format}
#'   \item{entrez}{Entrez ID}
#' }
"grch37_partners"

#' hg38 Gene Coordinates.
#'
#' All gene coordinates in respect to hg38.
#'
#' @format ## `hg38_gene_coordinates`
#' A data frame with 63,763 rows and 6 columns.
#' \describe{
#'   \item{ensembl_gene_id}{Ensembl gene ID}
#'   \item{chromosome}{The chromosome that the gene is residing on}
#'   \item{start}{The start coordinates for the gene}
#'   \item{end}{The end coordinates for the gene}
#'   \item{gene_name}{The gene name}
#'   \item{hugo_symbol}{Gene symbol in Hugo format}
#' }
"hg38_gene_coordinates"


#' Lymphoma Genes (hg38).
#'
#' Lymphoma associated genes in respect to hg38.
#'
#' @format ## `hg38_lymphoma_genes_bed`
#' A data frame with 195 rows and 4 columns.
#' \describe{
#'   \item{chromosome_name}{The chromosome for which the gene is residing on}
#'   \item{start_position}{The start coordinate for the gene}
#'   \item{end_position}{The end coordinate for the gene}
#'   \item{hgnc_symbol}{Gene symbol in Hugo format}
#' }
"hg38_lymphoma_genes_bed"


#' hg38 Partner Genes.
#'
#' Translocation partners for oncogenes in with coordinates in respect to hg38.
#'
#' @format ## `hg38_partners`
#' A data frame with 31 rows and 5 columns.
#' \describe{
#'   \item{chrom}{The chromsome for which the gene is residing on}
#'   \item{start}{The start coordinate for the gene}
#'   \item{end}{The end coordinate for the gene}
#'   \item{gene}{Gene symbol in Hugo format}
#'   \item{entrez}{Entrez ID}
#' }
"hg38_partners"


#' grch37 Hotspot Regions.
#'
#' Mutation hotspot regions in respect to grch37.
#'
#' @format ## `hotspot_regions_grch37`
#' A data frame with 6 rows and 3 columns.
#' \describe{
#'   \item{chrom}{Chromosome for the described region}
#'   \item{start}{The start coordinate for the region}
#'   \item{end}{The end coordinate for the region}
#' }
"hotspot_regions_grch37"


#' hg38 Hotspot Regions.
#'
#' Mutation hotspot regions in respect to hg38.
#'
#' @format ## `hotspot_regions_hg38`
#' A data frame with 6 rows and 3 columns.
#' \describe{
#'   \item{chrom}{Chromosome for the described region}
#'   \item{start}{The start coordinate for the region}
#'   \item{end}{The end coordinate for the region}
#' }
"hotspot_regions_hg38"


#' Lymphoma Genes Comprehensive.
#'
#' A detailed data frame with lymphoma genes, annotated with evidence from literature and aSHM.
#'
#' @format ## `lymphoma_genes_comprehensive`
#' A data frame with 127 rows and 9 columns.
#' \describe{
#'   \item{ensemble_gene_id}{Gene in Ensemble format}
#'   \item{Gene}{Gene symbol in Hugo format.}
#'   \item{Chapuy}{Boolean flag, TRUE if gene verified by the stated study (Chapuy)}
#'   \item{Reddy}{Boolean flag, TRUE if gene verified by the stated study (Reddy)}
#'   \item{LymphGen}{Boolean flag, TRUE if lymphGen}
#'   \item{curated}{Boolean flag, describing if the gene has been curated or not}
#'   \item{other}{Other support for the described gene}
#'   \item{Lacy}{Boolean flag, TRUE if gene verified by the stated study (Lacy)}
#'   \item{aSHM}{Boolean flag for annotating aSHM}
#' }
"lymphoma_genes_comprehensive"


#' Reddy Genes.
#'
#' Genes identified as significantly mutated in DLBCL by the study of Reddy et al.
#'
#' @format ## `reddy_genes`
#' A data frame with 150 rows and 4 columns.
#' \describe{
#'   \item{Input}{Input}
#'   \item{hgnc_symbol}{The HGNC symbol}
#'   \item{Approved name}{Approved name}
#'   \item{HGNC ID}{HGNC ID}
#' }
"reddy_genes"


#' Target Regions grch37.
#'
#' Target regions in respect to grch37.
#'
#' @format ## `target_regions_grch37`
#' A data frame with 295994 rows and 3 columns.
#' \describe{
#'   \item{chrom}{Chromosome for the descriebd region}
#'   \item{start}{Start coordinate of the region}
#'   \item{end}{End coordiante of the region}
#' }
"target_regions_grch37"


#' Target Regions hg38.
#'
#' Target regions in respect to hg38.
#'
#' @format ## `target_regions_hg38`
#' A data frame with 296453 rows and 3 columns.
#' \describe{
#'   \item{chrom}{Chromosome for the descriebd region}
#'   \item{start}{Start coordinate of the region}
#'   \item{end}{End coordiante of the region}
#' }
"target_regions_hg38"


#' Wright Genes With Weights.
#'
#' Description.
#'
#' @format ## `wright_genes_with_weights`
#' A data frame with 210 rows and 3 columns.
#' \describe{
#'   \item{Ensembl_ID}{Gene in Ensembl ID format}
#'   \item{Hugo_Symbol}{Gene symbol in Hugo format}
#'   \item{Weight_tValue}{Weight Value for the specified gene}
#' }
"wright_genes_with_weights"


#' Default mapping table between mutation type (aka, variant classification) to mutation class
#'
#' A dataset containing the mapping table between genomic mutation type (aka, variant classification) to mutation class.
#' This dataset comes from the g3viz package and was obtained via this URL:
#' https://github.com/morinlab/g3viz/tree/master/data
#'
#' @format A data frame with three columns:
#' \describe{
#'   \item{Mutation_Type}{Mutation type, aka, variant classification}
#'   \item{Mutation_Class}{mutation class}
#'   \item{Short_Name}{short name of mutation type}
#' }
#' @examples
#' mutation.table.df
"mutation.table.df"

#' Mapping table between gene.symbol, uniprot.id, and pfam
#'
#' A dataset containing the mapping table between Hugo symbol, UniProt ID, and
#' Pfam ACC. This dataset comes from the g3viz package and was obtained via this URL:
#' https://github.com/morinlab/g3viz/tree/master/data
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{symbol}{Gene symbol}
#'   \item{uniprot}{UniProt ID}
#'   \item{length}{protein length}
#'   \item{start}{starting position of Pfam domain}
#'   \item{end}{ending position of Pfam domain}
#'   \item{hmm.acc}{Pfam accession number}
#'   \item{hmm.name}{Pfam name}
#'   \item{type}{Pfam type, i.e., domain/family/motif/repeat/disordered/coiled-coil}
#' }
#' @examples
#' hgnc2pfam.df
#' @source Pfam (v31.0) and UniProt
"hgnc2pfam.df"
