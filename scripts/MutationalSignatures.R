library(MutationalPatterns)
library(BSgenome)
library(GenomicRanges)
library(dplyr)
library(rempsyc)

##### PREPROCESSING #####
#set the reference genome assembly for the study
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

#The datasets are previously filtered and saved as VCF files
#Import VCF files and set the annotations
vcf_files <- c("C:\\Users\\fazel\\...\\mutSigns\\AML_all.vcf",
               "C:\\Users\\fazel\\...\\CH_all.vcf",
               "C:\\Users\\fazel\\...\\AML_DNMT3A.vcf",
               "C:\\Users\\fazel\\...\\CH_DNMT3A.vcf",
               "C:\\Users\\fazel\\...\\AML_TET2.vcf",
               "C:\\Users\\fazel\\...\\CH_TET2.vcf")
sample_names <- c("AML_all","CH_all",
                  "AML_DNMT3A","CH_DNMT3A",
                  "AML_TET2","CH_TET2")
genes <- c("ALL","ALL",
           "DNMT3A","DNMT3A",
           "TET2","TET2")

#convert VCF files into IRange format
grl <- read_vcfs_as_granges(vcf_files,sample_names, ref_genome,change_seqnames = F, group = 'none', type = "snv")
chromosomes = paste0('chr', c(1:22,'X'))
seqlevels(grl, pruning.mode = 'tidy') = chromosomes

#A function for getting overall information from the IRange dataset 
# Define the function
process_mutations <- function(grl, ref_genome) {
  # Initialize lists to store results
  all_muts <- list()
  all_types <- list()
  all_contexts <- list()
  all_type_contexts <- list()
  # Loop through each element in the grl list
  for (i in seq_along(grl)) {
    all_muts[[i]] <- mutations_from_vcf(grl[[i]])
    all_types[[i]] <- mut_type(grl[[i]])
    all_contexts[[i]] <- mut_context(grl[[i]], ref_genome)
    all_type_contexts[[i]] <- type_context(grl[[i]], ref_genome)
  }
  return(list(all_muts = all_muts, all_types = all_types, all_contexts = all_contexts, all_type_contexts = all_type_contexts))
}

#getting the results
results <- process_mutations(grl, ref_genome)


##### PROCESSING ######
#mutational signature spectrum through all VCFs 
type_occurrences <- mut_type_occurrences(grl, ref_genome)
plot_spectrum(type_occurrences, CT = TRUE, 
                    indv_points = TRUE, legend = FALSE) #results in Figure 12 of the thesis

#signatures contribution to the reference signatures
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
plot_96_profile(mut_mat)
mut_mat_ext_context <- mut_matrix(grl, ref_genome, extension = 2)
head(mut_mat_ext_context)
plot_profile_heatmap(mut_mat_ext_context) #results in Figure 13 of the thesis

#deconvolution of finding contribution to the reference signatures for each gene
signatures = get_known_signatures()
fit_res <- fit_to_signatures(mut_mat, signatures)
plot_contribution(fit_res$contribution,
  coord_flip = FALSE,
  mode = "relative"
)
#Find the linear combination of mutation signatures that most closely reconstructs 
#the mutation matrix by solving the non-negative least-squares constraints problem.
sign_deconv<- fit_res$contribution

#find the top 6 mutational signatures contributing to our mutation data sets 
# Define the function to store the processed data frames
process_deconv <- function(sign_deconv) {
  # Define the column names and labels
  column_labels <- c("AML_all", "CH_all", "AML_DNMT3A", "CH_DNMT3A", "AML_TET2", "CH_TET2")
  
  processed_deconv <- list()
  
  for (i in seq_along(column_labels)) {
    df <- data.frame(sign_deconv[, i])
    colnames(df) <- column_labels[i]
    df$SBS <- rownames(df)
    
    df <- df %>% arrange(desc(!!sym(column_labels[i])))
    df <- df[!df[[column_labels[i]]] == 0, ]
    
    processed_deconv[[column_labels[i]]] <- df
  }
  
  top6signs <- data.frame(
    AML_all = head(processed_deconv[["AML_all"]]$SBS),
    CH_all = head(processed_deconv[["CH_all"]]$SBS),
    AML_DNMT3A = head(processed_deconv[["AML_DNMT3A"]]$SBS),
    CH_DNMT3A = head(processed_deconv[["CH_DNMT3A"]]$SBS),
    AML_TET2 = head(processed_deconv[["AML_TET2"]]$SBS),
    CH_TET2 = head(processed_deconv[["CH_TET2"]]$SBS)
  )
  
  return(top6signs)
}

top6signs <- process_deconv(sign_deconv)

#export the table of top 6 signatures
signtable<- nice_table(
  top6signs, 
  title = c("Table 1", "Top six COMSIC mutational signatures reconstructed in the mutation profiles"),
  note = c("By solving the nonnegative least-squares constraints problem, mutation signatures are deconvolved and the mutation matrix is reconstructed."))
flextable::save_as_docx(signtable, path = "C:\\Users\\fazel\\...\\mutSigns\\signtable.docx")
