
library(vegan)
library(viridis)
library(patchwork)

arg_file <- "AMR-Matrix.tsv"

metadata_file <- "MasterMetadataFinal.csv"

#==============================================================
# Output folders
#==============================================================

dir.create(
  "Results_ARG_2",
  showWarnings = FALSE
)

dir.create(
  "Results_ARG_2/Objects",
  recursive = TRUE,
  showWarnings = FALSE
)

#==============================================================
# Read ARG matrix
#==============================================================

arg <- read.delim(
  arg_file,
  sep="\t",
  header=TRUE,
  check.names=FALSE,
  stringsAsFactors=FALSE
)

###############################################################
## Remove ARGs with Unclassified antibiotic class
###############################################################

cat("Initial ARGs :", nrow(arg), "\n")

arg <- arg %>%
  filter(!is.na(.[[2]])) %>%                    # remove missing classes
  filter(trimws(.[[2]]) != "") %>%              # remove blank classes
  filter(tolower(trimws(.[[2]])) != "unclassified")

cat("After removing 'Unclassified' ARGs :", nrow(arg), "\n\n")

cat("\n")

cat("Rows :",nrow(arg),"\n")

cat("Columns :",ncol(arg),"\n\n")

#==============================================================
# Save ARG annotation separately
#==============================================================

arg_annotation <- arg[,1:2]

colnames(arg_annotation) <- c(
  "ARG",
  "ARG_Class"
)

#==============================================================
# Create abundance matrix
#==============================================================

arg_matrix <- arg[,-c(1,2)]

rownames(arg_matrix) <- make.unique(arg[[1]])

#==============================================================
# Convert to numeric
#==============================================================

arg_matrix[] <- lapply(
  arg_matrix,
  as.numeric
)

arg_matrix <- as.data.frame(arg_matrix)

#==============================================================
# Read metadata
#==============================================================

meta <- read.csv(
  metadata_file,
  check.names=FALSE,
  stringsAsFactors=FALSE
)

#==============================================================
# Filter required year
#==============================================================

## For 2024
meta <- meta %>% filter(year %in% c(2024))

## For 2022-2023

# meta <- meta %>%
# filter(year %in% c(2022,2023))

#==============================================================
# Check metadata columns
#==============================================================

required_columns <- c(
  
  "Sample_Date_ID",
  
  "stp",
  
  "Month_Name",
  
  "season"
  
)

missing_columns <- setdiff(
  required_columns,
  colnames(meta)
)

if(length(missing_columns)>0){
  
  stop(
    
    paste(
      
      "Missing columns:",
      
      paste(
        missing_columns,
        collapse=", "
      )
      
    )
    
  )
  
}

#==============================================================
# Match samples
#==============================================================

common_samples <- intersect(
  
  colnames(arg_matrix),
  
  meta$Sample_Date_ID
  
)

cat("Matched Samples :",length(common_samples),"\n")

if(length(common_samples)==0){
  
  stop("No matching sample names found.")
  
}

arg_matrix <- arg_matrix[,
                         
                         common_samples,
                         
                         drop=FALSE]

meta <- meta %>%
  
  filter(
    
    Sample_Date_ID %in%
      
      common_samples
    
  )

meta <- meta[
  match(
    common_samples,
    meta$Sample_Date_ID
  ),
]

rownames(meta) <- meta$Sample_Date_ID

#==============================================================
# Factor ordering
#==============================================================

meta$Month_Name <- factor(
  
  meta$Month_Name,
  
  levels=c(
    
    "January",
    "February",
    "March",
    "April",
    "May",
    "June",
    "July",
    "August",
    "September",
    "October",
    "November",
    "December"
    
  )
  
)

meta$season <- factor(
  
  meta$season,
  
  levels=c(
    
    "Summer",
    
    "Monsoon",
    
    "Winter"
    
  )
  
)

#==============================================================
# Remove ARGs with zero abundance
#==============================================================

arg_matrix <- arg_matrix[
  
  rowSums(arg_matrix)>0,
  
]

cat("\n")

cat("Remaining ARGs :",nrow(arg_matrix),"\n")

cat("Remaining Samples :",ncol(arg_matrix),"\n\n")

#==============================================================
# Save cleaned objects
#==============================================================

saveRDS(
  
  arg_matrix,
  
  "Results_ARG_2/Objects/ARG_matrix.rds"
  
)

saveRDS(
  
  meta,
  
  "Results_ARG_2/Objects/Metadata.rds"
  
)

saveRDS(
  
  arg_annotation,
  
  "Results_ARG_2/Objects/ARG_annotation.rds"
  
)

write.csv(
  
  arg_annotation,
  
  "Results_ARG_2/Objects/ARG_annotation.csv",
  
  row.names=FALSE
  
)

cat("Objects saved successfully.\n")

###############################################################
## Script 2
## Alpha Diversity Calculation
###############################################################

#==============================================================
# Output folder
#==============================================================

dir.create(
  "Results_ARG_2/Alpha_Diversity",
  recursive = TRUE,
  showWarnings = FALSE
)

#==============================================================
# Load cleaned objects
#==============================================================

arg_matrix <- readRDS(
  "Results_ARG_2/Objects/ARG_matrix.rds"
)

meta <- readRDS(
  "Results_ARG_2/Objects/Metadata.rds"
)

#==============================================================
# Convert to matrix
#==============================================================

arg_matrix <- as.matrix(arg_matrix)

storage.mode(arg_matrix) <- "numeric"

#==============================================================
# Transpose matrix
# Samples become rows
#==============================================================

arg_t <- t(arg_matrix)
Observed <- rowSums(arg_t > 0)

Shannon <- vegan::diversity(
  arg_t,
  index = "shannon",
  MARGIN = 1
)

Simpson <- vegan::diversity(
  arg_t,
  index = "simpson",
  MARGIN = 1
)

#==============================================================
# Combine results
#==============================================================
alpha <- data.frame(
  
  Sample = rownames(arg_t),
  
  Observed = Observed,
  
  Shannon = Shannon,
  
  Simpson = Simpson,
  
  stringsAsFactors = FALSE
  
)

#==============================================================
# Merge with metadata
#==============================================================

alpha <- left_join(
  
  alpha,
  
  meta,
  
  by = c(
    "Sample" = "Sample_Date_ID"
  )
  
)

#==============================================================
# Arrange samples
#==============================================================

alpha <- alpha %>%
  
  arrange(
    
    Month_Name,
    
    stp
    
  )

#==============================================================
# Save results
#==============================================================

write.csv(
  
  alpha,
  
  "Results_ARG_2/Alpha_Diversity/Alpha_Diversity.csv",
  
  row.names = FALSE
  
)

#==============================================================
# Summary
#==============================================================

cat("\n")

cat("=============================\n")

cat("Alpha Diversity Summary\n")

cat("=============================\n")

cat("Samples :", nrow(alpha), "\n")

cat("Observed richness range :",
    
    round(min(alpha$Observed),2),
    
    "-",
    
    round(max(alpha$Observed),2),
    
    "\n")

cat("Shannon range :",
    
    round(min(alpha$Shannon),2),
    
    "-",
    
    round(max(alpha$Shannon),2),
    
    "\n")

cat("Simpson range :",
    
    round(min(alpha$Simpson),3),
    
    "-",
    
    round(max(alpha$Simpson),3),
    
    "\n")

cat("=============================\n\n")

cat("Alpha diversity file saved.\n")

###############################################################
## Script 3
## Kruskal-Wallis Tests for Alpha Diversity
###############################################################

#==============================================================
# Output folder
#==============================================================

dir.create(
  "Results_ARG_2/Statistics",
  recursive = TRUE,
  showWarnings = FALSE
)

#==============================================================
# Read Alpha Diversity
#==============================================================

alpha <- read.csv(
  "Results_ARG_2/Alpha_Diversity/Alpha_Diversity.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

#==============================================================
# Preserve factor order
#==============================================================

alpha$Month_Name <- factor(
  alpha$Month_Name,
  levels = c(
    "January","February","March","April",
    "May","June","July","August",
    "September","October","November","December"
  )
)

alpha$season <- factor(
  alpha$season,
  levels = c(
    "Summer",
    "Monsoon",
    "Winter"
  )
)

alpha$stp <- factor(alpha$stp)

#==============================================================
# Diversity indices
#==============================================================

indices <- c(
  "Observed",
  "Shannon",
  "Simpson"
)

#==============================================================
# Metadata groupings
#==============================================================

groupings <- c(
  "stp",
  "Month_Name",
  "season"
)

#==============================================================
# Run Kruskal-Wallis tests
#==============================================================

kw_results <- data.frame()

for(ind in indices){
  
  for(grp in groupings){
    
    formula <- as.formula(
      paste(ind, "~", grp)
    )
    
    kw <- kruskal.test(
      formula,
      data = alpha
    )
    
    kw_results <- rbind(
      kw_results,
      data.frame(
        Diversity_Index = ind,
        Grouping = grp,
        Chi_Square = unname(kw$statistic),
        DF = unname(kw$parameter),
        P_value = kw$p.value,
        stringsAsFactors = FALSE
      )
    )
    
  }
  
}

#==============================================================
# Adjust p-values (Benjamini-Hochberg FDR)
#==============================================================

kw_results$FDR <- p.adjust(
  kw_results$P_value,
  method = "BH"
)

#==============================================================
# Significance labels
#==============================================================

kw_results$Significance <- cut(
  kw_results$P_value,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "ns")
)

#==============================================================
# Save results
#==============================================================

write.csv(
  kw_results,
  "Results_ARG_2/Statistics/Alpha_Kruskal_Results.csv",
  row.names = FALSE
)

#==============================================================
# Print summary
#==============================================================

cat("\n")
cat("=====================================\n")
cat("Kruskal-Wallis Tests Completed\n")
cat("=====================================\n\n")

print(kw_results)

cat("\nResults saved to:\n")
cat("Results_ARG_2/Statistics/Alpha_Kruskal_Results.csv\n")



###############################################################
## Script 4
## Alpha Diversity Plots
###############################################################

#==============================================================
# Libraries
#==============================================================

library(viridis)
#==============================================================
# Output folders
#==============================================================

dir.create(
  "Results_ARG_2/Figures/Alpha_Diversity/Boxplots",
  recursive=TRUE,
  showWarnings=FALSE
)

dir.create(
  "Results_ARG_2/Figures/Alpha_Diversity/Violinplots",
  recursive=TRUE,
  showWarnings=FALSE
)

#==============================================================
# Read alpha diversity
#==============================================================

alpha <- read.csv(
  "Results_ARG_2/Alpha_Diversity/Alpha_Diversity.csv",
  check.names=FALSE
)

#==============================================================
# Factor ordering
#==============================================================

alpha$Month_Name <- factor(
  alpha$Month_Name,
  levels=c(
    "January","February","March","April",
    "May","June","July","August",
    "September","October","November","December"
  )
)

alpha$season <- factor(
  alpha$season,
  levels=c(
    "Summer",
    "Monsoon",
    "Winter"
  )
)

alpha$stp <- factor(alpha$stp)

#==============================================================
# Theme
# Legend suppressed here because the x-axis already
# labels every group in the alpha diversity plots
#==============================================================

theme_pub <- theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

#==============================================================
# Function
#==============================================================

plot_alpha <- function(group_var,
                       group_name){
  
  n_groups <- length(unique(alpha[[group_var]]))
  
  pal <- viridis(
    n_groups,
    option="D"
  )
  
  p1 <- ggplot(
    alpha,
    aes_string(
      x=group_var,
      y="Shannon",
      fill=group_var
    )
  )+
    geom_boxplot(
      width=0.7,
      outlier.shape=16
    )+
    scale_fill_manual(values=pal)+
    labs(
      title="Shannon",
      x=NULL,
      y="Shannon Index"
    )+
    theme_pub
  
  p2 <- ggplot(
    alpha,
    aes_string(
      x=group_var,
      y="Simpson",
      fill=group_var
    )
  )+
    geom_boxplot(
      width=0.7,
      outlier.shape=16
    )+
    scale_fill_manual(values=pal)+
    labs(
      title="Simpson",
      x=NULL,
      y="Simpson Index"
    )+
    theme_pub
  
  fig <- p1 + p2
  
  ggsave(
    paste0(
      "Results_ARG_2/Figures/Alpha_Diversity/Boxplots/Alpha_Boxplot_",
      group_name,
      ".png"
    ),
    fig,
    width=12,
    height=5,
    dpi=600
  )
  
  ggsave(
    paste0(
      "Results_ARG_2/Figures/Alpha_Diversity/Boxplots/Alpha_Boxplot_",
      group_name,
      ".pdf"
    ),
    fig,
    width=12,
    height=5
  )
  
  ############################################################
  
  v1 <- ggplot(
    alpha,
    aes_string(
      x=group_var,
      y="Shannon",
      fill=group_var
    )
  )+
    geom_violin(trim=FALSE)+
    geom_boxplot(
      width=0.12,
      fill="white"
    )+
    scale_fill_manual(values=pal)+
    labs(
      title="Shannon",
      x=NULL,
      y="Shannon Index"
    )+
    theme_pub
  
  v2 <- ggplot(
    alpha,
    aes_string(
      x=group_var,
      y="Simpson",
      fill=group_var
    )
  )+
    geom_violin(trim=FALSE)+
    geom_boxplot(
      width=0.12,
      fill="white"
    )+
    scale_fill_manual(values=pal)+
    labs(
      title="Simpson",
      x=NULL,
      y="Simpson Index"
    )+
    theme_pub
  
  fig2 <- v1 + v2
  
  ggsave(
    paste0(
      "Results_ARG_2/Figures/Alpha_Diversity/Violinplots/Alpha_Violin_",
      group_name,
      ".png"
    ),
    fig2,
    width=12,
    height=5,
    dpi=600
  )
  
  ggsave(
    paste0(
      "Results_ARG_2/Figures/Alpha_Diversity/Violinplots/Alpha_Violin_",
      group_name,
      ".pdf"
    ),
    fig2,
    width=12,
    height=5
  )
  
}

#==============================================================
# Generate figures
#==============================================================

plot_alpha(
  "stp",
  "STP"
)

plot_alpha(
  "Month_Name",
  "Month"
)

plot_alpha(
  "season",
  "Season"
)

cat("\nAlpha diversity plots completed.\n")

###############################################################
## Script 5
## Bray-Curtis Distance
###############################################################

#==============================================================
# Output folder
#==============================================================

dir.create(
  "Results_ARG_2/Beta_Diversity",
  recursive = TRUE,
  showWarnings = FALSE
)

#==============================================================
# Read cleaned objects
#==============================================================

arg_matrix <- readRDS(
  "Results_ARG_2/Objects/ARG_matrix.rds"
)

meta <- readRDS(
  "Results_ARG_2/Objects/Metadata.rds"
)

#==============================================================
# Convert to numeric matrix
#==============================================================

arg_matrix <- as.matrix(arg_matrix)

storage.mode(arg_matrix) <- "numeric"

#==============================================================
# Samples should be rows
#==============================================================

arg_comm <- t(arg_matrix)

cat("\n")
cat("Community Matrix Dimensions\n")
cat("---------------------------\n")
cat("Samples :", nrow(arg_comm), "\n")
cat("ARGs    :", ncol(arg_comm), "\n\n")

#==============================================================
# Check sample order
#==============================================================

if(!all(rownames(arg_comm) == rownames(meta))){
  
  stop("Sample order mismatch between ARG matrix and metadata.")
  
}

cat("Sample order verified.\n\n")

###############################################################
## Bray-Curtis distance
## Abundances are already depth-normalised (ARG hits per Gbp)
###############################################################

bray_dist <- vegan::vegdist(
  arg_comm,
  method = "bray"
)

#==============================================================
# Read PERMANOVA results (from Script 6)
#==============================================================

permanova_results <- read.csv(
  "Results_ARG_2/Statistics/PERMANOVA_Betadisper.csv",
  stringsAsFactors = FALSE
)
###############################################################
## Save distance object
###############################################################

saveRDS(
  bray_dist,
  "Results_ARG_2/Beta_Diversity/BrayCurtis_Distance.rds"
)

###############################################################
## Save distance matrix
###############################################################

write.csv(
  as.matrix(bray_dist),
  "Results_ARG_2/Beta_Diversity/BrayCurtis_Distance.csv"
)

cat(
  "Bray-Curtis distance calculated.\n"
)
###############################################################
## Script 6
## PERMANOVA + Betadisper
###############################################################

#==============================================================
# Output folder
#==============================================================

dir.create(
  "Results_ARG_2/Statistics",
  recursive = TRUE,
  showWarnings = FALSE
)

#==============================================================
# Read files
#==============================================================

meta <- readRDS(
  "Results_ARG_2/Objects/Metadata.rds"
)

bray <- readRDS(
  "Results_ARG_2/Beta_Diversity/BrayCurtis_Distance.rds"
)

#==============================================================
# Factor ordering
#==============================================================

meta$Month_Name <- factor(
  meta$Month_Name,
  levels = c(
    "January","February","March","April",
    "May","June","July","August",
    "September","October","November","December"
  )
)

meta$season <- factor(
  meta$season,
  levels = c(
    "Summer",
    "Monsoon",
    "Winter"
  )
)

meta$stp <- factor(meta$stp)

#==============================================================
# Variables
#==============================================================

groups <- c(
  "stp",
  "Month_Name",
  "season"
)

results <- data.frame()

#==============================================================
# Loop through metadata variables
#==============================================================

for(g in groups){
  
  cat("\nRunning:", g, "\n")
  
  ############################################################
  ## PERMANOVA
  ############################################################
  
  form <- as.formula(
    paste("bray ~", g)
  )
  
  permanova <- adonis2(
    form,
    data = meta,
    permutations = 9999
  )
  
  ############################################################
  ## Betadisper
  ############################################################
  
  bd <- betadisper(
    bray,
    meta[[g]]
  )
  
  bd_test <- permutest(
    bd,
    permutations = 9999
  )
  
  ############################################################
  ## Store results
  ############################################################
  
  results <- rbind(
    results,
    data.frame(
      Group = g,
      
      PERMANOVA_F =
        permanova$F[1],
      
      PERMANOVA_R2 =
        permanova$R2[1],
      
      PERMANOVA_P =
        permanova$`Pr(>F)`[1],
      
      Betadisper_F =
        bd_test$tab[1,"F"],
      
      Betadisper_P =
        bd_test$tab[1,"Pr(>F)"],
      
      stringsAsFactors = FALSE
    )
  )
  
}

#==============================================================
# FDR correction
#==============================================================

results$PERMANOVA_FDR <- p.adjust(
  results$PERMANOVA_P,
  method = "BH"
)

results$Betadisper_FDR <- p.adjust(
  results$Betadisper_P,
  method = "BH"
)

#==============================================================
# Significance
#==============================================================

results$PERMANOVA_Significance <- cut(
  results$PERMANOVA_P,
  breaks = c(-Inf,0.001,0.01,0.05,Inf),
  labels = c("***","**","*","ns")
)

results$Betadisper_Significance <- cut(
  results$Betadisper_P,
  breaks = c(-Inf,0.001,0.01,0.05,Inf),
  labels = c("***","**","*","ns")
)

#==============================================================
# Save
#==============================================================

write.csv(
  results,
  "Results_ARG_2/Statistics/PERMANOVA_Betadisper.csv",
  row.names = FALSE
)

#==============================================================
# Print summary
#==============================================================

cat("\n")
cat("========================================\n")
cat("PERMANOVA + Betadisper Completed\n")
cat("========================================\n\n")

print(results)

cat("\nResults saved to:\n")
cat("Results_ARG_2/Statistics/PERMANOVA_Betadisper.csv\n")

###############################################################
## Script 7
## PCoA + Beta Diversity Plots
###############################################################

#==============================================================
# Output folder
#==============================================================

dir.create(
  "Results_ARG_2/Figures/Beta_Diversity",
  recursive = TRUE,
  showWarnings = FALSE
)

#==============================================================
# Read files
#==============================================================

meta <- readRDS(
  "Results_ARG_2/Objects/Metadata.rds"
)

bray <- readRDS(
  "Results_ARG_2/Beta_Diversity/BrayCurtis_Distance.rds"
)
#==============================================================
# Factor ordering
#==============================================================

meta$Month_Name <- factor(
  meta$Month_Name,
  levels = c(
    "January","February","March","April",
    "May","June","July","August",
    "September","October","November","December"
  )
)

meta$season <- factor(
  meta$season,
  levels = c(
    "Summer",
    "Monsoon",
    "Winter"
  )
)

meta$stp <- factor(meta$stp)

# setwd("D:/new_diversity_rf_amr/amr/24/")

###############################################################
## PCoA
###############################################################

n_samples <- attr(bray, "Size")

pcoa <- cmdscale(
  bray,
  eig = TRUE,
  add = TRUE,
  k = min(10, n_samples - 1)
)

###############################################################
## Calculate percentage variance explained
###############################################################

positive_eigenvalues <- pcoa$eig[pcoa$eig > 0]

variance <- round(
  100 * positive_eigenvalues /
    sum(positive_eigenvalues),
  2
)

###############################################################
## Create PCoA coordinate dataframe
###############################################################

scores <- as.data.frame(pcoa$points)

# Rename PCoA axes
colnames(scores) <- paste0(
  "PCoA",
  seq_len(ncol(scores))
)

# Add Sample ID
scores$SampleID <- rownames(scores)

###############################################################
## Add metadata to PCoA coordinates
###############################################################

metadata_for_pcoa <- meta[
  rownames(scores),
  ,
  drop = FALSE
]

# Add selected metadata columns if present
metadata_columns <- intersect(
  c(
    "stp",
    "Month_Name",
    "season",
    "Month",
    "year",
    "Sample_Date_ID"
  ),
  colnames(metadata_for_pcoa)
)

scores <- cbind(
  scores,
  metadata_for_pcoa[, metadata_columns, drop = FALSE]
)

###############################################################
## Reorder columns
###############################################################

coordinate_columns <- grep(
  "^PCoA",
  colnames(scores),
  value = TRUE
)

other_columns <- setdiff(
  colnames(scores),
  c(coordinate_columns, "SampleID")
)

scores <- scores[
  ,
  c(
    "SampleID",
    coordinate_columns,
    other_columns
  ),
  drop = FALSE
]

###############################################################
## Save PCoA coordinates as CSV
###############################################################

write.csv(
  scores,
  "Results_ARG_2/Beta_Diversity/PCoA_Coordinates.csv",
  row.names = FALSE
)

cat("\n")
cat("PCoA coordinates saved to:\n")
cat(
  "Results_ARG_2/Beta_Diversity/PCoA_Coordinates.csv\n"
)

###############################################################
## Save PCoA object
###############################################################

saveRDS(
  pcoa,
  "Results_ARG_2/Beta_Diversity/PCoA_BrayCurtis.rds"
)

#==============================================================
# Theme for PCoA
# Legend is kept here because group identity is shown
# only by colour and shape, not by the axes
#==============================================================

theme_pcoa <- theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.key = element_blank(),
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 11
    )
  )
#==============================================================
# Pretty legend titles
#==============================================================

pretty_label <- function(v){
  
  if(is.null(v)) return(NULL)
  
  switch(
    v,
    "stp"        = "STP",
    "Month_Name" = "Month",
    "season"     = "Season",
    v
  )
  
}

#==============================================================
# PERMANOVA annotation
#==============================================================

permanova_label <- function(vars){
  
  out <- sapply(vars, function(v){
    
    row <- permanova_results[
      permanova_results$Group == v,
    ]
    
    if(nrow(row) == 0) return(NA_character_)
    
    p_txt <- if(row$PERMANOVA_P[1] < 0.001){
      "p < 0.001"
    } else {
      paste0(
        "p = ",
        format(round(row$PERMANOVA_P[1], 3), nsmall = 3)
      )
    }
    
    paste0(
      "PERMANOVA (", pretty_label(v), "): R\u00b2 = ",
      format(round(row$PERMANOVA_R2[1], 3), nsmall = 3),
      ", ", p_txt
    )
    
  })
  
  paste(stats::na.omit(out), collapse = "\n")
  
}

#==============================================================
# Generic plotting function
#==============================================================
plot_pcoa <- function(color_var,
                      shape_var = NULL,
                      title,
                      filename,
                      stat_vars = NULL){
  
  if(is.null(stat_vars)){
    stat_vars <- c(color_var, shape_var)
  }
  
  subtitle_text <- permanova_label(stat_vars)
  
  p <- ggplot(
    scores,
    aes_string(
      x = "PCoA1",
      y = "PCoA2",
      colour = color_var,
      shape = shape_var
    )
  ) +
    geom_point(size = 3) +
    stat_ellipse(
      aes_string(group = color_var),
      linewidth = 0.8,
      alpha = 0.3,
      show.legend = FALSE
    ) +
    scale_colour_viridis_d(option = "D") +
    labs(
      title = title,
      subtitle = subtitle_text,
      x = paste0(
        "PCoA1 (",
        variance[1],
        "%)"
      ),
      y = paste0(
        "PCoA2 (",
        variance[2],
        "%)"
      ),
      colour = pretty_label(color_var),
      shape  = pretty_label(shape_var)
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(size = 3.5)
      )
    ) +
    theme_pcoa
  
  ggsave(
    paste0(
      "Results_ARG_2/Figures/Beta_Diversity/",
      filename,
      ".png"
    ),
    p,
    width = 9,
    height = 6,
    dpi = 600
  )
  
  ggsave(
    paste0(
      "Results_ARG_2/Figures/Beta_Diversity/",
      filename,
      ".pdf"
    ),
    p,
    width = 9,
    height = 6
  )
  
}
#==============================================================
# STP
#==============================================================

plot_pcoa(
  color_var = "stp",
  title = "PCoA - STP",
  filename = "PCoA_STP"
)

#==============================================================
# Month
#==============================================================

plot_pcoa(
  color_var = "Month_Name",
  title = "PCoA - Month",
  filename = "PCoA_Month"
)

#==============================================================
# Season
#==============================================================

plot_pcoa(
  color_var = "season",
  title = "PCoA - Season",
  filename = "PCoA_Season"
)

#==============================================================
# STP + Season
#==============================================================

plot_pcoa(
  color_var = "stp",
  shape_var = "season",
  title = "PCoA - STP + Season",
  filename = "PCoA_STP_Season"
)

cat("\n")
cat("=====================================\n")
cat("PCoA plots generated successfully.\n")
cat("=====================================\n")
