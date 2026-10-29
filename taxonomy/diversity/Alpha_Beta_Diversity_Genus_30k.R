
library(phyloseq)
library(vegan)
library(permute)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rstatix)

# ---------- INPUT ----------
biom_file     <- "STP_24.biom"          # <-- point at the new BIOM
metadata_file <- "MasterMetadataFinal.csv"

EXCLUDE_SAMPLES <- c("PR-293-22-Jul-2024", "PR-294-29-Jul-2024")
EXPECTED_N      <- 132                  # after PR-280 resolution + exclusions

GENUS_RANK   <- "Rank6"                 # kraken-biom: Rank6 = genus
RAREFY_DEPTH <- 70000

OUT <- "Results_Genus_3"
dir.create(OUT, showWarnings=FALSE)
dir.create(file.path(OUT,"Rarefied"),        recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT,"Statistics"),      recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT,"Alpha_Diversity"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT,"Beta_Diversity"),  recursive=TRUE, showWarnings=FALSE)

# ============================================================
# IMPORT AND SAMPLE RESOLUTION
# ============================================================

ps <- import_biom(biom_file)
cat("\n=== IMPORT ===\n")
cat("Samples in BIOM at import:", nsamples(ps), "\n")

meta <- read.csv(metadata_file, check.names=FALSE, stringsAsFactors=FALSE)
meta <- meta %>% filter(year %in% c(2024))
rownames(meta) <- meta$Sample_Date_ID
cat("Metadata rows for 2024:", nrow(meta), "\n\n")

# ---------- PR-280: keep the re-sequenced replicate, drop the original ----------
# Done BEFORE metadata attachment. The replicate is renamed to the original's
# ID so it inherits the correct metadata row (its own ID is formatted
# differently and carries a 2026 date, so it has no metadata row of its own).
cat("=== PR-280 RESOLUTION ===\n")
pr280_all  <- grep("^PR-280", sample_names(ps), value=TRUE)
if (length(pr280_all) == 0) cat("No PR-280 samples found.\n") else print(pr280_all)

# '^PR-280-2[-_]' matches PR-280-2-<date> but NOT PR-280-22-<mon>
pr280_rep  <- grep("^PR-280-2[-_]", pr280_all, value=TRUE)
pr280_orig <- setdiff(pr280_all, pr280_rep)

if (length(pr280_rep) == 1 && length(pr280_orig) == 1) {
  cat("  dropping:", pr280_orig, "\n")
  cat("  keeping :", pr280_rep, "-> renamed to", pr280_orig, "\n")
  ps <- prune_samples(sample_names(ps) != pr280_orig, ps)
  nm <- sample_names(ps)
  nm[nm == pr280_rep] <- pr280_orig
  sample_names(ps) <- nm
} else if (length(pr280_rep) == 1 && length(pr280_orig) == 0) {
  cat("  only the replicate is present; renaming to PR-280 form.\n")
  nm <- sample_names(ps)
  nm[nm == pr280_rep] <- sub("^PR-280-2", "PR-280", pr280_rep)
  sample_names(ps) <- nm
} else if (length(pr280_rep) == 0 && length(pr280_orig) == 1) {
  cat("  only one PR-280 sample present; nothing to swap.\n")
} else if (length(pr280_all) > 0) {
  warning("Unexpected PR-280 pattern (", length(pr280_orig), " original, ",
          length(pr280_rep), " replicate). NOT modified - inspect manually.")
}
cat("Samples remaining:", nsamples(ps), "\n\n")

# ---------- named exclusions ----------
cat("=== NAMED EXCLUSIONS ===\n")
found_excl <- intersect(EXCLUDE_SAMPLES, sample_names(ps))
not_found  <- setdiff(EXCLUDE_SAMPLES, sample_names(ps))
if (length(found_excl) > 0) cat("Excluding:", paste(found_excl, collapse=", "), "\n")
if (length(not_found)  > 0)
  warning("Requested for exclusion but NOT found: ",
          paste(not_found, collapse=", "))
ps <- prune_samples(!(sample_names(ps) %in% EXCLUDE_SAMPLES), ps)
cat("Samples remaining:", nsamples(ps), "\n\n")

# ---------- metadata: fail loudly rather than dropping silently ----------
missing_meta <- setdiff(sample_names(ps), rownames(meta))
if (length(missing_meta) > 0) {
  stop("No metadata row for: ", paste(missing_meta, collapse=", "),
       "\n  phyloseq would drop these silently. Fix the metadata or the ",
       "sample name before continuing.")
}
n_before <- nsamples(ps)
sample_data(ps) <- sample_data(meta[sample_names(ps), ])
stopifnot(nsamples(ps) == n_before)
cat("Samples after metadata attachment:", nsamples(ps), "\n")

if (nsamples(ps) != EXPECTED_N)
  warning("Sample count is ", nsamples(ps), ", expected ", EXPECTED_N,
          ". Check the output above.")

# ---------- Bacteria only ----------
ps <- subset_taxa(ps, grepl("Bacteria", Rank1))
ps <- prune_taxa(taxa_sums(ps)>0, ps)
cat("Taxa after Bacteria filter:", ntaxa(ps), "\n")
reads_species_level <- sum(sample_sums(ps))

# ---------- clean bare rank prefixes to NA ----------
# kraken-biom writes unassigned levels as "g__" rather than NA;
# without this, tax_glom lumps every unassigned taxon into one bin.
tt <- as.data.frame(tax_table(ps), stringsAsFactors=FALSE)
tt[] <- lapply(tt, function(x){
  x <- trimws(as.character(x))
  x[x %in% c("","NA","na")] <- NA
  x[grepl("^[a-z]__$", x)]  <- NA
  x
})
tax_table(ps) <- tax_table(as.matrix(tt))

if (!GENUS_RANK %in% colnames(tax_table(ps)))
  stop("Rank column '", GENUS_RANK, "' not found. Columns present: ",
       paste(colnames(tax_table(ps)), collapse=", "))

# ---------- collapse to genus ----------
cat("\nCollapsing to genus (", GENUS_RANK, ") ...\n")
ps.genus <- tax_glom(ps, taxrank=GENUS_RANK, NArm=TRUE)
ps.genus <- prune_taxa(taxa_sums(ps.genus)>0, ps.genus)

reads_genus_level <- sum(sample_sums(ps.genus))
cat("  genera:", ntaxa(ps.genus), "\n")
cat("  reads retained at genus level:",
    round(100*reads_genus_level/reads_species_level, 1),
    "% (reads unassigned at genus were dropped)\n")

# readable taxa names
gn <- as.character(tax_table(ps.genus)[, GENUS_RANK])
gn <- gsub("^g__", "", gn)
taxa_names(ps.genus) <- make.unique(gn)

saveRDS(ps.genus, file.path(OUT,"Rarefied","ps_genus_prerarefy.rds"))

# ---------- depth AFTER collapsing (this is what RAREFY_DEPTH applies to) ----------
gd <- sample_sums(ps.genus)
cat("\n=== GENUS-ASSIGNED READ DEPTH ===\n")
print(summary(as.numeric(gd)))
depth_tab <- data.frame(Sample=names(gd), Genus_Reads=as.numeric(gd))
depth_tab <- depth_tab[order(depth_tab$Genus_Reads), ]
write.csv(depth_tab, file.path(OUT,"Rarefied","Genus_Read_Depths.csv"),
          row.names=FALSE)

dropped <- depth_tab$Sample[depth_tab$Genus_Reads < RAREFY_DEPTH]
cat("Samples below", RAREFY_DEPTH, "genus reads:", length(dropped),
    "of", length(gd), "\n")
if (length(dropped) > 0) {
  lowdf <- depth_tab[depth_tab$Genus_Reads < RAREFY_DEPTH, ]
  lowdf <- cbind(lowdf,
                 data.frame(sample_data(ps.genus))[lowdf$Sample,
                                                   c("stp","Month_Name","season")])
  print(lowdf, row.names=FALSE)
  cat("  -> check these are not concentrated in one STP or season\n")
}
cat("\n")

# ---------- rarefy ----------
set.seed(123)
ps.rarefied <- rarefy_even_depth(ps.genus, sample.size=RAREFY_DEPTH,
                                 rngseed=123, replace=FALSE, verbose=TRUE)
saveRDS(ps.rarefied, file.path(OUT,"Rarefied","Rarefied_genus_phyloseq.rds"))
cat("Retained:", nsamples(ps.rarefied), "samples,",
    ntaxa(ps.rarefied), "genera\n\n")
write.csv(data.frame(Sample=sample_names(ps.rarefied)),
          file.path(OUT,"Rarefied","Samples_Retained.csv"), row.names=FALSE)

# ---------- order factors ----------
sd <- data.frame(sample_data(ps.rarefied))
sd$Month_Name <- factor(sd$Month_Name,
                        levels=c("January","February","March","April","May","June",
                                 "July","August","September","October","November","December"))
sd$season <- factor(sd$season, levels=c("Summer","Monsoon","Winter"))
sample_data(ps.rarefied) <- sample_data(sd)

# ============================================================
# ALPHA
# ============================================================
alpha <- estimate_richness(ps.rarefied,
                           measures=c("Observed","Chao1","Shannon","Simpson"))
alpha <- cbind(alpha, data.frame(sample_data(ps.rarefied)))
write.csv(alpha, file.path(OUT,"Alpha_Diversity","Alpha_Diversity_Genus.csv"),
          row.names=TRUE)

idx  <- c("Observed","Chao1","Shannon","Simpson")
grps <- c("stp","Month_Name","season")

out <- data.frame()
for(i in idx){
  for(g in grps){
    k <- kruskal.test(alpha[[i]] ~ alpha[[g]])
    out <- rbind(out, data.frame(Index=i, Grouping=g,
                                 Statistic=unname(k$statistic),
                                 Pvalue=k$p.value))
  }}
out$FDR <- p.adjust(out$Pvalue, method="BH")
write.csv(out, file.path(OUT,"Statistics","Alpha_Kruskal_Results.csv"),
          row.names=FALSE)
cat("=== ALPHA (Kruskal-Wallis) ===\n"); print(out, row.names=FALSE, digits=4)

# pairwise Wilcoxon -> CSV only (not plotted)
for(g in grps){
  res_g <- do.call(rbind, lapply(idx, function(i){
    df2 <- alpha[!is.na(alpha[[g]]) & !is.na(alpha[[i]]), ]
    df2[[g]] <- droplevels(as.factor(df2[[g]]))
    if(nlevels(df2[[g]]) < 2) return(NULL)
    f  <- as.formula(paste(i, "~", g))
    pw <- tryCatch(as.data.frame(
      rstatix::pairwise_wilcox_test(df2, f, p.adjust.method="BH")),
      error=function(e) NULL)
    if(is.null(pw) || nrow(pw)==0) return(NULL)
    pw$Index <- i; pw$Grouping <- g
    pw
  }))
  if(!is.null(res_g) && nrow(res_g) > 0)
    write.csv(res_g, file.path(OUT,"Statistics",
                               paste0("Pairwise_Wilcoxon_", g, ".csv")), row.names=FALSE)
}

############################################################
## Boxplots
############################################################
plot_box <- function(group_var, filename){
  
  if(group_var == "Month_Name"){
    alpha[[group_var]] <- factor(alpha[[group_var]],
                                 levels=c("January","February","March","April","May","June",
                                          "July","August","September","October","November","December"))
  } else if(group_var == "season"){
    alpha[[group_var]] <- factor(alpha[[group_var]],
                                 levels=c("Summer","Monsoon","Winter"))
  } else {
    alpha[[group_var]] <- factor(alpha[[group_var]],
                                 levels=unique(alpha[[group_var]]))
  }
  
  plots <- list()
  for(i in idx){
    p <- ggplot(alpha, aes_string(x=group_var, y=i, fill=group_var)) +
      geom_boxplot(width=0.65, outlier.shape=NA, colour="black",
                   linewidth=0.6, alpha=0.9) +
      geom_jitter(width=0.15, size=2.2, alpha=0.85, colour="black") +
      scale_fill_viridis_d(option="turbo") +
      labs(x="", y=i) +
      theme_classic(base_size=15) +
      theme(legend.position="none",
            axis.text.x=element_text(angle=45,hjust=1,face="bold",colour="black"),
            axis.text.y=element_text(face="bold",colour="black"),
            axis.title.y=element_text(face="bold"))
    plots[[i]] <- p
  }
  
  final_plot <- (plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]])
  ggsave(file.path(OUT,"Alpha_Diversity",paste0(filename,"f.png")),
         final_plot, width=13, height=10, dpi=600)
  ggsave(file.path(OUT,"Alpha_Diversity",paste0(filename,"f.pdf")),
         final_plot, width=13, height=10)
}

plot_box("stp","Alpha_Box_STP")
plot_box("Month_Name","Alpha_Box_Month")
plot_box("season","Alpha_Box_Season")

############################################################
## Violin plots
############################################################
plot_violin <- function(group_var, filename){
  
  if(group_var == "Month_Name"){
    alpha[[group_var]] <- factor(alpha[[group_var]],
                                 levels=c("January","February","March","April","May","June",
                                          "July","August","September","October","November","December"))
  } else if(group_var == "season"){
    alpha[[group_var]] <- factor(alpha[[group_var]],
                                 levels=c("Summer","Monsoon","Winter"))
  } else {
    alpha[[group_var]] <- factor(alpha[[group_var]],
                                 levels=unique(alpha[[group_var]]))
  }
  
  plots <- list()
  for(i in idx){
    p <- ggplot(alpha, aes_string(x=group_var, y=i, fill=group_var)) +
      geom_violin(trim=FALSE, colour="black", alpha=0.70, linewidth=0.5) +
      geom_boxplot(width=0.12, fill="white", colour="black", outlier.shape=NA) +
      geom_jitter(width=0.15, size=2, alpha=0.8, colour="black") +
      scale_fill_viridis_d(option="turbo") +
      labs(x="", y=i) +
      theme_classic(base_size=15) +
      theme(legend.position="none",
            axis.text.x=element_text(angle=45,hjust=1,face="bold",colour="black"),
            axis.text.y=element_text(face="bold",colour="black"),
            axis.title.y=element_text(face="bold"))
    plots[[i]] <- p
  }
  
  final_plot <- (plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]])
  ggsave(file.path(OUT,"Alpha_Diversity",paste0(filename,"fb.png")),
         final_plot, width=13, height=10, dpi=600)
  ggsave(file.path(OUT,"Alpha_Diversity",paste0(filename,"fb.pdf")),
         final_plot, width=13, height=10)
}

plot_violin("stp","Alpha_Violin_STP")
plot_violin("Month_Name","Alpha_Violin_Month")
plot_violin("season","Alpha_Violin_Season")

# ============================================================
# BETA
# ============================================================
bray <- phyloseq::distance(ps.rarefied, "bray")
write.csv(as.matrix(bray),
          file.path(OUT,"Statistics","BrayCurtis_Distance_Genus.csv"))

ord <- ordinate(ps.rarefied, "PCoA", "bray")

eig <- ord$values$Relative_eig
write.csv(data.frame(Axis=seq_along(eig), Percent_Variance=100*eig),
          file.path(OUT,"Statistics","PCoA_Variance_Explained.csv"),
          row.names=FALSE)
cat("\nPCoA axis 1:", round(100*eig[1],1),
    "%   axis 2:", round(100*eig[2],1), "%\n")

coords <- plot_ordination(ps.rarefied, ord, type="samples", justDF=TRUE)
write.csv(coords, file.path(OUT,"Statistics","PCoA_Coordinates.csv"),
          row.names=FALSE)

# ---- confound check: is axis 1 tracking sequencing depth? ----
pre_depth <- sample_sums(ps.genus)[rownames(coords)]
ct <- suppressWarnings(cor.test(coords$Axis.1, as.numeric(pre_depth),
                                method="spearman"))
cat("\n=== CONFOUND CHECK ===\n")
cat("PCoA axis 1 vs pre-rarefaction genus depth: rho =",
    round(unname(ct$estimate), 3), " p =", signif(ct$p.value, 3), "\n")
cat("  (a strong correlation would mean axis 1 is technical, not biological)\n")
write.csv(data.frame(rho=unname(ct$estimate), p=ct$p.value),
          file.path(OUT,"Statistics","Axis1_vs_Depth.csv"), row.names=FALSE)

# ---- PERMANOVA + PERMDISP ----
perm <- data.frame(); perm_labels <- list()

for(g in grps){
  f  <- as.formula(paste("bray ~", g))
  ad <- adonis2(f, data=sd, permutations=999)
  bd <- betadisper(bray, sd[[g]])
  pt <- permutest(bd)
  
  perm <- rbind(perm, data.frame(
    Group=g,
    PERMANOVA_F  = ad$F[1],
    PERMANOVA_R2 = ad$R2[1],
    PERMANOVA_P  = ad$`Pr(>F)`[1],
    Betadisper_F = pt$tab[1,"F"],
    Betadisper_P = pt$tab[1,"Pr(>F)"]))
  
  pv <- ad$`Pr(>F)`[1]
  perm_labels[[g]] <- sprintf(
    "PERMANOVA: R\u00B2 = %.3f, p = %s",
    ad$R2[1],
    ifelse(pv < 0.001, "< 0.001", formatC(pv, format="f", digits=3)))
}

write.csv(perm, file.path(OUT,"Statistics","PERMANOVA_Betadisper.csv"),
          row.names=FALSE)
cat("\n=== PERMANOVA (genus level) ===\n"); print(perm, row.names=FALSE, digits=4)

# ---- STP-blocked PERMANOVA for the temporal terms ----
# Samples are repeated measures within STP; restricting permutations to
# within-STP blocks respects that structure. Not applicable to stp itself.
blk <- data.frame()
for(g in c("Month_Name","season")){
  f  <- as.formula(paste("bray ~", g))
  ctrl <- permute::how(nperm = 999, blocks = sd$stp)
  adb <- adonis2(f, data=sd, permutations=ctrl)
  blk <- rbind(blk, data.frame(Group=g, F=adb$F[1], R2=adb$R2[1],
                               P_blocked=adb$`Pr(>F)`[1]))
}
write.csv(blk, file.path(OUT,"Statistics","PERMANOVA_Blocked_by_STP.csv"),
          row.names=FALSE)
cat("\n=== PERMANOVA, permutations blocked within STP ===\n")
print(blk, row.names=FALSE, digits=4)

# ---- pairwise PERMANOVA ----
pairwise_adonis <- function(dist_obj, groups, perm_n=999, adj="BH"){
  groups <- droplevels(as.factor(groups)); dm <- as.matrix(dist_obj)
  combos <- combn(levels(groups), 2, simplify=FALSE)
  res <- do.call(rbind, lapply(combos, function(cc){
    keep <- groups %in% cc
    dsub <- as.dist(dm[keep, keep])
    dfs  <- data.frame(g = droplevels(groups[keep]))
    a <- adonis2(dsub ~ g, data=dfs, permutations=perm_n)
    data.frame(Group1=cc[1], Group2=cc[2],
               F=a$F[1], R2=a$R2[1], p=a$`Pr(>F)`[1])
  }))
  res$p.adj <- p.adjust(res$p, method=adj)
  res[order(res$p.adj), ]
}

for(g in grps){
  pa <- pairwise_adonis(bray, sd[[g]])
  write.csv(pa, file.path(OUT,"Statistics",
                          paste0("Pairwise_PERMANOVA_", g, ".csv")), row.names=FALSE)
}

# ---- PCoA plots with PERMANOVA annotation ----
mk_pcoa <- function(colour_var, subtitle, shape_var=NULL, ellipse=TRUE){
  p <- if(is.null(shape_var))
    plot_ordination(ps.rarefied, ord, color=colour_var)
  else
    plot_ordination(ps.rarefied, ord, color=colour_var, shape=shape_var)
  p <- p + geom_point(size=3)
  if(ellipse) p <- p + stat_ellipse()
  p + labs(subtitle=subtitle) +
    theme_classic(base_size=14) +
    theme(plot.subtitle=element_text(size=11, face="italic"),
          axis.text=element_text(colour="black"),
          axis.title=element_text(face="bold"))
}

p1 <- mk_pcoa("stp",        perm_labels[["stp"]])
p2 <- mk_pcoa("Month_Name", perm_labels[["Month_Name"]])
p3 <- mk_pcoa("season",     perm_labels[["season"]])
p4 <- mk_pcoa("stp",        perm_labels[["stp"]], shape_var="season",
              ellipse=FALSE)

bp <- file.path(OUT,"Beta_Diversity")
ggsave(file.path(bp,"PCOA_STP.png"),        p1, width=9, height=6.5, dpi=600)
ggsave(file.path(bp,"PCOA_STP.pdf"),        p1, width=9, height=6.5)
ggsave(file.path(bp,"PCOA_Month.png"),      p2, width=9, height=6.5, dpi=600)
ggsave(file.path(bp,"PCOA_Month.pdf"),      p2, width=9, height=6.5)
ggsave(file.path(bp,"PCOA_Season.png"),     p3, width=9, height=6.5, dpi=600)
ggsave(file.path(bp,"PCOA_Season.pdf"),     p3, width=9, height=6.5)
ggsave(file.path(bp,"PCOA_STP_Season.png"), p4, width=9, height=6.5, dpi=600)
ggsave(file.path(bp,"PCOA_STP_Season.pdf"), p4, width=9, height=6.5)

writeLines(capture.output(sessionInfo()),
           file.path(OUT,"Statistics","sessionInfo.txt"))
cat("\nDone. Outputs in", OUT, "\n")