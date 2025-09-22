############################################################
# RNA-seq Analysis Pipeline
# Step1: Annotazione GTF
# Step2: PCA
# Step3: DE analisi (Model vs ComBat_seq)
# Step4: Heatmap and Clustering
# Step5: Genomic proximity
############################################################

# --- Librerie necessarie ---
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(DESeq2)
library(sva)
library(VennDiagram)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(cluster)
library(mclust)
library(grid)
library(reshape2)
############################################################
# STEP1 - GTF Annotation
############################################################
project_path <- "/data"
gtf_path <- file.path(project_path, "Homo_sapiens.GRCh38.114.gtf.gz")
gtf <- import(gtf_path)

genes <- gtf[gtf$type == "gene"]
mcols(genes)$gene_label <- ifelse(
  mcols(genes)$gene_biotype == "protein_coding",
  paste0(mcols(genes)$gene_id, ":", mcols(genes)$gene_name),
  mcols(genes)$gene_id
)
genes <- genes[!is.na(mcols(genes)$gene_id)]

gene_gr <- GRanges(
  seqnames = seqnames(genes),
  ranges = ranges(genes),
  strand = strand(genes),
  gene_id   = mcols(genes)$gene_id,
  gene_name = mcols(genes)$gene_name,
  gene_label = mcols(genes)$gene_label,
  biotype = mcols(genes)$gene_biotype
)
print(gene_gr[1:5])

############################################################
# STEP2 - PCA
############################################################

project_path <- "/data"
save_dir <- "/results"

file_gse1 <- file.path(project_path, "GSE244486_raw_counts.tsv.gz")
file_gse2 <- file.path(project_path, "GSE244485_raw_counts.tsv.gz")

raw_counts_gse1 <- read.delim(file_gse1, header = TRUE, sep = "\t", row.names = 1)
raw_counts_gse2 <- read.delim(file_gse2, header = TRUE, sep = "\t", row.names = 1)

colnames(raw_counts_gse1) <- paste0("GSE244486_", colnames(raw_counts_gse1))
colnames(raw_counts_gse2) <- paste0("GSE244485_", colnames(raw_counts_gse2))

common_genes <- intersect(rownames(raw_counts_gse1), rownames(raw_counts_gse2))
all_raw_counts <- cbind(
  raw_counts_gse1[common_genes, ],
  raw_counts_gse2[common_genes, ]
)

all_raw_counts <- all_raw_counts[rowSums(all_raw_counts) > 0, ]
all_raw_counts <- round(as.matrix(all_raw_counts))

sample_info <- data.frame(
  sample_id = colnames(all_raw_counts),
  condition = rep(c("mock", "bystander", "infected"), each = 6, times = 2),
  batch = rep(c("GSE244486", "GSE244485"), each = 18)
)
rownames(sample_info) <- sample_info$sample_id
sample_info$condition <- relevel(factor(sample_info$condition), ref = "mock")
sample_info$batch <- factor(sample_info$batch)

#PCA
run_pca <- function(count_matrix, sample_info, title_suffix, save_dir = "/results") {
  library(ggplot2)
  
  # Crea cartella se non esiste
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  # Trasformazione log2
  log_data <- log2(count_matrix + 1)
  
  # Rimuovi geni con varianza = 0
  log_data <- log_data[apply(log_data, 1, var) > 0, ]
  
  # Rimuovi campioni con varianza = 0
  sample_var <- apply(log_data, 2, var)
  log_data <- log_data[, sample_var > 0]
  
  # Aggiorna sample_info per i campioni rimasti
  sample_info <- sample_info[colnames(log_data), ]
  
  # PCA
  pca <- prcomp(t(log_data), scale. = TRUE)
  variance_explained <- (pca$sdev^2) / sum(pca$sdev^2)
  
  # Dataframe per plotting
  pca_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    batch = sample_info$batch,
    condition = sample_info$condition
  )
  
  # Plot
  p <- ggplot(pca_df, aes(PC1, PC2, color = condition, shape = batch)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(
      title = paste("PCA - Condition + Batch", title_suffix),
      x = paste0("PC1 (", round(variance_explained[1] * 100), "%)"),
      y = paste0("PC2 (", round(variance_explained[2] * 100), "%)")
    ) +
    theme_minimal()
  
  # Salvataggio PNG (senza spazi nel nome file)
  outfile <- file.path(save_dir, paste0("PCA_", gsub(" ", "_", title_suffix), ".png"))
  message(">>> Salvo PCA in: ", outfile)
  
  ggsave(
    filename = outfile,
    plot = p,
    width = 6,
    height = 5
  )
}


run_pca(all_raw_counts, sample_info, "Raw data")


############################################################
# STEP3 - Batch effect correction + DE analysis
############################################################

# Correzione batch con ComBat_seq
genes_to_keep <- rowSums(all_raw_counts[, 1:18]) > 0 & rowSums(all_raw_counts[, 19:36]) > 0
filtered_counts <- all_raw_counts[genes_to_keep, ]
filtered_sample_info <- sample_info[colnames(filtered_counts), ]

combat_corrected_data <- ComBat_seq(
  filtered_counts,
  batch = filtered_sample_info$batch,
  group = filtered_sample_info$condition
)
combat_corrected_data <- round(combat_corrected_data)

# --- PCA dopo ComBat ---
run_pca(combat_corrected_data, filtered_sample_info, "ComBat corrected data")

# Parametri filtro DEGs
padj_cutoff <- 0.05
lfc_cutoff <- 1

# Strategia A: DESeq2 con batch
dds_model <- DESeqDataSetFromMatrix(
  countData = all_raw_counts,
  colData = sample_info,
  design = ~ batch + condition
)
dds_model <- DESeq(dds_model)
# Estrai dati normalizzati con batch incluso nel modello
vsd_model <- vst(dds_model, blind = FALSE)
vsd_mat <- assay(vsd_model)
# PCA su questi dati
run_pca(vsd_mat, sample_info, "DESeq2 con batch")

res_infected_model <- results(dds_model, contrast = c("condition", "infected", "mock"))
res_bystander_model <- results(dds_model, contrast = c("condition", "bystander", "mock"))

de_infected_model <- rownames(res_infected_model)[which(res_infected_model$padj < padj_cutoff &
                                                          abs(res_infected_model$log2FoldChange) > lfc_cutoff)]
de_bystander_model <- rownames(res_bystander_model)[which(res_bystander_model$padj < padj_cutoff &
                                                            abs(res_bystander_model$log2FoldChange) > lfc_cutoff)]

# Strategia B: ComBat + DESeq2
dds_combat <- DESeqDataSetFromMatrix(
  countData = combat_corrected_data,
  colData = filtered_sample_info,
  design = ~ condition
)
dds_combat <- DESeq(dds_combat)

res_infected_combat <- results(dds_combat, contrast = c("condition", "infected", "mock"))
res_bystander_combat <- results(dds_combat, contrast = c("condition", "bystander", "mock"))

de_infected_combat <- rownames(res_infected_combat)[which(res_infected_combat$padj < padj_cutoff &
                                                            abs(res_infected_combat$log2FoldChange) > lfc_cutoff)]
de_bystander_combat <- rownames(res_bystander_combat)[which(res_bystander_combat$padj < padj_cutoff &
                                                              abs(res_bystander_combat$log2FoldChange) > lfc_cutoff)]

############################################################
# STEP3.2 - Confronto, falsi positivi e heatmap
############################################################

# Confronti
infected_overlap <- intersect(de_infected_model, de_infected_combat) #geni trovati da entrambi i metodi
infected_unique_model <- setdiff(de_infected_model, de_infected_combat) #i geni trovati solo da model
infected_unique_combat <- setdiff(de_infected_combat, de_infected_model) #geni trovati solo da combat

bystander_overlap <- intersect(de_bystander_model, de_bystander_combat) #stessa cosa per il confronto bystander
bystander_unique_model <- setdiff(de_bystander_model, de_bystander_combat)
bystander_unique_combat <- setdiff(de_bystander_combat, de_bystander_model)

#VennDiagrams
# Venn plot 1: Infected vs Mock
venn.diagram(
  x = list(Model = de_infected_model, ComBat = de_infected_combat),
  filename = file.path(save_dir, "venn_infected_vs_mock.png"),
  imagetype = "png",   # <--- forza solo il PNG
  fill = c("lightblue", "lightgreen"),
  main = "DEGs Infected vs Mock (Model vs ComBat")


# Venn plot 2: Bystander vs Mock
venn.diagram(
  x = list(Model = de_bystander_model, ComBat = de_bystander_combat),
  filename = file.path(save_dir, "venn_bystander_vs_mock.png"),
  imagetype = "png",  
  fill = c("lightpink", "lightyellow"),
  main = "DEGs Bystander vs Mock (Model vs ComBat)"
)


# Heatmap function
plot_heatmap <- function(gene_list, count_matrix, sample_info, title) {
  if (length(gene_list) > 1) {
    # log2 transform
    mat <- log2(count_matrix[gene_list, , drop = FALSE] + 1)
    
    # controlla che sample_info e count_matrix abbiano sample_id coerenti
    rownames(sample_info) <- sample_info$sample_id
    sample_info <- sample_info[colnames(mat), , drop = FALSE]
    
    # forza l'ordine dei livelli
    sample_info$condition <- factor(sample_info$condition,
                                    levels = c("mock", "bystander", "infected"))
    sample_info$batch <- factor(sample_info$batch,
                                levels = c("GSE244485", "GSE244486"))
    
    # ordina prima per condition, poi per batch
    sample_info <- sample_info[order(sample_info$condition, sample_info$batch), ]
    mat <- mat[, rownames(sample_info), drop = FALSE]
    
    # calcola gaps_col per disegnare le linee divisorie tra condizioni
    cond_tab <- as.integer(table(sample_info$condition))
    gaps_col <- head(cumsum(cond_tab), -1)
    
    # definisci colori ben distinti
    ann_colors <- list(
      condition = c(mock = "#E41A1C",       # rosso
                    bystander = "#377EB8",  # blu
                    infected = "#4DAF4A"),  # verde
      batch = c(GSE244485 = "#984EA3",      # viola
                GSE244486 = "#FF7F00")      # arancione
    )
    
    # plot
    pheatmap(mat,
             scale = "row",
             annotation_col = sample_info[, c("condition", "batch"), drop = FALSE],
             annotation_colors = ann_colors,
             show_rownames = FALSE,
             cluster_cols = FALSE,    # mantieni ordine forzato
             gaps_col = gaps_col,     # linee di separazione tra condizioni
             main = title,
             angle_col = "45")        # nota: stringa e non numero
  } else {
    message("No unique genes for: ", title)
  }
}

# Heatmap per condizioni
# Salva l'heatmap per i geni differenzialmente espressi (DEGs) del modello DESeq2
png(file.path(save_dir, "heatmap_infected_model_DEGs.png"), width = 1000, height = 800)
plot_heatmap(infected_unique_model, all_raw_counts, sample_info, "Unique Model DEGs - Infected vs Mock")
dev.off()
# Salva l'heatmap per i geni DEGs corretti con ComBat
png(file.path(save_dir, "heatmap_infected_combat_DEGs.png"), width = 1000, height = 800)
plot_heatmap(infected_unique_combat, combat_corrected_data, filtered_sample_info, "Unique ComBat DEGs - Infected vs Mock")
dev.off()
# Salva l'heatmap per i geni DEGs del modello DESeq2
png(file.path(save_dir, "heatmap_bystander_model_DEGs.png"), width = 1000, height = 800)
plot_heatmap(bystander_unique_model, all_raw_counts, sample_info, "Unique Model DEGs - Bystander vs Mock")
dev.off()
# Salva l'heatmap per i geni DEGs corretti con ComBat
png(file.path(save_dir, "heatmap_bystander_combat_DEGs.png"), width = 1000, height = 800)
plot_heatmap(bystander_unique_combat, combat_corrected_data, filtered_sample_info, "Unique ComBat DEGs - Bystander vs Mock")
dev.off()


############################################################
# STEP4 - Dendrogramma dei campioni basato sui DEGs ComBat
############################################################
# Seleziona tutti i DEGs trovati con ComBat
deg_all <- unique(c(de_infected_combat, de_bystander_combat))
length(deg_all)

# Estrai matrice log2 e calcola z-score
log_mat <- log2(combat_corrected_data[deg_all, ] + 1)
mat_z <- t(scale(t(log_mat)))
mat_z[is.na(mat_z)] <- 0

# ==================================================
# Dendrogramma semplice dei campioni con tutti i DEGs
# ==================================================
mat_samples_all <- t(mat_z)
d <- dist(mat_samples_all, method = "euclidean")
hc <- hclust(d, method = "ward.D2")
png(file.path(save_dir, "dendrogram_top_DEGs_ComBat.png"), width = 1000, height = 800)
plot(hc,
     main = "Dendrogram of samples (top DEGs ComBat)",
     xlab = "", sub = "",
     cex = 0.8,
     hang = -1)
dev.off()

# ==================================================
# Top 500 geni più variabili
# ==================================================
gene_var <- apply(log_mat, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:500]
mat_z_small <- mat_z[top_genes, ]

# ==================================================
# Heatmap dei top 500 geni
# ==================================================
# Prepara annotazioni campioni
anno_col <- filtered_sample_info[, c("condition", "batch")]
anno_col$condition <- factor(anno_col$condition, levels = c("mock", "bystander", "infected"))
anno_col$batch     <- factor(anno_col$batch, levels = c("GSE244485", "GSE244486"))

# Colori annotazioni (ora definiti)
ann_colors <- list(
  condition = c("mock" = "#E41A1C", "bystander" = "#377EB8", "infected" = "#4DAF4A"),
  batch     = c("GSE244485" = "#984EA3", "GSE244486" = "#FF7F00")
)

# Ordina colonne per condition e batch
anno_col <- anno_col[order(anno_col$condition, anno_col$batch), ]
mat_heat <- mat_z_small[, rownames(anno_col), drop = FALSE]
cat("Dimensioni mat_z_small:", dim(mat_z_small), "\n")
cat("Dimensioni mat_heat:", dim(mat_heat), "\n")

# Genera heatmap

png(file.path(save_dir, "heatmap_top500_DEGs.png"), width = 1200, height = 1000, res = 150)
pheatmap(mat_heat,
         scale = "row",
         annotation_col = anno_col,
         annotation_colors = ann_colors,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         main = "Heatmap of top 500 DE genes (clustered genes, ordered samples)",
         fontsize = 8)
dev.off()


# ==================================================
# Confronto metriche diverse + ARI
# ==================================================
mat_samples <- t(mat_z_small)
common_samples <- intersect(rownames(mat_samples), rownames(filtered_sample_info))
mat_samples <- mat_samples[common_samples, , drop=FALSE]
sample_info_sub <- filtered_sample_info[common_samples, , drop=FALSE]
stopifnot(nrow(mat_samples) == nrow(sample_info_sub))

k <- length(levels(sample_info_sub$condition))
metrics <- c("euclidean","manhattan","pearson","spearman")
linkages <- c("complete","average","ward.D2")
ari_table <- data.frame()

for (m in metrics) {
  for (l in linkages) {
    if (m %in% c("pearson", "spearman")) {
      d <- as.dist(1 - cor(t(mat_samples), method = m))
    } else {
      d <- dist(mat_samples, method = m)
    }
    hc <- hclust(d, method = l)
    clusters <- cutree(hc, k)
    
    cond_vector <- as.integer(factor(sample_info_sub$condition,
                                     levels = levels(sample_info_sub$condition)))
    names(cond_vector) <- rownames(sample_info_sub)
    clusters <- clusters[names(cond_vector)]
    
    ari <- adjustedRandIndex(clusters, cond_vector)
    ari_table <- rbind(ari_table, data.frame(metric=m, linkage=l, ARI=ari))
  }
}

ari_table$metric <- factor(ari_table$metric, levels = c("euclidean","manhattan","pearson","spearman"))
ari_table$linkage <- factor(ari_table$linkage, levels = c("complete","average","ward.D2"))
ari_table$metric_link <- paste(ari_table$metric, ari_table$linkage, sep="_")

write.csv(ari_table, file.path(save_dir, "ARI_comparison_500genes_samples.csv"), row.names = FALSE)

# ==================================================
# Boxplot ARI per metrica
# ==================================================
png(file.path(save_dir, "ARI_boxplot_500genes_samples.png"),
    width = 1200, height = 800, res = 150)
boxplot(ARI ~ metric, data = ari_table,
        main = "ARI across distance metrics (top 500 genes, samples)",
        ylab = "Adjusted Rand Index",
        xlab = "Metric")
dev.off()

# ==================================================
# Boxplot combinato metric × linkage
# ==================================================
png(file.path(save_dir, "ARI_boxplot_500genes_samples_combined_horiz.png"),
    width = 3000, height = 1200, res = 150)
par(mar = c(5, 10, 4, 4))
boxplot(ARI ~ metric_link, data = ari_table,
        horizontal = TRUE,
        las = 1,
        cex.axis = 1,
        ylab = "",
        xlab = "Adjusted Rand Index",
        main = "ARI across metric and linkage (top 500 genes, samples)")
dev.off()

############################################################
# STEP 5 - Coding vs Noncoding DE genes (±5kb)
############################################################

# --- Setup directory ---
pairs_file <- file.path(save_dir, "DE_coding_noncoding_pairs_investigation.csv")

# --- 1. Uniforma nomi dei geni ---
mcols(gene_gr)$gene_id <- sub("\\..*", "", mcols(gene_gr)$gene_id)
rownames(combat_corrected_data) <- sub("\\..*", "", rownames(combat_corrected_data))

# --- 2. DE coding e noncoding genes ---
coding_genes    <- mcols(gene_gr)$gene_id[mcols(gene_gr)$biotype == "protein_coding"]
noncoding_genes <- mcols(gene_gr)$gene_id[mcols(gene_gr)$biotype != "protein_coding"]

deg_all <- sub("\\..*", "", deg_all)

res_coding    <- deg_all[deg_all %in% coding_genes]
res_noncoding <- deg_all[deg_all %in% noncoding_genes]

cat("INVESTIGATION: #DE coding =", length(res_coding), 
    "; #DE noncoding =", length(res_noncoding), "\n")

# --- 3. Subset GenomicRanges DEGs ---
de_coding_gr    <- gene_gr[mcols(gene_gr)$gene_id %in% res_coding]
de_noncoding_gr <- gene_gr[mcols(gene_gr)$gene_id %in% res_noncoding]
names(de_noncoding_gr) <- mcols(de_noncoding_gr)$gene_id

# --- 4. Espandi coding genes ±5kb ---
flank <- 5000
coding_flank <- resize(de_coding_gr, width = width(de_coding_gr) + 2*flank, fix = "center")
names(coding_flank) <- mcols(coding_flank)$gene_id

# --- 5. Trova overlaps entro ±5kb ---
hits <- findOverlaps(de_noncoding_gr, coding_flank, ignore.strand = TRUE)

hits_df <- data.frame(
  noncoding_idx   = queryHits(hits),
  coding_idx      = subjectHits(hits),
  noncoding_id    = names(de_noncoding_gr)[queryHits(hits)],
  coding_id       = names(coding_flank)[subjectHits(hits)],
  noncoding_strand= as.character(strand(de_noncoding_gr)[queryHits(hits)]),
  coding_strand   = as.character(strand(coding_flank)[subjectHits(hits)]),
  noncoding_chr   = as.character(seqnames(de_noncoding_gr))[queryHits(hits)],
  noncoding_start = start(de_noncoding_gr)[queryHits(hits)],
  noncoding_end   = end(de_noncoding_gr)[queryHits(hits)],
  coding_chr      = as.character(seqnames(coding_flank))[subjectHits(hits)],
  coding_start    = start(coding_flank)[subjectHits(hits)],
  coding_end      = end(coding_flank)[subjectHits(hits)]
)

# --- 6. Filtra per strand opposti (se disponibili) ---
hits_df$strand_check <- hits_df$noncoding_strand != hits_df$coding_strand
hits_opposite <- hits_df[hits_df$strand_check, ]

if (nrow(hits_opposite) > 0) {
  cat("Trovate", nrow(hits_opposite), "coppie coding-noncoding entro ±5kb con strand opposti\n")
  hits_to_use <- hits_opposite
} else {
  cat("Nessuna coppia entro ±5kb con strand opposti ??? useremo tutte le coppie entro ±5kb\n")
  hits_to_use <- hits_df
}

# --- 7. Tabella coppie finali ---
pairs_df <- unique(data.frame(
  noncoding_id     = hits_to_use$noncoding_id,
  coding_id        = hits_to_use$coding_id,
  noncoding_chr    = hits_to_use$noncoding_chr,
  noncoding_start  = hits_to_use$noncoding_start,
  noncoding_end    = hits_to_use$noncoding_end,
  coding_chr       = hits_to_use$coding_chr,
  coding_start     = hits_to_use$coding_start,
  coding_end       = hits_to_use$coding_end,
  noncoding_strand = hits_to_use$noncoding_strand,
  coding_strand    = hits_to_use$coding_strand,
  stringsAsFactors = FALSE
))

# Salva tabella su file
write.csv(pairs_df, pairs_file, row.names = FALSE)
cat("Tabella coppie salvata in:", pairs_file, "\n")

# --- 8. Violin plots ---
if (nrow(pairs_df) > 0) {
  all_coding    <- unique(pairs_df$coding_id)
  all_noncoding <- unique(pairs_df$noncoding_id)
  
  all_coding    <- all_coding[all_coding %in% rownames(combat_corrected_data)]
  all_noncoding <- all_noncoding[all_noncoding %in% rownames(combat_corrected_data)]
  
  ## --- 8.1 Media per campione ---
  expr_coding <- data.frame(
    sample     = colnames(combat_corrected_data),
    condition  = filtered_sample_info$condition,
    gene_type  = "coding",
    expression = as.numeric(colMeans(combat_corrected_data[all_coding, , drop=FALSE]))
  )
  
  expr_noncoding <- data.frame(
    sample     = colnames(combat_corrected_data),
    condition  = filtered_sample_info$condition,
    gene_type  = "noncoding",
    expression = as.numeric(colMeans(combat_corrected_data[all_noncoding, , drop=FALSE]))
  )
  
  expr_mean <- rbind(expr_coding, expr_noncoding)
  
  p_mean <- ggplot(expr_mean, aes(x=condition, y=expression, fill=condition)) +
    geom_violin(trim=FALSE, scale="width") +
    geom_jitter(width=0.2, size=1, alpha=0.6) +
    facet_wrap(~gene_type, scales="free_y") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    ggtitle("Coding vs Noncoding DE genes (±5kb) - Media per campione")
  
  ggsave(file.path(save_dir, "violinplot_coding_vs_noncoding_means.png"),
         p_mean, width=8, height=6, dpi=300)
  
  ## --- 8.2 Tutti i geni ---
  expr_coding_full <- data.frame(
    sample     = rep(colnames(combat_corrected_data), each=length(all_coding)),
    condition  = rep(filtered_sample_info$condition, each=length(all_coding)),
    gene_type  = "coding",
    expression = as.vector(as.matrix(combat_corrected_data[all_coding, , drop=FALSE]))
  )
  
  expr_noncoding_full <- data.frame(
    sample     = rep(colnames(combat_corrected_data), each=length(all_noncoding)),
    condition  = rep(filtered_sample_info$condition, each=length(all_noncoding)),
    gene_type  = "noncoding",
    expression = as.vector(as.matrix(combat_corrected_data[all_noncoding, , drop=FALSE]))
  )
  
  expr_full <- rbind(expr_coding_full, expr_noncoding_full)
  
  p_full <- ggplot(expr_full, aes(x=condition, y=expression, fill=condition)) +
    geom_violin(trim=FALSE, scale="width") +
    theme_bw() +
    facet_wrap(~gene_type, scales="free_y") +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    ggtitle("Coding vs Noncoding DE genes (±5kb) - Tutti i geni")
  
  ggsave(file.path(save_dir, "violinplot_coding_vs_noncoding_melt.png"),
         p_full, width=8, height=6, dpi=300)
  
  cat("Entrambi i plot salvati in:", save_dir, "\n")
} else {
  cat("Nessuna coppia disponibile per plot\n")
}


