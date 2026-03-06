#Modul: Profil Ekspresi Gen Unit Neurovaskular pada Diabetes Melitus Tipe 2
#Dataset: GSE161355 (DIabetes Tipe 2 vs Normal)
#Platform: Microarray (Affymetrix Human Genome U133 Plus 2.0 Array - GPL570)
#Tujuan: Mengidentifikasi Differentially Expressed Genes (DEG) 

#PART A. PENGANTAR KONSEP 

#Analisis ekspresi gen bertujuan untuk membandingkan tingkat ekspresi gen 
#antara dua kondisi biologis (misalnya kanker vs normal) 
#Pada modul ini kita menggunakan pendekatan statistik limma (Linear Models
#for Microarray Data), yang merupakan standar emas untuk data microarray. 

#PART B. PERSIAPAN LINGKUNGAN KERJA (INSTALL & LOAD PACKAGE) 

#Apa itu package? 
#Package adalah kumpulan fungsi siap pakai di R
#Bioinformatika di R sangat bergantung pada package dari CRAN dan Bioconductor 

#1. Install BiocManager (manajer paket Bioconductor) 
#IF adalah struktur logika : “jika kondisi terpenuhi, lakukan aksi”

if (!require("BiocManager", quietly = TRUE))  {
  install.packages("BiocManager") 
}

# 2. Install paket Bioconductor (GEOquery & limma) 
#GEOquery: mengambil data dari database GEO 
#limma: analisis statistik ekspresi gen 

BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update = FALSE) 

#Install annotation package sesuai platform
#GPL570 = Affymetrix Human Genome U133 Plus 2.0 Array 
BiocManager::install("hgu133plus2.db", ask = FALSE, update = FALSE)

#3. Install paket CRAN untuk visualisasi dan manipulasi data 
#phetmap: heatmap ekspresi gen 
#ggplot2: grafik (volcano plot)
#dplyr: manipulasi tabel data 

install.packages(c("pheatmap", "ggplot2", "dplyr"))

#umap: grafik (plot UMAP) 
if (!requireNamespace("umap", quietly = TRUE)) {
  install.packages("umap")
}

#4. Memanggil library 
#library() digunakan agar fungsi di dalam package bisa digunakan 
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(hgu133plus2.db)
library(AnnotationDbi)
library(umap)

------------------------------------------------------------------------------------

#PART C. PENGAMBILAN DATA DARI GEO 


#GEO (Gene Expression Omnibus) adalah database publik milik NCBI
#getGEO(): fungsi untuk mengunduh dataset berdasarkan ID GEO
#GSEMatrix = TRUE -> data diambil dalam format ExpressionSet
#AnnotGPL  = TRUE -> anotasi gen (Gene Symbol) ikut diunduh

gset <- getGEO("GSE161355", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

#ExpressionSet berisi:
# - exprs() : matriks ekspresi gen
# - pData() : metadata sampel
# - fData() : metadata fitur (probe / gen)


--------------------------------------------------------------------------------
  

#PART D. PRE-PROCESSING DATA EKSPRESI 

# exprs(): mengambil matriks ekspresi gen
# Baris  = probe/gen
# Kolom  = sampel
ex <- exprs(gset)

#Mengapa perlu log2 transformasi?
#Data microarray mentah memiliki rentang nilai sangat besar.
#Log2 digunakan untuk:
#1. Menstabilkan varians
#2. Mendekati asumsi model linear
#3. Memudahkan interpretasi log fold change

#quantile(): menghitung nilai kuantil (persentil)
#as.numeric(): mengubah hasil quantile (yang berupa named vector)
#menjadi vektor numerik biasa agar mudah dibandingkan
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))

#LogTransform adalah variabel logika (TRUE / FALSE)
#Operator logika:
#>  : lebih besar dari
#|| : OR (atau)
#&& : AND (dan)
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

#IF statement:
#Jika LogTransform = TRUE, maka lakukan log2
#Nilai <= 0 tidak boleh di-log, maka diubah menjadi NA

if (LogTransform) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}


---------------------------------------------------------------------------------

#PART E. DEFINISI KELOMPOK SAMPEL 

#Dari ExpressionSet, perlu tahu mana kolom spesifik mana yang berisikan info biologis
View(pData(gset))
unique(pData(gset))["source_name_ch1"]

#pData(): metadata sampel
#source_name_ch1 berisi informasi kondisi biologis sampel
group_info <- pData(gset)[["source_name_ch1"]]

#make.names(): mengubah teks menjadi format valid untuk R
groups <- make.names(group_info)

#factor():
#Mengubah data kategorik menjadi faktor
#Faktor sangat penting untuk analisis statistik di R
gset$group <- factor(groups)

#levels(): melihat kategori unik dalam faktor
nama_grup <- levels(gset$group)
print(nama_grup)


-----------------------------------------------------------------------------------


# PART F. FILTERING SAMPEL (OPSI 2: NEURON & ENDOTEL)

# Kita hanya ingin sampel Neuron dan Endotel (Abaikan untuk Astrosit)
# Kita buat vektor nama grup yang ingin dianalisis (Neuron & Endotel)
grup_fokus <- c("Control.Endothelial.cells", "Diabetic.Endothelial.cells", 
                "Control.Neurones", "Diabetic.Neurones")

# Filter ExpressionSet agar hanya berisi sampel dari 4 grup tersebut
gset_filtered <- gset[, gset$group %in% grup_fokus]

# Refresh faktor agar level "Astrocytes" hilang dari memori
gset_filtered$group <- factor(gset_filtered$group, levels = grup_fokus)


--------------------------------------------------------------------------------


#  PART G. DESIGN MATRIX 

# model.matrix():
# Membuat matrix yang memetakan setiap sampel ke grupnya
# Mengubah data kategori (teks kelompok) menjadi data numerik (angka biner; 0 dan 1)
# ~0 berarti TANPA intercept (best practice LIMMA)
design <- model.matrix(~0 + gset_filtered$group)

# colnames(): pemberian nama kolom agar mudah dibaca
# Memberi nama kolom angka biner menggunakan nama asli kelompok
colnames(design) <- levels(gset_filtered$group)


---------------------------------------------------------------------------------


# PART H. MEMBUAT KONTRAS (LOGIKA PERBANDINGAN) 

# Membuat matriks perbandingan
# Di sini kita tentukan apa yang mau kita bandingkan
cont.matrix <- makeContrasts(
  # 1. Dampak Diabetes pada Endotel (Pembuluh Darah)
  Diabetes_Efek_Endotel = Diabetic.Endothelial.cells - Control.Endothelial.cells,
  
  # 2. Dampak Diabetes pada Neuron (Sel Saraf)
  Diabetes_Efek_Neuron = Diabetic.Neurones - Control.Neurones,
  
  levels = design
)


---------------------------------------------------------------------------------


# PART I. LINEAR MODEL & EBAYES

# lmFit():
# Membangun model linear untuk setiap gen
# Menghitung rata-rata ekspresi untuk setiap gen pada tiap grup
fit <- lmFit(gset_filtered, design)


# contrasts.fit(): Menerapkan kontras perbandingan ke model
# Tahapan eksekusi
fit2 <- contrasts.fit(fit, cont.matrix)


# eBayes(): 
# Uji validitas statistik dengan menggunakan eBayes 
# Empirical Bayes untuk menstabilkan estimasi variasi
# Untuk mendapatkan perhitungan p-value yang lebih akurat
fit2 <- eBayes(fit2)


# topTable()
# Menampilkan hasil akhir DEG dalam bentuk table
# Tabel gen yang berbeda signifikan pada Endotel
res_endotel <- topTable(
  fit2, 
  coef="Diabetes_Efek_Endotel", 
  adjust="fdr", 
  sort.by="B", 
  number=Inf 
)

head(res_endotel)

# Tabel gen yang berbeda signifikan pada Neuron
res_neuron <- topTable(
  fit2,
  coef="Diabetes_Efek_Neuron", 
  adjust="fdr", 
  sort.by="B", 
  number=Inf
)

head(res_neuron)

# Kita ambil gen yang P.Value-nya (mentah) di bawah 0.05
# Ini dilakukan untuk kedua sel (Endotel dan Neuron)
res_endotel_sig <- res_endotel[res_endotel$P.Value < 0.05, ]
res_neuron_sig <- res_neuron[res_neuron$P.Value < 0.05, ]

# Cek berapa banyak gen yang "terpilih"
nrow(res_endotel_sig)
nrow(res_neuron_sig)

# --- CATATAN ANALISIS: PEMILIHAN AMBANG BATAS SIGNIFIKANSI ---
# Pada dataset ini (GSE161355), jumlah sampel per kelompok terbatas (n=6).
# Penggunaan koreksi Multiple Testing (FDR/adj.P.Val) menghasilkan nilai yang 
# sangat konservatif (mendekati 1), sehingga tidak ada gen yang lolos filter ketat.
#
# Untuk tujuan eksplorasi biologis dan pembuatan profil ekspresi (Heatmap/Volcano), 
# saya menggunakan 'Raw P-Value < 0.05' sebagai ambang batas. 
# Strategi ini umum dilakukan pada studi awal (discovery phase) dengan ukuran 
# sampel kecil untuk mengidentifikasi tren perubahan ekspresi yang potensial 
# sebelum divalidasi lebih lanjut.


---------------------------------------------------------------------------------


#PART J. ANOTASI NAMA GEN 

#Penting:
#Pada data microarray Affymetrix, unit analisis awal adalah PROBE,
#bukan gen. Oleh karena itu, anotasi ulang diperlukan menggunakan
#database resmi Bioconductor.


# --- PROSES ANOTASI UNTUK ENDOTEL ---
  
# Mengambil ID probe dari hasil signifikan Endotel
probe_ids_endotel <- rownames(res_endotel_sig)

# Mapping probe -> gene symbol & gene name
# Gunakan library anotasi yang sesuai dengan platform (Affymetrix U133 Plus 2.0)
# Biasanya menggunakan package 'hgu133plus2.db'
gene_annotation_endotel <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = probe_ids_endotel,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)

# Gabungkan dengan hasil limma (Endotel)
res_endotel_sig$PROBEID <- rownames(res_endotel_sig)
res_endotel_sig <- merge(
  res_endotel_sig,
  gene_annotation_endotel,
  by = "PROBEID",
  all.x = TRUE
)

# --- PROSES ANOTASI UNTUK NEURON ---

# Mengambil ID probe dari hasil signifikan Neuron
probe_ids_neuron <- rownames(res_neuron_sig)

# Mapping probe -> gene symbol & gene name
# Gunakan library anotasi yang sesuai dengan platform (Affymetrix U133 Plus 2.0)
# Biasanya menggunakan package 'hgu133plus2.db'
gene_annotation_neuron <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = probe_ids_neuron,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)

# Gabungkan dengan hasil limma (Neuron)
res_neuron_sig$PROBEID <- rownames(res_neuron_sig)
res_neuron_sig <- merge(
  res_neuron_sig,
  gene_annotation_neuron,
  by = "PROBEID",
  all.x = TRUE
)

# --- CEK HASIL ANOTASI ---
print("Hasil Anotasi Endotel:")
head(res_endotel_sig[, c("PROBEID", "SYMBOL", "GENENAME")])

print("Hasil Anotasi Neuron:")
head(res_neuron_sig[, c("PROBEID", "SYMBOL", "GENENAME")])
  

----------------------------------------------------------------------------------  


#PART K.1 BOXPLOT DISTRIBUSI NILAI EKSPRESI 

# Boxplot digunakan untuk:
# - Mengecek distribusi nilai ekspresi antar sampel
# - Melihat apakah ada batch effect
# - Mengevaluasi apakah normalisasi/log-transform sudah wajar

# 1. Atur margin (bawah, kiri, atas, kanan)
par(mar = c(8, 4, 4, 12), xpd = TRUE)
  
# 2. Set warna agar lebih kontras antar grup
# Kita gunakan palette warna yang lebih soft tapi jelas
group_colors <- as.numeric(gset_filtered$group) + 1


# Mengambil matriks ekspresi dari data yang sudah difilter
ex_filtered <- exprs(gset_filtered)

boxplot(
  ex_filtered,
  col = group_colors,
  las = 2,
  outline = FALSE,
  cex.axis = 0.7,
  main = "Boxplot Distribusi Nilai Ekspresi per Sampel (Endotel & Neuron)",
  ylab = "Expression Value (log2)"
)

# 3. Letakkan legenda di luar area plot (inset koordinat x digunakan untuk geser ke kanan)
# Kita hitung jumlah sampel untuk menentukan titik X nya
posisi_x <- ncol(ex_filtered) + 2

legend(
  posisi_x, 14,
  legend = levels(gset_filtered$group),
  fill = unique(group_colors),
  cex = 0.6,
  bty = "n"
)

# Reset margin ke default setelah selesai, agar tidak mempengaruhi plot berikutnya
par(mar=c(5, 4, 4, 2), xpd=FALSE)

------------------------------------------------------------------------------------

#PART K.2 DISTRIBUSI NILAI EKSPRESI (DENSITY PLOT) 

#Density plot menunjukkan sebaran global nilai ekspresi gen
#Digunakan untuk:
#- Mengecek efek log-transform
#- Membandingkan distribusi antar grup

# Pastikan data frame dibuat dari data yang sudah difilter
expr_long <- data.frame(
  Expression = as.vector(ex_filtered),
  Group = rep(gset_filtered$group, each = nrow(ex_filtered))
)

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") + 
  labs(
    title = "Distribusi Global Ekspresi Gen (Endotel vs Neuron)",
    subtitle = "Dataset: GSE161355",
    x = "Expression Value (log2)",
    y = "Density"
  )

------------------------------------------------------------------------------------

#PART K.3 UMAP (VISUALISASI DIMENSI RENDAH)

# Catatan Analisis: 
# UMAP digunakan untuk melihat pengelompokan (clustering) sampel secara global.
# Kita ingin memastikan bahwa profil genetik sel Endotel dan Neuron memang 
# terpisah secara alami, serta melihat seberapa kuat pengaruh Diabetes 
# terhadap perubahan ekspresi gen di masing-masing tipe sel.
  

# 1. Persiapan Data
# Transpose matriks: UMAP membutuhkan baris sebagai Sampel dan kolom sebagai Gen.
# ex_filtered berisikan data ekspresi dari sampel Endotel & Neuron.
umap_input <- t(ex_filtered)

# 2. Eksekusi Algoritma UMAP
# set.seed(123) memastikan hasil koordinat titik selalu sama setiap kali di-run.
set.seed(123) 
umap_result <- umap(umap_input)

# 3. Re-formatting Hasil untuk Visualisasi
# Menyimpan koordinat UMAP1 dan UMAP2 ke dalam satu tabel (data frame)
# agar bisa diolah dengan ggplot2.
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset_filtered$group  # Menambahkan informasi label grup (4 kategori)
)

# 4. Visualisasi Scatter Plot UMAP
# Interpretasi: 
# - Jarak antar titik menunjukkan kemiripan profil ekspresi gen.
# - Kelompok sel yang berbeda (Endotel vs Neuron) idealnya terpisah jauh.
# - Sampel dalam grup yang sama (misal: Diabetic vs Control) diharapkan mengelompok.

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "UMAP Plot: Pengelompokan Sampel Global",
    subtitle = "Dataset: GSE161355 | Perbandingan Sel Endotel vs Neuron",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

# --- KESIMPULAN VISUAL ---
# Jika titik merah/biru (Endotel) terpisah jauh dari titik hijau/ungu (Neuron),
# itu menunjukkan identitas seluler yang terjaga dalam dataset ini.


--------------------------------------------------------------------------------------


#PART L.1 VISUALISASI VOLCANO PLOT 

# Volcano Plot menggabungkan Log Fold Change (sumbu X) dan Signifikansi Statistik (sumbu Y).
# Penggunaan Raw P-Value karena jumlah sampel yang kecil (n=5/6)
# - Gen di kanan atas (Merah): Gen yang naik signifikan (UP-regulated)
# - Gen di kiri atas (Biru): Gen yang turun signifikan (DOWN-regulated)
# - Gen di bawah (Abu-abu): Gen yang tidak berubah secara signifikan.
  

# 1. VOLCANO PLOT: SEL ENDOTEL (Diabetic vs Control)

# STEP 1: Persiapan Data Frame
# Mengekstrak kolom logFC, P.Value, dan SYMBOL dari hasil limma res_endotel
# Konversi ke Data Frame agar kolom mudah dipanggil
df_endotel <- as.data.frame(res_endotel)

# Persiapan Data (SESUAI NAMA KOLOM DATA)
v_data_endotel <- data.frame(
  logFC = df_endotel$logFC,
  P.Value = df_endotel$P.Value,
  Gene = df_endotel$Gene.symbol 
)
  
# STEP 2: Klasifikasi Gen Signifikan (Thresholding)
# Menandai gen berdasarkan ambang batas P < 0.05 dan |logFC| > 0.1
v_data_endotel$status <- "NO"
v_data_endotel$status[v_data_endotel$logFC > 0.1 & v_data_endotel$P.Value < 0.05] <- "UP"
v_data_endotel$status[v_data_endotel$logFC < -0.1 & v_data_endotel$P.Value < 0.05] <- "DOWN"


# STEP 3: Plotting dengan ggplot2
ggplot(v_data_endotel, aes(x = logFC, y = -log10(P.Value), color = status)) +
  geom_point(alpha = 0.4, size = 1.5) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(title = "Volcano Plot: Endothelial Cells",
       subtitle = "Threshold: P < 0.05 & |logFC| > 0.1",
       x = "Log2 Fold Change", y = "-log10 P-Value")



# 2. VOLCANO PLOT: SEL NEURON (Diabetic vs Control)

# STEP 1: Persiapan Data Frame
# Mengekstrak kolom logFC, P.Value, dan SYMBOL dari hasil limma res_neuron
# Konversi ke Data Frame agar kolom mudah dipanggil
df_neuron <- as.data.frame(res_neuron)

# Persiapan Data (SESUAI NAMA KOLOM DATA)
v_data_neuron <- data.frame(
  logFC = df_neuron$logFC,
  P.Value = df_neuron$P.Value,
  Gene = df_neuron$Gene.symbol 
)

# STEP 2: Klasifikasi Gen Signifikan (Thresholding)
v_data_neuron$status <- "NO"
v_data_neuron$status[v_data_neuron$logFC > 0.1 & v_data_neuron$P.Value < 0.05] <- "UP"
v_data_neuron$status[v_data_neuron$logFC < -0.1 & v_data_neuron$P.Value < 0.05] <- "DOWN"

# STEP 3: Plotting dengan ggplot2
ggplot(v_data_neuron, aes(x = logFC, y = -log10(P.Value), color = status)) +
  geom_point(alpha = 0.4, size = 1.5) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(title = "Volcano Plot: Neuronal Cells",
       subtitle = "Threshold: P < 0.05 & |logFC| > 0.1",
       x = "Log2 Fold Change", y = "-log10 P-Value")



-----------------------------------------------------------------------------------


# PART K: ANALISIS IRISAN GEN (VENN DIAGRAM)

# Catatan: Kita membandingkan simbol gen yang signifikan (P < 0.05) 
# antara sel Endotel dan sel Neuron untuk mencari biomarker universal.

# Install package jika belum ada
install.packages("VennDiagram")
library(VennDiagram)
grid.newpage()

# STEP 1: Ambil daftar simbol gen yang unik dari masing-masing grup
genes_endotel <- unique(res_endotel_sig$SYMBOL)
genes_neuron <- unique(res_neuron_sig$SYMBOL)

# STEP 2: Membuat Plot Venn
venn.plot <- draw.pairwise.venn(
  area1 = length(genes_endotel),
  area2 = length(genes_neuron),
  cross.area = length(intersect(genes_endotel, genes_neuron)),
  category = c("Endothelial", "Neuronal"),
  fill = c("skyblue", "pink"),
  lty = "blank",
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05)
)

# Menampilkan Venn di RStudio
grid.draw(venn.plot)

# Menambahkan judul di atasnya pakai grid.text
grid.text("Comparison of DEGs in Diabetic Endothelial vs Neuronal Cells", 
          y = unit(0.95, "npc"), # Mengatur posisi di atas (95% tinggi canvas)
          gp = gpar(fontsize = 16, fontface = "bold")) # Ukuran dan gaya tulisan


# STEP 3: Mencari tahu gen apa saja yang ada di irisan tersebut
common_genes <- intersect(genes_endotel, genes_neuron)

# Tampilkan gen yang sama di console
print("Gen yang beririsan di kedua sel:")
print(common_genes)  

# Hasil analisis ini memetakan dampak Diabetes (GSE161355):
# 1. Irisan (94 gen): Dampak Diabetes secara sistemik/umum pada jaringan otak.
# 2. Sisi Endotel (1018 gen): Kerusakan spesifik pada sistem vaskular otak.
# 3. Sisi Neuron (927 gen): Kerusakan spesifik pada sel saraf otak.
  

-----------------------------------------------------------------------------------


#PART L.3 VISUALISASI HEATMAP 

# Heatmap ini digunakan untuk menunjukkan pola ekspresi 94 gen yang 
# secara konsisten terpengaruh oleh diabetes baik di sel endotel maupun neuron.
  
# STEP 1: Identifikasi Simbol Gen yang Beririsan (94 Gen)
# Mengambil gen yang signifikan di kedua analisis (hasil intersect sebelumnya)
common_symbols <- intersect(res_endotel_sig$SYMBOL, res_neuron_sig$SYMBOL)


# STEP 2: Mengambil Probe ID dari Hasil Anotasi
# Karena ini irisan, ProbeID-nya pasti ada di kedua tabel (Endotel & Neuron)
idx <- which(res_endotel_sig$SYMBOL %in% common_symbols)
selected_probes <- res_endotel_sig$PROBEID[idx]
selected_symbols <- res_endotel_sig$SYMBOL[idx]


# STEP 3: Mengambil Matriks Ekspresi dari 'ex_filtered'
# ex_filtered berisikan data ekspresi dari sampel Endotel & Neuron.
mat_heatmap <- ex_filtered[selected_probes, ]


# STEP 4: Memberil nama baris menggunakan Gene Symbol (agar lebih informatif dibanding Probe ID)
# Jika Symbol kosong, tetap gunakan Probe ID sebagai cadangan
gene_label <- ifelse(
  is.na(selected_symbols) | selected_symbols == "",
  NA,                # Ganti ini jadi NA agar bisa difilter
  selected_symbols
)
rownames(mat_heatmap) <- gene_label
 ---------------------------
gene_label <- ifelse(
  is.na(selected_symbols) | selected_symbols == "",
  selected_probes,
  selected_symbols
)
rownames(mat_heatmap) <- gene_label


# STEP 5: Pembersihan data (Wajib agar tidak error saat Clustering)
# Hapus baris yang mengandung NA
mat_heatmap <- mat_heatmap[!is.na(rownames(mat_heatmap)), ]
mat_heatmap <- mat_heatmap[rowSums(is.na(mat_heatmap)) == 0, ]

# Hapus gen dengan varians nol (gen yang tidak berubah sama sekali di semua sampel)
gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]


# STEP 6: Anotasi kolom (Mengambil grup dari metadata gset)
# Di sini kita buat label untuk membedakan sampel Endotel vs Neuron 
# serta kondisi Diabetic vs Normal
annotation_col <- data.frame(Group = gset_filtered$group)
rownames(annotation_col) <- colnames(mat_heatmap)


# STEP 7: Visualisasi Heatmap
library(pheatmap)
pheatmap(
  mat_heatmap,
  scale = "row", 
  annotation_col = annotation_col,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  main = "Comparison of 94 Common DEGs in Diabetic Cells",
  color = colorRampPalette(c("blue", "white", "red"))(100)
)


--------------------------------------------------------------------------------------


#PART M. MENYIMPAN HASIL 

# 1. Menyimpan hasil lengkap analisis DEGs untuk Sel Endotel
write.csv(res_endotel_sig, "Hasil_DEGs_Endotel_Diabetes.csv", row.names = FALSE)

# 2. Menyimpan hasil lengkap analisis DEGs untuk Sel Neuron
write.csv(res_neuron_sig, "Hasil_DEGs_Neuron_Diabetes.csv", row.names = FALSE)

# 3. Menyimpan daftar 94 gen yang beririsan 
# Mengambil data dari res_endotel_sig yang masuk dalam daftar common_symbols
common_degs_final <- res_endotel_sig[res_endotel_sig$SYMBOL %in% common_symbols, ]
write.csv(common_degs_final, "Hasil_94_Common_DEGs_Diabetes.csv", row.names = FALSE)

# Memberikan pesan konfirmasi di console
message("Analisis selesai! Tiga file CSV hasil telah disimpan di folder kerja Anda.")
  
  
# write.csv(): menyimpan hasil analisis ke file CSV
write.csv(topTableResults, "Hasil_GSE10072_DEG.csv")

message("Analisis selesai. File hasil telah disimpan.")



