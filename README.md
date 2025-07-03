# scRNA-seq COVID vs Control Lung Analysis

![scRNA-seq](https://img.shields.io/badge/scRNA--seq-Analysis-blue)
![Status](https://img.shields.io/badge/status-active-brightgreen)

## 📌 Project Overview

This repository contains a **comprehensive single-cell RNA sequencing (scRNA-seq) analysis** comparing lung tissue from a COVID-19 patient with healthy control tissue. The analysis reveals **pathological cellular stress responses and immune dysfunction** characteristic of severe COVID-19, providing molecular insights into **why some patients develop life-threatening disease**.

Rather than showcasing protective immunity, this dataset uncovers **stress-dominated signatures, immune exhaustion, and epithelial dysfunction** that underlie severe COVID-19 pathology.

---

## 🧪 Dataset Summary

* High-quality scRNA-seq dataset with:

  * **\~20 distinct cell clusters** (respiratory epithelium, immune, stromal)
  * **Low mitochondrial gene expression (0-4%)**
  * **0–5000+ genes per cell**, **0–20,000+ molecules per cell**
* Comprehensive annotation:

  * Alveolar (AT1, AT2), airway epithelial, endothelial, fibroblasts
  * CD4+, CD8+ T cells, B cells, NK cells, macrophages, dendritic cells

---

## 🔍 Key Findings

### 1️⃣ Universal Cellular Stress

* **HSP90AA1** upregulated across all cell types
* Additional heat shock proteins (HSPA1A, HSPH1) prominent
* Indicates **crisis-level stress responses** over protective immunity

### 2️⃣ Immune Exhaustion Patterns

* CD4+ and CD8+ T cells: hyperactivation signatures without effective immunity
* B cells: minimal activation (immune suppression)
* Macrophages: inflammatory signatures without resolution

### 3️⃣ Epithelial and Structural Cell Dysfunction

* AT2 cells: disruption in surfactant production and repair
* Airway epithelium: compromised barrier genes
* Endothelium: signatures of vascular leak and damage

### 4️⃣ Pathway Insights

* **KEGG Enrichment** highlights:

  * Protein processing in ER
  * PI3K-Akt signaling
  * Antigen processing and presentation
* Pathways suggest **maladaptive stress responses** over immune protection

### Novel Contributions Highlight:
✅ Universal stress response as a marker for severe COVID-19.
✅ Cell-type-specific dysfunction hierarchy explaining severe disease outcomes.
✅ Identification of tissue-specific immune exhaustion in severe COVID-19 lungs.
✅ Data suggests immune failure, not hyperactivation, in severe COVID-19.

---

## ⚡ Technical Workflow

1. **Quality Control & Preprocessing**:

   * Filtering, normalization (SCTransform), mitochondrial filtering
   * Batch correction if applicable

2. **Clustering & Annotation**:

   * Unsupervised clustering (`Seurat`)
   * Marker gene identification
   * Manual annotation with canonical markers

3. **Differential Gene Expression**:

   * Per-cluster COVID vs Control DEGs (`FindMarkers`)
   * Volcano plots for each cell type

4. **Pathway Enrichment**:

   * KEGG enrichment on up/downregulated genes (`clusterProfiler`)
   * Dot plot and barplot visualization

5. **Visualization**:

   * UMAP plots
   * Volcano plots
   * KEGG pathway enrichment plots

---

## 📁 Repository Structure

```
.
├── data/                       # Raw or processed input data (large files may require LFS)
├── results/                    # Processed Seurat objects, DGE results
├── figures/                    # UMAPs, volcano plots, KEGG plots
├── scripts/                    # R scripts for each step (QC, clustering, DGE, KEGG)
└── README.md                   # Project overview and instructions
```

---

## 🚀 Usage

### Prerequisites:

* R >= 4.2
* Packages:

  * `Seurat`
  * `dplyr`, `ggplot2`
  * `clusterProfiler`, `org.Hs.eg.db`
  * `EnhancedVolcano`, `forcats`

### Run Analysis:

1. Clone the repository:

   ```bash
   git clone https://github.com/Maha6932/scRNAseq-COVID-vs-Control.git
   cd scRNAseq-COVID-vs-Control
   ```

2. Run analysis scripts in order:

   * `01_Data_Loading_and_intial_setup.R`
   * `02_Metadata_alignment_Quality_control.R`
   * `03_Normalization_Clustering_Marker_gene_identification.R`
   * `04_Enhanced_single_cell_RNA_seq_analysis_pipeline.R`
   * `05_marker_genes.R`
   * `06_DEGS_per_cluster.R`
   * `07_KEGG_Enrichment.R`

3. Outputs:

   * Processed Seurat objects in `results/`
   * Volcano and KEGG plots in `figures/`
   * Combined KEGG enrichment CSVs for analysis

---

## 📈 Key Results

* Universal **stress signatures** in severe COVID-19 lung tissue
* **Immune cell exhaustion patterns** rather than effective activation
* **Epithelial and endothelial dysfunction** consistent with respiratory failure
* **Pathological KEGG pathway enrichment** indicating maladaptive cellular responses

---

## ⚖️ Limitations

* Single patient-control comparison
* Lung-specific analysis
* Represents **severe COVID-19 pathology** rather than mild disease

---

### References

Mohammed, R.N., Tamjidifar, R., Rahman, H.S., et al. (2022). A comprehensive review about immune responses and exhaustion during coronavirus disease (COVID-19). Cell Communication and Signaling, 20, 79. https://doi.org/10.1186/s12964-022-00856-w

This paper contextualizes immune exhaustion in COVID-19, complementing our findings of stress-dominated responses and immune dysfunction at the single-cell level in lung tissue.

📌 This project complements Mohammed et al. (2022), which reviewed systemic immune exhaustion in COVID-19. By analyzing single-cell lung tissue, we demonstrate how severe COVID-19 involves cell-type-specific stress responses, immune dysfunction, and impaired repair, aligning with their findings but providing tissue-level, single-cell resolution insights. This connection strengthens our understanding of severe COVID-19 as a disease of immune exhaustion and cellular crisis, identifying potential biomarkers and therapeutic targets for intervention.

---

## ✏️ Citation

If you use this workflow or analysis structure, please cite:

> **Mahalakshmi Aranganathan. "Single-cell RNA-seq analysis of severe COVID-19 lung tissue reveals universal cellular stress and immune exhaustion patterns." GitHub, 2025.**
> [GitHub Repository](https://github.com/Maha6932/scRNAseq-COVID-vs-Control)

---

## 🤝 Contributions

Open to contributions including:

* Visualization utilities
* Workflow optimization
* Expansion to multi-sample or multi-condition analyses

---

## 📫 Contact

**Mahalakshmi Aranganathan**
Email: mahalakshmi.aranga@gmail.com
LinkedIn: https://www.linkedin.com/in/mahalakshmi-aranganathan/

---

## License

This project is licensed under the MIT License.

---

**This repository provides a clear, reusable pipeline for scRNA-seq analysis and insights into severe COVID-19 pathophysiology, supporting therapeutic exploration and academic learning.**

