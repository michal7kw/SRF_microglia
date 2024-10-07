## Balanced SET levels favor the correct enhancer repertoire during cell fate acquisition

[Balanced SET levels favor the correct enhancer repertoire during cell fate acquisition | Nature Communications](https://www.nature.com/articles/s41467-023-39043-x)

## Description and Specification:

This study contains RNA sequencing data from human neural progenitor cells (NPCs) and cerebral organoids modeling Schinzel-Giedion syndrome (SGS). Key details include:

1. Sample types:

   - iPSC-derived NPCs
   - Cerebral organoids
2. Experimental conditions:

   - SETBP1 wildtype (D868D, I871I)
   - SETBP1 SGS mutations (D868N, I871T)
3. Sample details:

   - 9 NPC samples (3 per genotype: D868D, D868N, I871T)
   - 2 cerebral organoid samples (WT D868D vs Mut D868N)
4. Sequencing details:

   - NPCs: Bulk RNA-seq on Illumina HiSeq 3000
   - Organoids: Single-cell RNA-seq on Illumina NovaSeq 6000
5. Key analyses:

   - Differential gene expression between wildtype and SGS mutant cells
   - Functional enrichment analysis
   - Single-cell transcriptomics of cerebral organoids
6. Data processing:

   - Bulk RNA-seq: STAR aligner, DESeq2 for differential expression
   - scRNA-seq: Cell Ranger pipeline, Seurat for analysis
7. Genome build: hg38

This dataset allows for analysis of transcriptional changes in neural progenitors and developing cortical cells due to SGS-associated SETBP1 mutations. The combination of bulk and single-cell approaches provides insights into both overall expression changes and cell type-specific effects in SGS.
