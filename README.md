1-MergePeak_calldif_DENV.r: Merge peak and call diff for DENV vs Mock.

2-DRIPseq_MergePeak_calldif_KO.r: Merge peak and call diff for KO_ZV vs KO_Mock.

3-DRIPseq_MergePeak_calldif_WT.r: Merge peak and call diff for WT_ZV vs WT_Mock.

4-DRIPseq_Heatmap_Fig2B_DENV.r: Genome-wide significantly upregulated and downregulated sites detected by DRIP-seq in DENV-infected SNB-19s (FDR < 0.05). DRIP-seq: Heatmap of comparison between DENV vs Mock. 

5-Scatter_Fig3A_mr.r: Scatter plots illustrates DRIP-seq reads in Mock- or ZIKV_mr-infected cells. Genes up/down regulated in both DRIP-seq and Bru-seq are highlighted.

6-Scatter_Fig3A_pr.r: Scatter plots illustrates DRIP-seq reads in Mock- or ZIKV_pr-infected cells. Genes up/down regulated in both DRIP-seq and Bru-seq are highlighted.

7-Scatter_Fig3D.r: Significantly upregulated ISGs or Non-ISGs in Bru-seq. Upregulated genes that are also detected by DRIP-seq are highlighted.

8-BoxPlot_Fig3E.r: DRIP-seq logFC for Bru-seq upregulated ISGs and Non-ISGs. **** indicates P ＜ 0.0001 (Unpaired t test).

9-GO_supplementFig3A_Four_up.r: Gene ontology of common upregulated genes detected by both Drip-seq and Bru-seq in ZIKVMR- or ZIKVPR-infected cells.
10-DRIPseq_Heatmap_supplementFig4D_KO.r: Genome-wide significantly upregulated and downregulated R-loop regions detected by DRIP-seq in ZIKVMR-infected U-251 MG/IFNAR1KO cells (FDR < 0.05)
11-DRIPseq_Heatmap_supplementFig4D_WT.r: Genome-wide significantly upregulated and downregulated R-loop regions detected by DRIP-seq in ZIKVMR-infected U-251 MG/IFNAR1WT cells (FDR < 0.05).
12-Bruseq_Heatmap_Rloop_processing_genes.r: Heatmap demonstrate the normalized Bru-seq reads in 5 R-loop processing genes in ZIKVMR and ZIKVPR-infected cells.
13-DRIPseq_Heatmap_Fig2B_IFN.r: Genome-wide significantly upregulated and downregulated sites detected by DRIP-seq in IFN-infected cells (FDR < 0.05). DRIP-seq: Heatmap of comparison between IFN vs Mock. 
14-GO_Bru_down_genes.r
15-Bruseq_ReadsCount_PRV_Mock.r: Count Bruseq reads in genebody for PRV and Mock.
16-Bruseq_ReadsCount_MR_Mock.r: Count Bruseq reads in gene body for MR and Mock.
17-Bruseq_DifferentialAnalysis.r: Bru-seq reads were counted in gene body and stored in a matrix called "mat.rda". Reads were normalized to total mapped reads and call diff by DESeq2
