
# NCOR1 and SMRT binding
computeMatrix reference-point -S ../../../bigWig/SMRT_0hr_9M.ucsc.bigWig ../../../bigWig/SMRT_CpG_9M.ucsc.bigWig ../../../bigWig/NCoR1_Uns_9M.ucsc.bigWig ../../../bigWig/NCoR1_CpG_9M.ucsc.bigWig -R NCoR1_Uns_unique_peaks.bed NCoR1_CpG_unique_peaks.bed  NCoR1_SMRT_Uns_CpG_unique_peaks.bed  NCoR1_SMRT_CpG_unique_peaks.bed SMRT_Uns_unique_peaks.bed SMRT_CpG_unique_peaks.bed SMRT_Uns_CpG_unique_peaks.bed  -b 2000 -a 2000 -o NCoR1_SMRT_Uns_CpG_unique_matrix.mat.gz -p 20

plotHeatmap --matrixFile NCoR1_SMRT_Uns_CpG_unique_matrix.mat.gz --outFileName NCoR1_SMRT_Uns_CpG_unique_matrix.pdf --plotFileFormat pdf --heatmapHeight 0.5 --missingDataColor 1 --verbose --zMin 0 0 --zMax 1 1 --heatmapWidth 2 --whatToShow 'heatmap and colorbar' --colorList "#DEB132,#058983" "#DEB132,#058983" "#058983,#DEB132" "#058983,#DEB132" --samplesLabel "SMRT_0h" "SMRT_6h_CpG" "NCoR1_0h" "NCoR1_CpG_6h"


# NCoR1_CpG NCoR1_SMRT_CPG and SMRT_CpG
computeMatrix reference-point -S ../../../bigWig/SMRT_0hr_9M.ucsc.bigWig ../../../bigWig/SMRT_CpG_9M.ucsc.bigWig ../../../bigWig/NCoR1_Uns_9M.ucsc.bigWig ../../../bigWig/NCoR1_CpG_9M.ucsc.bigWig -R NCoR1_CpG_unique_peaks.bed NCoR1_SMRT_CpG_unique_peaks.bed SMRT_CpG_unique_peaks.bed -b 2000 -a 2000 -o matrix.mat.gz -p 20


plotHeatmap --matrixFile matrix.mat.gz --outFileName matrix.pdf --plotFileFormat pdf --heatmapHeight 0.5 --missingDataColor 1 --verbose --zMin 0 0 --zMax 1 1 --heatmapWidth 2 --whatToShow 'heatmap and colorbar' --colorList "#DEB132,#058983" "#DEB132,#058983" "#058983,#DEB132" "#058983,#DEB132" --samplesLabel "SMRT_0h" "SMRT_6h_CpG" "NCoR1_0h" "NCoR1_CpG_6h"


# NoR1_CpG NCoR1_SMRT_CPG and SMRT_CpG
computeMatrix reference-point -S ../../../bigWig/SMRT_0hr_9M.ucsc.bigWig ../../../bigWig/SMRT_CpG_9M.ucsc.bigWig ../../../bigWig/NCoR1_Uns_9M.ucsc.bigWig ../../../bigWig/NCoR1_CpG_9M.ucsc.bigWig -R NCoR1_Uns_unique_peaks.bed NCoR1_SMRT_Uns_CpG_unique_peaks.bed  NCoR1_SMRT_CpG_unique_peaks.bed SMRT_Uns_unique_peaks.bed SMRT_Uns_CpG_unique_peaks.bed -b 2000 -a 2000 -o NCoR1_SMRT_Uns_CpG_binding_matrix1.mat.gz -p 20


plotHeatmap --matrixFile NCoR1_SMRT_Uns_CpG_binding_matrix1.mat.gz --outFileName NCoR1_SMRT_Uns_CpG_binding_matrix1.mat.pdf --plotFileFormat pdf --heatmapHeight 0.5 --missingDataColor 1 --verbose --zMin 0 0 --zMax 1 1 --heatmapWidth 2 --whatToShow 'heatmap and colorbar' --colorList "#DEB132,#058983" "#DEB132,#058983" "#058983,#DEB132" "#058983,#DEB132" --samplesLabel "SMRT_0h" "SMRT_6h_CpG" "NCoR1_0h" "NCoR1_CpG_6h" --regionsLabel "6446" "4271" "1678" "1063" "1067" --refPointLabel "0" --heatmapWidth 2







# H3K27ac acetylation 
computeMatrix reference-point -S ../../../../NCoR1_H3K27ac_analysis/bedGraph/Emp_0hr.bw ../../../../NCoR1_H3K27ac_analysis/bedGraph/NCoR1_0hr.bw ../../../../NCoR1_H3K27ac_analysis/bedGraph/Emp_2hr_CpG.bw ../../../../NCoR1_H3K27ac_analysis/bedGraph/NCoR1_2hr_CpG.bw ../../../../NCoR1_H3K27ac_analysis/bedGraph/Emp_6hr_CpG.bw ../../../../NCoR1_H3K27ac_analysis/bedGraph/NCoR1_6hr_CpG.bw -R NCoR1_Uns_unique_peaks.bed NCoR1_CpG_unique_peaks.bed  NCoR1_Uns_CpG_unique_peaks.bed NCoR1_SMRT_CpG_unique_peaks.bed SMRT_Uns_CpG_unique_peaks.bed  SMRT_CpG_unique_peaks.bed  SMRT_Uns_unique_peaks.bed -b 2000 -a 2000 --skipZeros -o matrix.mat.gz -p 20

plotHeatmap --matrixFile matrix.mat.gz --outFileName matrix_NCoR1_SMRT_0h_6h_CpG_H3k27ac.pdf --plotFileFormat pdf --heatmapHeight 0.5 --missingDataColor 1 --verbose --zMin 0 0 --zMax 1 1 --heatmapWidth 2 --whatToShow 'heatmap and colorbar' --colorList "red,blue" 
