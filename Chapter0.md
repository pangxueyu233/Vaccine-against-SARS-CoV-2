# 1. code of data pre-analysis

**Preprocessing and Alignment of scRNA-seq and BCR Data Using Cell Ranger**

**Description:**
 This pipeline performs preprocessing and alignment of single-cell transcriptomic (scRNA-seq), and B-cell receptor (BCR) sequencing data from Irina’s mouse study, using the Cell Ranger 7.2.0 suite. The script includes steps for three experimental groups — **Ad5_prowt**, **PBS**, and **STing** — covering both gene expression and immune repertoire data. It runs `cellranger count` to generate gene expression matrices with intronic reads included, and uses `cellranger vdj` to reconstruct immune repertoires based on V(D)J sequencing. All relevant outputs, including filtered gene expression matrices and clonotype annotations, are copied and organized into designated result directories for downstream integration and analysis.

~~~R
mkdir -p ./xiangyu/workshop/HWQ_scOmics/res/Ad5_prowt
cd ./xiangyu/workshop/HWQ_scOmics/res/Ad5_prowt
./xiangyu/programme/cellranger-7.2.0/bin/cellranger count \
--id=Ad5_prowt --localcores 20 --sample=Ad5_prowt --nosecondary --include-introns true \
--transcriptome=/local./xiangyu/programme/genome_index/refdata-gex-mm10-2020-A \
--fastqs=/workdir/Omics/HWQ_scOmics/'LC-X20240825003-10X免疫组学测序(10XVDJ)'/'2024-09-12_17:33:26'/Data/CleanData/mRNA/Ad5_prowt \
mkdir -p ./xiangyu/workshop/HWQ_scOmics/res/PBS
cd ./xiangyu/workshop/HWQ_scOmics/res/PBS
./xiangyu/programme/cellranger-7.2.0/bin/cellranger count \
--id=PBS --localcores 20 --sample=PBS --nosecondary --include-introns true \
--transcriptome=/local./xiangyu/programme/genome_index/refdata-gex-mm10-2020-A \
--fastqs=/workdir/Omics/HWQ_scOmics/'LC-X20240825003-10X免疫组学测序(10XVDJ)'/'2024-09-12_17:33:26'/Data/CleanData/mRNA/PBS
mkdir -p ./xiangyu/workshop/HWQ_scOmics/res/STing
cd ./xiangyu/workshop/HWQ_scOmics/res/STing
./xiangyu/programme/cellranger-7.2.0/bin/cellranger count \
--id=STing --localcores 20 --sample=STing --nosecondary --include-introns true \
--transcriptome=/local./xiangyu/programme/genome_index/refdata-gex-mm10-2020-A \
--fastqs=/workdir/Omics/HWQ_scOmics/'LC-X20240825003-10X免疫组学测序(10XVDJ)'/'2024-09-12_17:33:26'/Data/CleanData/mRNA/STing

mkdir ./xiangyu/workshop/HWQ_scOmics/res/Ad5_prowt_TCR
cd ./xiangyu/workshop/HWQ_scOmics/res/Ad5_prowt_TCR
./xiangyu/programme/cellranger-7.2.0/bin/cellranger vdj \
--id=Ad5_prowt --sample=Ad5_prowt --localcores=20 \
--reference=/local./xiangyu/programme/genome_index/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0/ \
--fastqs=/workdir/Omics/HWQ_scOmics/'LC-X20240825003-10X免疫组学测序(10XVDJ)'/'2024-09-12_17:33:26'/Data/CleanData/TCR/Ad5_prowt
mkdir ./xiangyu/workshop/HWQ_scOmics/res/PBS_TCR
cd ./xiangyu/workshop/HWQ_scOmics/res/PBS_TCR
./xiangyu/programme/cellranger-7.2.0/bin/cellranger vdj \
--id=PBS --sample=PBS --localcores=20 \
--reference=/local./xiangyu/programme/genome_index/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0/ \
--fastqs=/workdir/Omics/HWQ_scOmics/'LC-X20240825003-10X免疫组学测序(10XVDJ)'/'2024-09-12_17:33:26'/Data/CleanData/TCR/PBS
mkdir ./xiangyu/workshop/HWQ_scOmics/res/STing_TCR
cd ./xiangyu/workshop/HWQ_scOmics/res/STing_TCR
./xiangyu/programme/cellranger-7.2.0/bin/cellranger vdj \
--id=STing --sample=STing --localcores=20 \
--reference=/local./xiangyu/programme/genome_index/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0/ \
--fastqs=/workdir/Omics/HWQ_scOmics/'LC-X20240825003-10X免疫组学测序(10XVDJ)'/'2024-09-12_17:33:26'/Data/CleanData/TCR/STing

mkdir ./xiangyu/workshop/HWQ_scOmics/res/Ad5_prowt_BCR
cd ./xiangyu/workshop/HWQ_scOmics/res/Ad5_prowt_BCR
./xiangyu/programme/cellranger-7.2.0/bin/cellranger vdj \
--id=Ad5_prowt --sample=Ad5_prowt --localcores=20 \
--reference=/local./xiangyu/programme/genome_index/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0/ \
--fastqs=/workdir/Omics/HWQ_scOmics/'LC-X20240825003-10X免疫组学测序(10XVDJ)'/'2024-09-12_17:33:26'/Data/CleanData/BCR/Ad5_prowt
mkdir ./xiangyu/workshop/HWQ_scOmics/res/PBS_BCR
cd ./xiangyu/workshop/HWQ_scOmics/res/PBS_BCR
./xiangyu/programme/cellranger-7.2.0/bin/cellranger vdj \
--id=PBS --sample=PBS --localcores=20 \
--reference=/local./xiangyu/programme/genome_index/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0/ \
--fastqs=/workdir/Omics/HWQ_scOmics/'LC-X20240825003-10X免疫组学测序(10XVDJ)'/'2024-09-12_17:33:26'/Data/CleanData/BCR/PBS
mkdir ./xiangyu/workshop/HWQ_scOmics/res/STing_BCR
cd ./xiangyu/workshop/HWQ_scOmics/res/STing_BCR
./xiangyu/programme/cellranger-7.2.0/bin/cellranger vdj \
--id=STing --sample=STing --localcores=20 \
--reference=/local./xiangyu/programme/genome_index/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0/ \
--fastqs=/workdir/Omics/HWQ_scOmics/'LC-X20240825003-10X免疫组学测序(10XVDJ)'/'2024-09-12_17:33:26'/Data/CleanData/BCR/STing

mkdir -p ./xiangyu/workshop/HWQ_scOmics/res/res/Ad5_prowt/
mkdir -p ./xiangyu/workshop/HWQ_scOmics/res/res/Ad5_prowt_BCR/
mkdir -p ./xiangyu/workshop/HWQ_scOmics/res/res/Ad5_prowt_TCR/
mkdir -p ./xiangyu/workshop/HWQ_scOmics/res/res/PBS/
mkdir -p ./xiangyu/workshop/HWQ_scOmics/res/res/PBS_BCR/
mkdir -p ./xiangyu/workshop/HWQ_scOmics/res/res/PBS_TCR/
mkdir -p ./xiangyu/workshop/HWQ_scOmics/res/res/STing/
mkdir -p ./xiangyu/workshop/HWQ_scOmics/res/res/STing_BCR/
mkdir -p ./xiangyu/workshop/HWQ_scOmics/res/res/STing_TCR/

cp -r ./xiangyu/workshop/HWQ_scOmics/res/PBS/PBS/outs/filtered_feature_bc_matrix ./xiangyu/workshop/HWQ_scOmics/res/res/PBS/
cp -r ./xiangyu/workshop/HWQ_scOmics/res/STing/STing/outs/filtered_feature_bc_matrix ./xiangyu/workshop/HWQ_scOmics/res/res/STing/
cp -r ./xiangyu/workshop/HWQ_scOmics/res/Ad5_prowt/Ad5_prowt/outs/filtered_feature_bc_matrix ./xiangyu/workshop/HWQ_scOmics/res/res/Ad5_prowt/
cp -r ./xiangyu/workshop/HWQ_scOmics/res/PBS_TCR/PBS/outs/filtered_contig_annotations.csv ./xiangyu/workshop/HWQ_scOmics/res/res/PBS_TCR/
cp -r ./xiangyu/workshop/HWQ_scOmics/res/STing_TCR/STing/outs/filtered_contig_annotations.csv ./xiangyu/workshop/HWQ_scOmics/res/res/STing_TCR/
cp -r ./xiangyu/workshop/HWQ_scOmics/res/Ad5_prowt_TCR/Ad5_prowt/outs/filtered_contig_annotations.csv ./xiangyu/workshop/HWQ_scOmics/res/res/Ad5_prowt_TCR/
cp -r ./xiangyu/workshop/HWQ_scOmics/res/PBS_BCR/PBS/outs/filtered_contig_annotations.csv ./xiangyu/workshop/HWQ_scOmics/res/res/PBS_BCR/
cp -r ./xiangyu/workshop/HWQ_scOmics/res/STing_BCR/STing/outs/filtered_contig_annotations.csv ./xiangyu/workshop/HWQ_scOmics/res/res/STing_BCR/
cp -r ./xiangyu/workshop/HWQ_scOmics/res/Ad5_prowt_BCR/Ad5_prowt/outs/filtered_contig_annotations.csv ./xiangyu/workshop/HWQ_scOmics/res/res/Ad5_prowt_BCR/

cp -r ./xiangyu/workshop/HWQ_scOmics/res/PBS_TCR/PBS/outs/clonotypes.csv ./xiangyu/workshop/HWQ_scOmics/res/res/PBS_TCR/
cp -r ./xiangyu/workshop/HWQ_scOmics/res/STing_TCR/STing/outs/clonotypes.csv ./xiangyu/workshop/HWQ_scOmics/res/res/STing_TCR/
cp -r ./xiangyu/workshop/HWQ_scOmics/res/Ad5_prowt_TCR/Ad5_prowt/outs/clonotypes.csv ./xiangyu/workshop/HWQ_scOmics/res/res/Ad5_prowt_TCR/
cp -r ./xiangyu/workshop/HWQ_scOmics/res/PBS_BCR/PBS/outs/clonotypes.csv ./xiangyu/workshop/HWQ_scOmics/res/res/PBS_BCR/
cp -r ./xiangyu/workshop/HWQ_scOmics/res/STing_BCR/STing/outs/clonotypes.csv ./xiangyu/workshop/HWQ_scOmics/res/res/STing_BCR/
cp -r ./xiangyu/workshop/HWQ_scOmics/res/Ad5_prowt_BCR/Ad5_prowt/outs/clonotypes.csv ./xiangyu/workshop/HWQ_scOmics/res/res/Ad5_prowt_BCR/
~~~

