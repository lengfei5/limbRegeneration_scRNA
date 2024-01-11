##########################################################################
##########################################################################
# Project: axolotl WE (wound epidermis) trajectory 
# Script purpose: process and analyze trajectory of Leigh et al., 2018
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri May 19 11:01:01 2023
##########################################################################
##########################################################################
rm(list = ls())

library(Seurat)
library(decoupleR)
library(tictoc)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library("viridis")

version.analysis = '_axolotl_Leigh2018_20230519'
resDir = paste0("../results/scRNAseq", version.analysis)
RdataDir = paste0(resDir, '/Rdata')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts'
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_scRNAseq.R')
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_Visium.R')

dataDir = '/groups/tanaka/People/current/jiwang/projects/limbRegeneration_scRNA/raw_NGS/axolotl/Leigh_2018'
metaDir = paste0(dataDir, '/indrops/Leigh_et_al_2018_Supplementary_R_code')

########################################################
########################################################
# Section I: import processed Seurat object from Steven Blair (Whited lab)
# 
########################################################
########################################################
load(paste0(dataDir, '/processed_data/nLeigh_seuratConverted_natCom2018.rdata'))

## mature sample cell types
meta = read.table(file = paste0(metaDir, '/homeostasis.meta.data.txt'), header = TRUE)
aa = seu_intact
aa$celltypes = NA
aa$celltypes = meta$Cell_type[match(colnames(aa), meta$cell)]

DimPlot(aa, group.by = 'celltypes', label = TRUE, repel = FALSE)

seu_intact = aa

rm(aa)

# wound healing
meta = read.table(file = paste0(metaDir, '/wound_healing.meta.data.txt'), header = TRUE)
aa = seu_3dpa
aa$celltypes = NA
aa$celltypes = meta$cell_type[match(colnames(aa), meta$Cell)]

DimPlot(aa, group.by = 'celltypes', label = TRUE, repel = FALSE)
seu_3dpa = aa

rm(aa)

# 14dpa
meta = read.table(file = paste0(metaDir, '/early_bud_blastema.meta.data.txt'), header = TRUE)
aa = seu_14dpa
aa$celltypes = NA
aa$celltypes = meta$cell_type[match(colnames(aa), meta$cell)]

DimPlot(aa, group.by = 'celltypes', label = TRUE, repel = FALSE)
seu_14dpa = aa

rm(aa)

# 23dpa
meta = read.table(file = paste0(metaDir, '/medium_bud_blastema.meta.data.txt'), header = TRUE)
aa = seu_23dpa
aa$celltypes = NA
aa$celltypes = meta$Cell_type[match(colnames(aa), meta$cell)]

DimPlot(aa, group.by = 'celltypes', label = TRUE, repel = FALSE)
seu_23dpa = aa

rm(aa)
rm(meta)

saveRDS(seu_intact, file = paste0(RdataDir, '/seuratObj_processedWhited_0dpa.rds'))
saveRDS(seu_3dpa, file = paste0(RdataDir, '/seuratObj_processedWhited_3dpa.rds'))
saveRDS(seu_14dpa, file = paste0(RdataDir, '/seuratObj_processedWhited_14dpa.rds'))
saveRDS(seu_23dpa, file = paste0(RdataDir, '/seuratObj_processedWhited_23dpa.rds'))
