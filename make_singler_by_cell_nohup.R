#make 
library(Seurat)
library(qs)
library(RColorBrewer)
library(ComplexHeatmap)
library(dplyr)
library(reshape2)
library(patchwork)
#7. cell annotation by SingleR
library(SingleR)
#refdata<-ref_data
#testdata<-stem_vivo
singler_anno<-function(reffile="/data/wangrj/single_cell/ref_singlecell_data/embryo_named_230513.qs",testdata,method="cluster"){
  refdata=qread(reffile)
  ref_sce<-as.SingleCellExperiment(refdata)
  #trained<-trainSingleR(ref_sce,labels=ref_sce$comparetype)
  test_sce<-as.SingleCellExperiment(testdata)
  if(method=="cluster"){
    annodata<-SingleR(test=test_sce,ref=ref_sce,labels=ref_sce$comparetype,clusters=test_sce$RNA_snn_res.0.8)
    #rename idents
    Idents(testdata)<-testdata$RNA_snn_res.0.8
    Celltype<-annodata$pruned.labels
    names(Celltype)<-levels(testdata)
    testdata<-RenameIdents(testdata,Celltype)
    testdata$Celltype_SingleR<-Idents(testdata)
    return(testdata)
  }
  else if(method=="cell"){
    annodata<-SingleR(test=test_sce,ref=ref_sce,labels=ref_sce$comparetype)
    qsave(annodata,"invivo_singler_anno_by_cell.qs")
  }
}
#refdata<-qread("/data/wangrj/single_cell/ref_singlecell_data/embryo_named_230513.qs")
stem_vivo_no_outier<-qread("/data/wangrj/single_cell/all_matrix/RM_doublets/TC_InVivo/RM_outliers/stem_data_TimeCourse_add_InVivo_no_outlier.qs")
result<-singler_anno(testdata=stem_vivo_no_outier,method="cell")



