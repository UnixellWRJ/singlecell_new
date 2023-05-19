source("/data/wangrj/single_cell/all_matrix/RM_doublets/basic_singlecell_pipeline.R")
#example: in vivo analysis data
#These lines are changed according to sample
#directory to save the results
setwd("/data/wangrj/single_cell/all_matrix/RM_doublets/TC_InVivo/RM_outliers/")
#Seurat object load
stem_data<-qread("/data/wangrj/single_cell/all_matrix/RM_doublets/TC_InVivo/RM_outliers/stem_data_TimeCourse_add_InVivo_no_outlier.qs")
#marker files calcluated by FindAllMarkers
marker_file<-read.csv("/data/wangrj/single_cell/all_matrix/RM_doublets/TC_InVivo/RM_outliers/stem_data_all_markers_TimeCourse_add_InVivo_by_rna.csv",header=T,row.names=1)
#prefix used for output results
fname<-"InVivo_rm_outlier"
#cluster name used 
cluster_name="RNA_snn_res.0.8"

#dimplot making, merged and split by sample, output pdf files
dimplot_making(stem_data,fname=fname,cluster=cluster_name)

#calculate the percentage of each cluster and make lolliplog
#freey
percent_graph_make(stem_data,fname=fname,scale="free_y")
#fixed
percent_graph_make(stem_data,fname=fname)

#make heatmap in each cluster for top n genes sorted by avg_log2FC, default n is 5
make_dot_heatmap(stem_data,marker_file=marker_file,fname=fname,ntop=7,color_h="#FF3399")

#cell type annotation by GSEA with self-made reference.
gsea_result_selfmade<-make_gsea_celltype_analysis(marker_file,fname=fname,gmt_file="/data/wangrj/single_cell/gmt_files/self_made_celltype_anno_genelist.gmt",width_plot=22)

#filter for plasma membrane genes
marker_membrane<-is_plasma_membrane(marker_file)
write.csv(paste0(fname,"_plasma_membrane_diff_markers.csv"),quote=F)

#filter for transcription factor genes
marker_tf<-is_transcription_factor<-is_transcription_factor(marker_file)
write.csv(paste0(fname,"_transcription_factor_diff_markers.csv"),quote=F)


