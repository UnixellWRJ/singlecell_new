#main

#sub9
sub9<-qread("recluster_of_cluster9.qs")
pdf("circplot_of_sub9.pdf",height=7,width=9)
pall<-make_circled_featureplot(sub9,features=features,circ_cluster=1)
dev.off()

#in vivo
setwd("/data/wangrj/single_cell/all_matrix/RM_doublets/TC_InVivo/")
stem_vivo<-qread("stem_data_integrated_TimeCourse_add_InVivo_rna_no_vars.qs")
marker_file<-read.csv("stem_data_all_markers_TimeCourse_add_InVivo_by_rna.csv",header=T,row.names=1)
#make_celltype_enrichment_analysis(marker_file,fname="InVivo")
make_celltype_enrichment_analysis(marker_file,fname="InVivo_top50_gene",genenum=50)

make_gsea_celltype_analysis(marker_file,fname="InVivo_GSEA_analysis",width_plot=22)
#cell type annotation
CellType<-c("Stromal","hGaba","Stromal","hNbM","Olf_NE_secretory","hNbM","Stromal","Stromal","hRgl2b","hDA2","hSert","Stromal","Ensheathing_glia","hNbML1","hSert","hRgl3","hRgl1","hGaba","hNbML1","Olf_NE_ciliated","Astrocytes","hSert","None","Neuroendocrine cells","hRgl1/3","OPC")
names(CellType)<-levels(stem_vivo$RNA_snn_res.0.8)
Idents(stem_vivo)<-stem_vivo$RNA_snn_res.0.8
stem_vivo<-RenameIdents(stem_vivo,CellType)
stem_vivo$CellType<-Idents(stem_vivo)
qsave(stem_vivo,"stem_data_integrated_TimeCourse_add_InVivo_rna_no_vars.qs")
dimplot_making(stem_vivo,cluster="CellType",fname="stem_vivo_by_CellType",colorset="Set1",colornum=20,height_m=6,width_m=7,height_s=9)

FeaturePlot(stem_vivo,features="MBP",split.by="Sample")

#in vivo without outlier
setwd("/data/wangrj/single_cell/all_matrix/RM_doublets/TC_InVivo/RM_outliers/")
stem_vivo_no_outlier<-qread("/data/wangrj/single_cell/all_matrix/RM_doublets/TC_InVivo/RM_outliers/stem_data_TimeCourse_add_InVivo_no_outlier.qs")
Idents(stem_vivo)<-stem_vivo$RNA_snn_res.0.8
plot1<-DimPlot(stem_vivo,label=T)
plot2<-DimPlot(stem_vivo_no_outier,label=T)
wrap_plots(plot1,plot2)
marker_file<-read.csv("/data/wangrj/single_cell/all_matrix/RM_doublets/TC_InVivo/RM_outliers/stem_data_all_markers_TimeCourse_add_InVivo_by_rna.csv",header=T,row.names=1)

dimplot_making(stem_vivo_no_outlier,fname="In_Vivo_no_outlier",colorset="Accent",colornum=30)
percent_graph_make(stem_vivo_no_outlier,fname="TC_InVivo_no_outlier_freey",scale="free_y")
percent_graph_make(stem_vivo_no_outlier,fname="TC_InVivo_no_outlier")
make_dot_heatmap(stem_vivo_no_outlier,marker_file=marker_file,fname="TC_InVivo_no_outlier",ntop=7,color_h="#FF3399")
stem_vivo_no_outlier_singler<-singler_anno(testdata=stem_vivo_no_outier)
DimPlot(stem_vivo_no_outlier_singler,group.by="Celltype_SingleR",label=T)
gsea_result<-make_gsea_celltype_analysis(marker_file,fname="InVivo_GSEA_analysis_without_outlier",width_plot=22)
write.csv(gsea_result,"gsea_result_rm_outliers_invivo.csv")
#percentage graph by CellType
percent_graph_make(stem_vivo_no_outlier,fname="TC_InVivo_no_outlier_freey_CellType_by_GSEA",clustername="CellType_by_GSEA",
                   scale="free_y",txt_size=11)
#percent_graph_make(stem_vivo_no_outlier,fname="TC_InVivo_no_outlier")




#in vitro
setwd("/data/wangrj/single_cell/all_matrix/RM_doublets/TC_InVitro/RM_outliers_in_vitro/")
marker_file<-read.csv("stem_data_all_markers_TimeCourse_add_TC_InVitro_by_rna.csv",header=T,row.names=1)
stem_vitro<-qread("stem_data_integrated_TimeCourse_add_TC_InVitro_rna.qs")
#make_celltype_enrichment_analysis(marker_file,fname="InVitro")
make_celltype_enrichment_analysis(marker_file,fname="InVitro_top50_gene",genenum=50)
make_gsea_celltype_analysis(marker_file,fname="InVitro_GSEA_analysis",width_plot=22)
#ES_IPS













