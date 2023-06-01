library(Seurat)
library(qs)
library(RColorBrewer)
library(ComplexHeatmap)
library(dplyr)
library(reshape2)
library(patchwork)

#colorlist
#18 colors for refdata
colorset_18<-c("grey80","grey80","grey80","grey80","grey80","#FF00CC","#FF0033","#66CCFF","#66FF00","#66CCCC","grey80","grey80",
                       "grey80","grey80","grey80","grey80","grey80","grey80")

#celltype annotation features
##vivo
features_fib<-c("COL1A1","COL1A2","VIM","FN1","PDGFRB","PDGFRA","LAMA1","ANPEP")
features_peri<-c("PDGFRB", "CSPG4", "DES", "KCNJ8", "ABCC9", "ANPEP")
features_vlmc<-c("PDGFRB", "ACTA2", "ANPEP", "CSPG4", "MCAM", "DES")
features_dural_fib<-c("FXYD5", "FOXP1", "SIX1")
features_ara_fib<-c("CRABP2", "ALDH1A2", "SLC6A13")
features_pial_fib<-c("S100A6","NGFR")
features_DA<-"TH"
features_DA_A9<-c("ALDH1A1","KCNJ6","NR4A2")
features_DA_A10<-"CALB1"
features_Glut<-c("SLC17A6","PITX2")
#SLC32A1 is a specific marker
features_gaba<-c("SLC32A1","GAD1","GAD2")
features_Sero<-c('SLC6A4','SLC17A8','FEV','TPH2','GATA3')
features_OTP<-"OTP"
#AQP4 is the specific marker
features_Astrocyte<-c('AQP4','SLC1A3','GFAP','GJA1')
features_astro_A1<-"S100A10"
features_astro_A2<-"C3"
features_od<-c('PLP1','MOG','MBP')
features_opc<-c('PDGFRA','SOX10','CSPG4','OLIG2','OLIG1')
features_peric<-"PDGFRB"
features_Nb<-c('NEUROD1','NEUROD2','ASCL1','INSM1')
features_neuron<-"STMN2"
features_prog<-"SOX2"
##vitro
features_position<-c("OTX2","EN1","PLP1","LMX1A","FOXA2")
features_all<-c(features_fib,features_peri,features_vlmc,features_dural_fib,features_ara_fib,features_pial_fib,features_DA,features_DA_A9,features_DA_A10,features_Glut,features_gaba,features_Sero,features_OTP,features_Astrocyte,features_astro_A1,features_astro_A2,features_od,features_opc,features_peric,features_Nb,features_neuron,features_prog,features_position)


#1. dimplot making, include merged and split
dimplot_making<-function(sdata,fname="stem.data",cluster="RNA_snn_res.0.8",colorset="Pastel1",colornum=50,height_m=5,width_m=6,height_s=8,
                         Legend="No"){
  merge_name<-paste0("dimplot_",fname,"_merged.pdf")
  split_name<-paste0("dimplot_",fname,"_split.pdf")
  Idents(sdata)<-sdata[[cluster]]
  if(Legend=="No"){
    pdf(merge_name,height=height_m,width=width_m)
    print(DimPlot(sdata,label=T)+NoLegend()+scale_color_manual(values=colorRampPalette(brewer.pal(8,colorset))(colornum)))
    dev.off()
    sample_num<-length(unique(sdata$Sample))
    colnum<-ceiling(sample_num/2)
    width_s=colnum/2*height_s+2
    pdf(split_name,height=height_s,width=width_s)
    print(DimPlot(sdata,split.by="Sample",label=T,ncol=colnum)+NoLegend()+
            scale_color_manual(values=colorRampPalette(brewer.pal(8,colorset))(colornum)))
    dev.off()
  }
  else if(Legend=="Yes"){
    width_m1=width_m+2
    pdf(merge_name,height=height_m,width=width_m1)
    print(DimPlot(sdata,label=T)+NoLegend()+scale_color_manual(values=colorRampPalette(brewer.pal(8,colorset))(colornum)))
    dev.off()
    sample_num<-length(unique(sdata$Sample))
    colnum<-ceiling(sample_num/2)
    width_s1=colnum/2*height_s+4
    pdf(split_name,height=height_s,width=width_s1)
    print(DimPlot(sdata,split.by="Sample",label=T,ncol=colnum)+NoLegend()+
            scale_color_manual(values=colorRampPalette(brewer.pal(8,colorset))(colornum)))
    dev.off()
  }
}

#2. percentage calculation and graph making
#three types of graph to choose: lollipop, bar, point-line
cal_percentage<-function(samplename,clustername,metadata){
  samlist<-unique(metadata[[samplename]])
  clist<-unique(metadata[[clustername]])
  metadata$Sample<-metadata[[samplename]]
  metadata$Cluster<-metadata[[clustername]]
  perdata<-c()
  for(sa in samlist){
    sub1<-subset(metadata,Sample==sa)
    num1<-length(sub1$Sample)
    clsub<-unique(sub1$Cluster)
    for(cl in clist){
      if(cl %in% clsub){
        sub2<-subset(sub1,Cluster==cl)
        num2<-length(sub2$Cluster)
        per<-num2/num1
      }
      else{
        per<-0
      }
      tmpdata<-data.frame(sa,cl,per)
      perdata<-rbind(perdata,tmpdata)
    }
  }
  colnames(perdata)<-c("Sample","Cluster","Percentage")
  perdata$Sample<-factor(perdata$Sample,levels=unique(metadata$Sample))
  return(perdata)
}
percent_graph_make<-function(sdata,fname="stem.data",clustername="RNA_snn_res.0.8",graph_type="lollipop",scale="fixed",
                             txt_size=15){
  metadata<-sdata@meta.data
  color_num<-length(unique(metadata$Sample))
  perdata<-cal_percentage("Sample",clustername,metadata)
  perdata$Sample<-factor(perdata$Sample,unique(metadata$Sample))
  if(clustername=="RNA_snn_res.0.8"){
    level_max<-max(as.numeric(as.character(perdata$Cluster)))
    perdata$Cluster<-factor(perdata$Cluster,levels=0:level_max)
  }
  clusternum<-length(unique(perdata$Cluster))
  #num_w<-ceiling(sqrt(clusternum))+1
  num_w<-ceiling(sqrt(clusternum))
  num_h<-ceiling((clusternum+1)/num_w)
  pdfname=paste0("percentage_of_samples_",fname,"_",graph_type,"_",scale,"_plot.pdf")
  height_p<-num_h*1.6
  width_p<-num_w*1.6+2
  if(graph_type=="point-line"){
    p1<-ggplot(perdata,aes(x=Sample,y=Percentage,group=Cluster,color=Sample))+
      geom_line(linewidth=1)+geom_point(size=2,aes(shape=Stage))+
      scale_color_manual(values=colorRampPalette(brewer.pal(8,"Set2"))(color_num))+
      theme_classic()+theme(text=element_text(size=txt_size),axis.text.x = element_text(angle=60,hjust=1))+
      facet_wrap(~Cluster,ncol=num_w,scale=scale)
  }
  else if(graph_type=="bar"){
    p1<-ggplot(perdata,aes(x=Sample,y=Percentage,fill=Sample))+geom_bar(stat="identity",position="dodge",color="black")+
      scale_fill_manual(values=colorRampPalette(brewer.pal(8,"Set2"))(color_num))+
      theme_classic()+theme(text=element_text(size=txt_size),axis.text.x = element_text(angle=90,hjust=1))+
      facet_wrap(~Cluster,ncol=num_w,scale=scale)
  }
  else{
    p1<-ggplot(perdata,aes(x=Sample,y=Percentage,group=Cluster,color=Sample))+
      geom_point(size=3)+
      geom_segment(aes(x=Sample,xend=Sample,y=0,yend=Percentage),linewidth=1)+
      scale_color_manual(values=colorRampPalette(brewer.pal(8,"Set2"))(color_num))+
      theme_classic()+theme(text=element_text(size=txt_size),axis.text.x = element_text(angle=60,hjust=1))+
      facet_wrap(~Cluster,ncol=num_w,scale=scale)
  }
  pdf(pdfname,height=height_p,width=width_p)
  print(p1)
  dev.off()
}

#3. heatmap making
#3-1 heatmap making by marker_file
make_dot_heatmap<-function(sdata,marker_file,fname="stem.data",ntop=5,mcenter=T,mscale=T,
                           color_l="grey20",color_m="white",color_h="red"){
  if(ntop=="all"){
    topn_markers<-marker_file
  }
  else{
    marker_file %>% group_by(cluster) %>% top_n(n=ntop,wt=avg_log2FC) -> topn_markers
  }
  genelist<-unique(topn_markers$gene)
  avg_list<-AverageExpression(sdata,assays="RNA",features=genelist)
  avg<-avg_list[["RNA"]]
  avg_scale<-data.frame(t(scale(t(avg),center=mcenter,scale=mscale)))
  avg_scale$Gene<-rownames(avg_scale)
  #scale the data by all
  #Colname<-colnames(avg)
  #Rowname<-rownames(avg)
  #avg_melt<-data.frame(melt(avg,value.name="Exp"))
  #Exp_scaled<-scale(avg_melt$Exp,center=mcenter,scale=mscale)
  #avg_melt$Exp_scaled<-Exp_scaled
  #avg_scale<-dcast(avg_melt,Var1~Var2)
  plot_data<-melt(avg_scale,id.vars="Gene",value.name="Avg_exp",variable.name="Cluster")
  plot_data$Cluster<-gsub("X","",plot_data$Cluster)
  #make pct list
  Pct.exp<-c()
  clist<-unique(plot_data$Cluster)
  for(cl in clist){
    markersub<-subset(marker_file,cluster==cl)
    gesub<-unique(markersub$gene)
    for(gl in genelist){
      if(gl %in% gesub){
        pct<-markersub[gl,"pct.1"]
      }
      else{
        pct<-0
      }
      Pct.exp<-c(Pct.exp,pct)
    }
  }
  plot_data$Pct_exp<-Pct.exp
  plot_data$Cluster<-factor(plot_data$Cluster,levels=unique(plot_data$Cluster))
  plot_data$Gene<-factor(plot_data$Gene,levels=unique(plot_data$Gene))
  pdfname=paste0("heatmap_of_top",ntop,"_gene_",fname,".pdf")
  height_h<-5+(length(unique(plot_data$Cluster))-18)/8
  width_h<-ceiling(length(genelist)/6)
  p1<-ggplot(data = plot_data, mapping = aes_string(x = "Gene", y = "Cluster")) +
    geom_tile(fill="white",color="grey50")+
    geom_point(mapping = aes_string(size = "Pct_exp", color = "Avg_exp")) +
    scale_size_area(max_size = 4.8)+
    scale_color_gradient2(low=color_l,mid=color_m,high=color_h)+
    theme_classic()+
    theme(text=element_text(size=12),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  pdf(pdfname,height=height_h,width=width_h)
  print(p1)
  dev.off()
}
#3-2 heatmap making of selected genes
make_dot_heatmap_features<-function(sdata,features,fname="stem.data",mcenter=T,mscale=T,
                                    color_l="grey20",color_m="white",color_h="red"){
  avg_list<-AverageExpression(sdata,assays="RNA",features=features)
  avg<-avg_list[["RNA"]]
  avg_scale<-data.frame(t(scale(t(avg),center=mcenter,scale=mscale)))
  avg_scale$Gene<-rownames(avg_scale)
  #Colname<-colnames(avg)
  #Rowname<-rownames(avg)
  #avg_melt<-data.frame(melt(avg,value.name="Exp"))
  #Exp_scaled<-scale(avg_melt$Exp,center=mcenter,scale=mscale)
  #avg_melt$Exp_scaled<-Exp_scaled
  #avg_scale<-dcast(avg_melt,Var1~Var2)
  #colnames(avg_scale)[1]<-"Gene"
  plot_data<-melt(avg_scale,id.vars="Gene",value.name="Avg_exp",variable.name="Cluster")
  plot_data$Cluster<-gsub("X","",plot_data$Cluster)
  #make pct data
  #calculate fold change of selected genes
  features_sub<-unique(plot_data$Gene)
  cluster_list<-unique(plot_data$Cluster)
  fcdata<-c()
  for(cl in cluster_list){
    fctmp<-FoldChange(sdata,features=features_sub,ident.1=cl)
    fctmp$Cluster<-cl
    fctmp$Gene<-features_sub
    fcdata<-rbind(fcdata,fctmp)
  }
  plot_data$Pct_exp<-fcdata$pct.1
  plot_data$Cluster<-factor(plot_data$Cluster,levels=unique(plot_data$Cluster))
  plot_data$Gene<-factor(plot_data$Gene,levels=unique(plot_data$Gene))
  pdfname=paste0("heatmap_of_celltype_features_",fname,".pdf")
  height_h<-5+(length(unique(plot_data$Cluster))-18)/8
  width_h<-5+ceiling(length(features_sub)/7.5)
  p1<-ggplot(data = plot_data, mapping = aes_string(x = "Gene", y = "Cluster")) +
    geom_tile(fill="white",color="grey50")+
    geom_point(mapping = aes_string(size = "Pct_exp", color = "Avg_exp")) +
    scale_size_area(max_size = 4.8)+
    scale_color_gradient2(low=color_l,mid=color_m,high=color_h,midpoint = mean(plot_data$Avg_exp))+
    theme_classic()+
    theme(text=element_text(size=12),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  pdf(pdfname,height=height_h,width=width_h)
  print(p1)
  dev.off()
}

#4. make featureplot with circle in specific genes
library(ggplot2)
library(ggraph)
library(ggalt)
#features<-c("TH","CALB2","KCNJ6")
make_circled_featureplot<-function(sdata,features,circ_cluster,cluster_name="RNA_snn_res.0.8",circ_scale=0.9,nmads=5){
  exp_matrix<-FetchData(sdata,vars=features)
  Cluster<-sdata@meta.data[[cluster_name]]
  cell_embeddings<-sdata@reductions$umap@cell.embeddings
  circ_data<-data.frame(Cluster,cell_embeddings)
  circ_select_raw<-subset(circ_data,Cluster==circ_cluster)
  #remove outliers first
  thres_umap1<-c(median(circ_select_raw$UMAP_1)-nmads*mad(circ_select_raw$UMAP_1),median(circ_select_raw$UMAP_1)+nmads*mad(circ_select_raw$UMAP_1))
  thres_umap2<-c(median(circ_select_raw$UMAP_2)-nmads*mad(circ_select_raw$UMAP_2),median(circ_select_raw$UMAP_2)+nmads*mad(circ_select_raw$UMAP_2))
  circ_select1<-subset(circ_select_raw,UMAP_1>thres_umap1[1] & UMAP_1<thres_umap1[2] & 
                         UMAP_2>thres_umap2[1] & UMAP_2<thres_umap2[2] )
  circ_select<-subset(circ_select1,UMAP_1<(circ_scale*max(circ_select_raw$UMAP_1)) & UMAP_1>(circ_scale*min(circ_select_raw$UMAP_1)) & UMAP_2<(circ_scale*max(circ_select_raw$UMAP_2)) &
    UMAP_2>(circ_scale*min(circ_select_raw$UMAP_2)))
  plots<-c()
  for(i in 1:length(exp_matrix)){
    genename<-colnames(exp_matrix)[i]
    Exp<-exp_matrix[,i]
    plot_data<-data.frame(cell_embeddings,Exp)
    p1<-ggplot(plot_data, aes(UMAP_1, UMAP_2, color=Exp)) +
      geom_point(size=1.5) + scale_color_gradient2(low="grey20",mid="grey80",high="#FF0099")+
      geom_encircle(aes(x=UMAP_1,y=UMAP_2),data=circ_select,color="#0000CC",size=2,expand=0.05)+
      theme_light(base_size = 15)+labs(title = genename)+
      theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
      theme(plot.title = element_text(hjust = 0.5))
    plots<-c(plots,list(p1))
  }
  pall<-wrap_plots(plots,ncol=ceiling(sqrt(length(features))))
  return(pall)
}

#4-1 make ellipse circle plot, dimplot
make_ellipse_dimplot<-function(sdata,group.by="RNA_snn_res.0.8",alpha_value=1/10,linesize=1,circsub="all"){
  Cluster<-sdata@meta.data[[group.by]]
  cell_embeddings<-sdata@reductions$umap@cell.embeddings
  circ_data<-data.frame(Cluster,cell_embeddings)
  if(length(circsub)==1){
    if(circsub=="all"){
      subdata<-circ_data
    }
    else{
      subdata<-subset(circ_data,Cluster==circsub)
    }
  }
  else{
    subdata<-subset(circ_data,Cluster %in% circsub)
  }
  p1<-DimPlot(sdata,label=T,group.by=group.by)+
    stat_ellipse(data=subdata,aes(x=UMAP_1,y=UMAP_2,color=Cluster,fill=Cluster),level = 0.9,geom = "polygon", alpha = alpha_value,linewidth=linesize)+
    NoLegend()
  return(p1)
}

#5. convert to scanpy anndata
convert_to_scanpy<-function(sdata,fname){
  dyn.load('/home/software/hdf5-1.14.0/hdf5/lib/libhdf5_hl.so.310')
  library(SeuratDisk)
  SaveH5Seurat(sdata, filename = paste0(fname,".h5Seurat"))
  Convert(paste0(fname,".h5Seurat"), dest = "h5ad")
}

#6. enrichment analysis by gsea cell type database
library(clusterProfiler)
library(dplyr)
library(msigdbr)
library(enrichplot)
library(ggplot2)
#library(DOSE)
library(fgsea)
library(tibble)
#em <- enricher(genelist, TERM2GENE=m_t2g)
#p1 <- cnetplot(em, foldChange=gsea_list)
#pdf("tmp_heatplot.pdf",height=5,width=8)
#heatplot(em, foldChange=gsea_list, showCategory=10)
#dev.off()
#6-1 enrichment analysis
enrich_analysis<-function(markersub,title,species="Homo sapiens",catgory_type="C8",catgory_num=10,cols=c("#FF3399","grey80")){
  genelist<-markersub$gene
  logFC<-markersub$avg_log2FC
  names(logFC)<-genelist
  m_t2g <- msigdbr(species = species, category = catgory_type) %>%
    dplyr::select(gs_name, gene_symbol)
  em <- enricher(genelist, TERM2GENE=m_t2g)
  p1<-barplot(em, showCategory=catgory_num)+scale_fill_gradient(low=cols[1],high=cols[2])+labs(title=title)+
    theme(plot.title=element_text(size=25))
  p2<-heatplot(em,foldChange=logFC,showCategory = catgory_num)+scale_fill_gradient(low=cols[2],high=cols[1])+labs(title=title)+
    theme(plot.title=element_text(size=25))
  result<-list(em,p1,p2)
  return(result)
}
make_celltype_enrichment_analysis<-function(marker_file,genenum=20,fname,height_bar=6,width_bar=20,height_heat=6,width_heat=9){
  marker_file %>%
    group_by(cluster) %>%
    top_n(n=genenum,wt=avg_log2FC)->
    topn_markers
  clist<-unique(marker_file$cluster)
  enrich_list=list()
  plot_list1=list()
  plot_list2=list()
  pdfname_bar<-paste0("Celltype_enrichment_analysis_msigDB_",fname,"_barplot.pdf")
  pdfname_heat<-paste0("Celltype_enrichment_analysis_msigDB_",fname,"_heatplot.pdf")
  for(cl in clist){
    markersub<-subset(topn_markers,cluster==cl)
    titlename=paste0("Cluster ",cl," enrichment")
    result_list<-enrich_analysis(markersub,titlename)
    em<-result_list[[1]]
    p1<-result_list[[2]]
    p2<-result_list[[3]]
    enrich_list<-c(enrich_list,list(em))
    plot_list1<-c(plot_list1,list(p1))
    plot_list2<-c(plot_list2,list(p2))
  }
  pdf(pdfname_bar,height=height_bar,width=width_bar)
  for(i in 1:length(plot_list1)){
    print(plot_list1[[i]])
  }
  dev.off()
  pdf(pdfname_heat,height=height_heat,width=width_heat)
  for(i in 1:length(plot_list2)){
    print(plot_list2[[i]])
  }
  dev.off()
  result_list<-list(enrich_list,plot_list1,plot_list2)
  return(result_list)
}

#6-2 gsea analysis by clusterProfiler
gsea_analysis_by_clusterProfiler<-function(markersub,title,gene_num=20,catgory_num=10){
  #make gsea analysis input list
  genelist<-markersub$gene
  gsea_list<-markersub$avg_log2FC
  names(gsea_list)<-as.character(genelist)
  #if genelist is prepared by dplyr top_n function, this step can be omitted
  gsea_list = sort(gsea_list, decreasing = TRUE)
  #make gsea with downloaded gmt file
  ct_gmt<-read.gmt("/data/wangrj/single_cell/gmt_files/c8.all.v2023.1.Hs.symbols.gmt")
  gsea_result<-GSEA(gsea_list,TERM2GENE=ct_gmt,minGSSize=0)
  return(gsea_result)
}

#6-3 gsea analysis by fgsea, with barplot and heatmap
#markersub<-subset(marker_file,cluster==1 & p_val_adj<0.05 & avg_log2FC>0)
construct_gsea_heatmap_data<-function(markersub,fg_result,select_num=100){
  pathway_list<-fg_result$Pathway
  rownames(markersub)<-markersub$gene
  fg_genelist<-fg_result$leadingEdge
  ge_data<-c()
  for(i in 1:length(fg_genelist)){
    pway<-fg_result$Pathway[i]
    fgene<-fg_genelist[i]
    fglist<-unlist(fgene)
    for(ge in fglist){
      dtmp<-data.frame(ge,pway)
      ge_data<-rbind(ge_data,dtmp)
    }
  }
  geneall<-unique(ge_data$ge)
  ge_new<-c()
  for(gl in geneall){
    tmp1<-c(gl,markersub[gl,'avg_log2FC'],markersub[gl,'pct.1'])
    tmp2<-c(gl,0,0)
    for(pl in pathway_list){
      sub1<-subset(ge_data,pway==pl)
      gsub<-sub1$ge
      if(gl %in% gsub){
        tmp_p<-c(tmp1,pl)

      }
      else{
        tmp_p<-c(tmp2,pl)
      }
      ge_new<-rbind(ge_new,tmp_p)
    }
  }
  ge_new<-data.frame(ge_new)
  colnames(ge_new)<-c("Gene","Exp_avg_logFC","Pct","Pathway")
  ge_new$Exp_avg_logFC<-as.numeric(ge_new$Exp_avg_logFC)
  ge_new$Pct<-as.numeric(ge_new$Pct)
  #select top 150 gene by fold change
  if(length(ge_new$Gene)<=select_num){
    ge_select<-ge_new
  }
  else{
    ge_new %>%
      group_by(Gene) %>%
      summarise(Sum_result=sum(Exp_avg_logFC)) ->
      gene_summarized
    genes_summarized<-gene_summarized[order(gene_summarized$Sum_result,decreasing=T),]
    genes_select<-genes_summarized$Gene[1:select_num]
    ge_select<-subset(ge_new,Gene %in% genes_select)
  }
  return(ge_select)
}
gsea_analysis<-function(markersub,title,catgory_num=10,pval_used="padj",color_l="#FFCCCC",color_h="#FF3366",gmt_file){
  ct_pathway<-gmtPathways(gmt_file)
  Genes<-markersub$gene
  avg_logFC<-markersub$avg_log2FC
  diff_frame<-data.frame(Genes,avg_logFC)
  ranks<-deframe(diff_frame)
  fgseaRes <- fgseaMultilevel(pathways = ct_pathway, stats=ranks,minSize=3,maxSize=500,scoreType = "pos",nPermSimple = 1000)
  fgdata<-fgseaRes[order(fgseaRes$NES,decreasing=T),]
  if(pval_used=="padj"){
    fgsub<-subset(fgdata,padj<0.05 & NES>0)
  }
  else if(pval_used=="pval"){
    fgsub<-subset(fgdata,pval<0.05 & NES>0)
  }
  if(length(fgsub$NES)>catgory_num){
    fg_final<-fgsub[1:10,]
  }
  else if(length(fgsub$pathway)>0){
    fg_final<-fgsub
  }
  else if(length(fgsub$pathway)==0){
    fg_final<-fgdata[1:10,]
    color_l<-"grey100"
    color_h<-"grey80"
  }
  fg_final<-fg_final[order(fg_final$NES),]
  colnames(fg_final)[1]<-"Pathway"
  fg_final$Significance<-log2(fg_final[[pval_used]])*(-1)
  fg_final$Pathway<-factor(fg_final$Pathway,levels=as.character(fg_final$Pathway))
  #barplot
  p1<-ggplot(fg_final,aes(x=Pathway,y=NES,fill=Significance))+
    geom_bar(stat="identity")+
    scale_fill_gradient(low=color_l,high=color_h)+
    coord_flip()+theme_classic()+labs(title="",x="")+
    theme(plot.title=element_text(hjust=1,size=18),axis.title=element_text(size=15),axis.text.y=element_blank())
  #heatmap
  plot_data<-construct_gsea_heatmap_data(markersub,fg_final)
  plot_data$Pathway<-factor(plot_data$Pathway,levels=as.character(fg_final$Pathway))
  p2<-ggplot(data = plot_data, mapping = aes_string(x = "Gene", y = "Pathway")) +
    geom_tile(fill="white",color="grey50")+
    geom_point(mapping = aes_string(size = "Pct", color = "Exp_avg_logFC")) +
    scale_size_area(max_size = 4.8)+
    scale_color_gradient(low=color_l,high=color_h)+
    labs(title=title)+theme_classic()+
    theme(plot.title=element_text(hjust=0,size=18),text=element_text(size=10),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),axis.title.y=element_text(size=13))
  genelist<-fgdata$leadingEdge
  Gene<-c()
  for(i in 1:length(genelist)){
    ge<-genelist[[i]]
    ge_merge<-paste(ge,collapse=" | ")
    Gene<-c(Gene,ge_merge)
  }
  fgdata$Gene<-Gene
  gsea_new<-select(fgdata,-leadingEdge)
  result<-list(p2,p1,gsea_new)
  #design<-"AAAAAAB"
  #wrap_plots(A=p2,B=p1,design=design)
  return(result)
}
make_gsea_celltype_analysis<-function(marker_file,fname,genenum="all",height_plot=5,width_plot=20,
                                      gmt_file="/data/wangrj/single_cell/gmt_files/c8.all.v2023.1.Hs.symbols.gmt"){
  gsea_result<-c()
  markers_select<-subset(marker_file,p_val_adj<0.05 & avg_log2FC>0)
  if(genenum!="all"){
    markers_select %>%
      group_by(cluster) %>%
      top_n(n=genenum,wt=avg_log2FC)->
      topn_markers
  }
  else{
    topn_markers<-markers_select
  }
  clist<-unique(marker_file$cluster)
  pdfname<-paste0("Celltype_GSEA_analysis_msigDB_",fname,"_plots.pdf")
  pdf(pdfname,height=height_plot,width=width_plot)
  for(cl in clist){
    markersub<-subset(topn_markers,cluster==cl)
    titlename=paste0("Cluster ",cl," GSEA analysis")
    plot_list<-gsea_analysis(markersub,titlename,gmt_file=gmt_file)
    result_sub<-plot_list[[3]]
    result_sub$Cluster<-cl
    gsea_result<-rbind(gsea_result,result_sub)
    design<-"AAAAAAB"
    print(wrap_plots(A=plot_list[[1]],B=plot_list[[2]],design=design))
  }
  dev.off()
  return(gsea_result)
}

#6-4 gsva analysis for celltype annotation
library(GSVA)
library(GSEABase)
gsva_analysis_cell_type<-function(sdata,marker_file,cluster_name="RNA_snn_res.0.8"){
  diff_gene<-unique(as.vector(marker_file$gene))
  Idents(sdata)<-sdata@meta.data[[cluster_name]]
  avg_list<-AverageExpression(sdata,assays="RNA",features=diff_gene)
  avg<-avg_list[["RNA"]]
  avg<-as.matrix(avg)
  gmt<-getGmt("/data/wangrj/single_cell/gmt_files/c8.all.v2023.1.Hs.symbols.gmt")
  GSVA_result <- gsva(avg, gmt, min.sz=10, max.sz=500, verbose=FALSE, parallel.sz=1)
  return(GSVA_result)
}

#7. cell annotation by SingleR
library(SingleR)
singler_anno<-function(reffile="/data/wangrj/single_cell/ref_singlecell_data/embryo_named_230513.qs",testdata,method="cluster"){
  refdata=qread(reffile)
  ref_sce<-as.SingleCellExperiment(refdata)
  test_sce<-as.SingleCellExperiment(testdata)
  if(method=="cluster"){
    annodata<-SingleR(test=test_sce,ref=ref_sce,labels=ref_sce$comparetype,clusters=test_sce$RNA_snn_res.0.8)
    #rename idents
    Idents(testdata)<-testdata$RNA_snn_res.0.8
    Celltype<-annodata$pruned.labels
    names(Celltype)<-levels(testdata)
    testdata<-RenameIdents(testdata,Celltype)
    testdata$Celltype_by_SingleR<-Idents(testdata)
  }
  else if(method=="cell"){
    annodata<-SingleR(test=test_sce,ref=ref_sce,labels=ref_sce$comparetype)
    testdata$Celltype_by_SingleR_cell<-annodata@listData$pruned.labels
  }
  result<-list(testdata,annodata)
  return(result)
}

#plotScoreHeatmap(annodata)
#plotDeltaDistribution(annodata, ncol = 3)


#8. cell type annotation by Seurat
seurat_celltype<-function(reffile="/data/wangrj/single_cell/ref_singlecell_data/embryo_named_230513.qs",testdata,defaultassay="RNA"){
  refdata=qread(reffile)
  DefaultAssay(refdata)<-defaultassay
  DefaultAssay(testdata)<-defaultassay
  trans_anchors<-FindTransferAnchors(reference=refdata,query=testdata,dims=1:30,reference.reduction="pca")
  predictions<-TransferData(anchorset=trans_anchors,refdata=refdata$comparetype,dims=1:30)
  testdata<-AddMetaData(testdata,metadata=predictions)
  #testdata<-MapQuery(anchorset=trans_anchors,reference=refdata,query=testdata,refdata=list(celltype="comparetype"),
  #                   reference.reduction="pca",reduction.model="umap")
  return(testdata)
}

#9. subset marker files for plasma membrane or transcription factor genes
#9-1 plasma membrane
is_plasma_membrane<-function(marker_file){
  pm_list<-read.delim("/data/wangrj/single_cell/gmt_files/GO_plasma_membrane_select_genes.txt",header=F)
  colnames(pm_list)<-c("Gene_ID","Gene_symbol")
  pm_genes<-unique(pm_list$Gene_symbol)
  marker_sub<-subset(marker_file,gene %in% pm_genes)
  return(marker_sub)
}
#9-2 transcription factor
is_transcription_factor<-function(marker_file){
  tf_list<-read.delim("/data/wangrj/single_cell/gmt_files/TF_names_v_1.01.txt",header=F)
  tf_genes<unique(tf_list$V1)
  marker_sub<-subset(marker_file,gene %in% tf_genes)
  return(marker_sub)
}

#10. calculate the expression percentage of certain genes, default threshold is 0
cal_exp_percentage<-function(stem.data,features,thres=0){
  exp_data<-FetchData(stem.data,vars=features)
  colnames<-c()
  for(i in 1:length(exp_data)){
    gname<-paste0("Gene",i)
    colnames<-c(colnames,gname)
  }
  colnames(exp_data)<-colnames
  samlist=unique(stem.data$Sample)
  metadata<-stem.data@meta.data
  metadata<-cbind(metadata,exp_data)
  #make the combination list of different features
  combn_list<-list()
  numall<-length(features)
  for(i in 1:numall){
    com_each<-combn(numall,i)
    for(j in 1:length(com_each[1,])){
      com_li<-com_each[,j]
      combn_list<-c(combn_list,list(com_li))
    }
  }
  dataout<-c()
  for(sl in samlist){
    sub1<-subset(metadata,Sample==sl)
    for(cl in combn_list){
      subsub<-sub1
      for(k in cl){
        subsub$tmp<-subsub[[paste0("Gene",k)]]
        sub_tmp<-subset(subsub,tmp>thres)
        if(length(sub_tmp[,1])==0){
          subsub<-sub_tmp
          break
        }
        else{
          subsub<-sub_tmp
        }
      }
      num_each<-length(subsub[,1])
      per_each<-num_each/length(sub1[,1])*100
      fsub<-features[cl]
      fname<-paste(fsub,collapse="&")
      dtmp<-data.frame(sl,fname,per_each)
      dataout<-rbind(dataout,dtmp)
    }
  }
  colnames(dataout)<-c("Sample","Gene","Percentage(%)")
  return(dataout)
}

#11. make blend featureplot by dimplot function
make_blend_featureplot<-function(sdata,features){
  exp_matrix<-FetchData(sdata,vars=features)
  exp_type<-c()
  exp_num<-c()
  for(i in 1:length(exp_matrix[,1])){
    if(exp_matrix[i,1]>0 & exp_matrix[i,2]>0){
      exp_type<-c(exp_type,"Yes")
      exp_num<-c(exp_num,(exp_matrix[i,1]+exp_matrix[i,2]))
    }
    else{
      exp_type<-c(exp_type,"No")
      exp_num<-c(exp_num,0)
    }
  }
  sdata$exp_type<-exp_type
  ptitle=paste0(features[1]," & ",features[2])
  p1<-DimPlot(sdata,group.by="exp_type",cols=c("lightgrey","#FF0099"),split.by="Sample",pt.size=0.5)+labs(title=ptitle)
  return(p1)
}

#12. Add cluster name to the dataset by Seurat function RenameIdents by loading cell type csv file
#The file has two columns:Cluster,CellType
add_cluster_name<-function(sdata,CelltypeFileName="Celltype.csv",clustername="Cluster_name",cluster_to_annotate="RNA_snn_res.0.8"){
  ctdata<-read.csv(CelltypeFileName,header=T)
  cluster<-ctdata$Cluster
  celltype<-ctdata$CellType
  Idents(sdata)<-sdata[[cluster_to_annotate]]
  names(celltype)<-cluster
  sdata<-RenameIdents(sdata,celltype)
  sdata[[clustername]]<-Idents(sdata)
  return(sdata)
}



