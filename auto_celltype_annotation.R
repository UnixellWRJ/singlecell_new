feature_position<-c("OTX2","EN1","FGF8","STMN2","LMX1A","FOXA2","BARHL1","BARHL2","CRH","PLP1","OLIG1","OLIG2","OTP","HOXA2")
feature_midbrain<-c("LMX1A","FOXA2","PLP1","PAX7","POU4F1","NKX2-2","SIX3","SIM1")
feature_neuron<-c("SLC32A1","GAD1","GAD2","TH","ALDH1A1","KCNJ6","NR4A2","CALB1","SLC17A6","NEUROD1","NEUROD2")
celltype_names<-c("Neuron&Nb","P_dien","P_dien_PLP1","P_meten","P_meten_PLP1","P_mesenBP",
                  "P_mesenFP","P_MHB")
celltype_names<-data.frame(celltype_names)
variable_names<-c("nnb_clusters","dien_clusters","dien_PLP1_clusters","hind_clusters","hind_PLP1_clusters",
                  "mBP_clusters","mFP_clusters","mhb_clusters")
rownames(celltype_names)<-variable_names
auto_cell_annotation_vitro<-function(sdata,clustername="RNA_snn_res.0.8",clustername_auto="Cluster_name_auto_res0.8"){
  #merge the used features
  featureall<-unique(c(feature_position,feature_midbrain,feature_neuron))
  exp_data<-FetchData(sdata,vars=featureall)
  exp_data$cluster_tmp<-sdata@meta.data[[clustername]]
  #calculate the expression percentage in each cluster
  pct_data<-c()
  cnames_all=levels(sdata@meta.data[[clustername]])
  cnames=subset(cnames_all,cnames_all %in% unique(sdata@meta.data[[clustername]]))
  for(cm in cnames){
    sub1<-subset(exp_data,cluster_tmp==cm)
    exp_sub<-sub1[,1:(length(sub1)-1)]
    pct.exp <- apply(X = exp_sub, MARGIN = 2, FUN = PercentAbove, threshold = 0)
    pct_data<-rbind(pct_data,pct.exp)
  }
  rownames(pct_data)<-cnames
  pct_data<-data.frame(pct_data)
  #calculate the average expression result of selected genes
  Idents(sdata)<-sdata@meta.data[[clustername]]
  avg_list<-AverageExpression(sdata,assays="RNA",features=featureall)
  avg<-avg_list[["RNA"]]
  avg_scale<-data.frame(t(scale(t(avg),center=T,scale=T)))
  avg_scale$Gene<-rownames(avg_scale)
  colnames(avg_scale)<-gsub("X","",colnames(avg_scale))
  avg_exp<-data.frame(t(avg))
  CellType<-rep("NA",length(cnames))
  CellType<-data.frame(CellType)
  rownames(CellType)<-cnames
  
  #step1: STMN2 for Neuron&Nb
  nnb_clusters<-rownames(subset(pct_data,STMN2>0.3))
  #step2: FGF8 for MHB
  pct_data_tmp<-pct_data
  pct_data_tmp$exp_FGF8<-avg_exp$FGF8
  mhb_clusters<-rownames(subset(pct_data_tmp,FGF8>0.3 & exp_FGF8>1))
  #step3: dien, midbrain, hind, annotated by the OTX2/EN1 value.
  rowname_sub<-!rownames(pct_data) %in% c(nnb_clusters,mhb_clusters)
  pct_sub<-pct_data[rowname_sub,]
  avg_sub<-avg_exp[rowname_sub,]
  pct_sub$OTX2_to_EN1<-pct_sub$OTX2/pct_sub$EN1
  pct_sub$exp_OTX<-avg_sub$OTX
  #dien: BARHL1/BARHL2/BARHL3/CRH+, OTX2+, EN1-
  dien_sub1<-subset(pct_sub,OTX2_to_EN1>2)
  dien_sub2<-subset(dien_sub1,BARHL1>0.1 & BARHL2>0.3 & CRH>0.1 )
  dien_clusters<-rownames(dien_sub2)
  dien_PLP1_clusters<-rownames(subset(dien_sub1,PLP1>0.3))
  #hind, OTX2-EN1+
  hind_clusters<-rownames(subset(pct_sub,OTX2_to_EN1<0.5 & exp_OTX<3.5 & PLP1<0.3))
  hind_PLP1_clusters<-rownames(subset(pct_sub,OTX2_to_EN1<0.5 & exp_OTX<3.5 & PLP1>0.3))
  #mid,OTX2+EN1+
  mid_rownames<-!rownames(pct_sub) %in% c(dien_clusters,dien_PLP1_clusters,hind_clusters,hind_PLP1_clusters)
  mid_clusters<-rownames(pct_sub)[mid_rownames]
  #mid_clusters<-rownames(subset(pct_sub,OTX2_to_EN1>0.5 & OTX2_to_EN1<2))
  #step 4: annotation of mid_clusters
  pct_mid<-pct_data[mid_clusters,]
  pct_mid$FOXA2_to_LMX1A<-pct_mid$FOXA2/pct_mid$LMX1A
  mFP_clusters<-rownames(subset(pct_mid,FOXA2_to_LMX1A<2 | LMX1A>0.2))
  mBP_clusters<-rownames(subset(pct_mid,FOXA2_to_LMX1A>2 & LMX1A<0.2))
  
  #mFP_clusters<-rownames(subset(pct_mid,FOXA2_to_LMX1A<2 | LMX1A>0.2 & PLP1<0.3))
  #mBP_clusters<-rownames(subset(pct_mid,FOXA2_to_LMX1A>2 & LMX1A<0.2 & PLP1<0.3))
  #mFP_PLP1_clusters<-rownames(subset(pct_mid,FOXA2_to_LMX1A<2 | LMX1A>0.2 & PLP1>0.3))
  #mBP_PLP1_clusters<-rownames(subset(pct_mid,FOXA2_to_LMX1A>2 & LMX1A<0.2 & PLP1>0.3))
  #ct_list<-list(get(variable_names))
  for(vn in variable_names){
    cl<-get(vn)
    CellType[cl,1]<-celltype_names[vn,1]
  }
  Idents(sdata)<-sdata@meta.data[[clustername]]
  celltype_anno<-CellType$CellType
  names(celltype_anno)<-levels(sdata)
  sdata<-RenameIdents(sdata,celltype_anno)
  sdata@meta.data[[clustername_auto]]<-Idents(sdata)
  return(list(celltype_anno,sdata))
}








