###
library(Seurat)
library(cowplot)
library(Matrix)
library(qs)
library(future)
library(DoubletFinder)

#own_dir: dir to save the output results
#batchname: the name used for output file
#filelist, namelist: files used for integration, and the names used to identify samples

#ES & IPS batch
es_ips_list=list("/home/wrj/wangrj/single_cell/all_matrix/RM_doublets/ES_IPS/","ES_IPS_with_add_no_doublets",c("d21_add","YW3D21_add","ZYYD21","time_D28","YW3D28","ZYYD28_add","00629","00630"),
                 c("ES_H9_D21","IPS_YW_D21","IPS_ZYY_D21","ES_H9_D28","IPS_YW_D28","IPS_ZYY_D28","IPS_006Line29_D28","IPS_006Line30_D28"))

fname_list=paste0("stem_list_",batchname,".qs")
fname_anchor=paste0("stem_anchors_",batchname,".qs")
fname_rna=paste0("stem_data_integrated_",batchname,"_rna.qs")
fname_integrated=paste0("stem_data_integrated_",batchname,"_integrated.qs")
fname_final=paste0("stem_data_final_",batchname,".qs")
fname_marker=paste0("stem_data_all_markers_",batchname,"_rna.qs")
options(future.globals.maxSize = 1000000 * 1024^2)

def_ratio<-function(number){
  if(number>=10000){
    ratio<-0.076
  }
  else if(number>=9000){
    ratio<-0.069
  }
  else if(number>=8000){
    ratio<-0.061
  }
  else if(number>=7000){
    ratio<-0.054
  }
  else if(number>=6000){
    ratio<-0.054
  }
  else if(number>=5000){
    ratio<-0.039
  }
  else if(number>=4000){
    ratio<-0.031
  }
  else if(number>=3000){
    ratio<-0.023
  }
  else if(number>=2000){
    ratio<-0.016
  }
  else if(number>=1000){
    ratio<-0.008
  }
  else if(number>=800){
    ratio<-0.004
  }
  else{
    ratio<-0
  }
  return(ratio)
}

# filtering doublets
rm_doublet <- function(sample,ndims=10){
  object<-sample
  #define the ratio number
  cellnums<-length(colnames(object))
  ratio<-def_ratio(cellnums)
  sweep.res.list <- paramSweep_v3(object, PCs = 1:ndims, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  nExp_poi <- round(ratio*nrow(object@meta.data))
  object <- doubletFinder_v3(object, PCs = 1:ndims, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  ob_meta<-object@meta.data
  sample$Def_doublets<-ob_meta[,length(ob_meta)]
  sample_sub<-subset(sample,Def_doublets=="Singlet")
  return(sample_sub)
}

#select proper dim number according to elbow results, i.e., the results from "Stdev" function
#choose dim where the elbow slope is relatively stable. The cutoff is set to 0.05 to ensure the 
#reuslt is relatively stable. This number can be changed for proper reason.
select_dim<-function(sample, avg_cutoff=0.05){
  data_std<-Stdev(sample,reduction="pca")
  data_tmp<-c(data_std[2:length(data_std)],data_std[length(data_std)])
  data_diff<-data_std-data_tmp
  num<-length(data_std)-5
  for(i in 1:num){
    avg<-mean(data_diff[i:(i+5)])
    if(avg<avg_cutoff){
      break
    }
  }
  dimnum=i+3
  return(dimnum)
}

#data reading and filtering
#transfer the 10X data to Seurat object, remove doublets, and filter low activate cells.
read_and_filter_10X<-function(dirname,samplename){
  sample.raw<-Read10X(dirname)
  sample<-CreateSeuratObject(counts=sample.raw,min.cells=3)
  sample<-NormalizeData(sample,verbose=FALSE,scale.factor = 40000)
  sample<-FindVariableFeatures(sample,selection.method="vst",nfeatures=2000)
  sample<-ScaleData(sample)
  sample<-RunPCA(sample)
  #Find proper dims
  dimnum<-select_dim(sample)
  sample<-RunUMAP(sample, dims = 1:dimnum)
  print(dimnum)
  sample_no_doublet<-rm_doublet(sample,ndims=dimnum)
  #hard filtering
  sample_no_doublet[["percent.mt"]] <- PercentageFeatureSet(sample_no_doublet, pattern = "^MT-")
  sample_sub<-subset(sample_no_doublet,subset = nCount_RNA > 1000 & nCount_RNA < 30000 & 
                       nFeature_RNA>500 & nFeature_RNA<6000 & 
                       percent.mt < 8)
  sample_sub$Sample<-samplename
  return(sample_sub)
}

#merge ES and IPS data
merge_data<-function(argslist){
  own_dir=argslist[1]
  batchname=argslist[[2]]
  filelist=argslist[[3]]
  namelist=argslist[[4]]
  stem.list<-list()
  for(i in 1:length(filelist)){
    fl=filelist1[i]
    dirname=paste0("/home/wrj/wangrj/single_cell/all_matrix/all_matrix_raw/",fl,"/outs/filtered_feature_bc_matrix/")
    samplename<-namelist[i]
    sample<-read_and_filter_10X(dirname,samplename)
	  stem.list<-c(stem.list,list(sample))
  }
  qsave(stem.list,paste0(own_dir,fname_list))
  features <- SelectIntegrationFeatures(object.list = stem.list)
  plan("multisession", workers = 4)
  #use RPCA for reduction instead of CCA. According to the tutorial, RPCA is more suited for samples with distinct cell types. 
  stem.anchors<-FindIntegrationAnchors(object.list=stem.list,anchor.features=features,dims=1:40,reduction = "rpca")
  qsave(stem.anchors,paste0(own_dir,fname_anchor))
  stem.data<-IntegrateData(anchorset=stem.anchors,dims=1:40)
  qsave(stem.data,paste0(own_dir,fname_rna))

  DefaultAssay(stem.data) <- "integrated"
  stem.data <- ScaleData(stem.data, verbose = FALSE)
  stem.data <- RunPCA(stem.data, verbose = FALSE, npcs=100)
  dim_all<-select_dim(stem.data)
  stem.data <- RunUMAP(stem.data, dims = 1:dim_all)
  stem.data <- FindNeighbors(stem.data, reduction = "pca", dims = 1:dim_all)
  stem.data <- FindClusters(stem.data, resolution = 0.8)
  qsave(stem.data,paste0(own_dir,fname_integrated))

  #filter data
  #according to the VariableFeaturePlot, there are not much difference between 2000, 4000 or 8000. 
  #A number of the proliferation-related genes need to be delete, so choose 4000. 
  DefaultAssay(stem.data)<-"RNA"
  stem.data<-FindVariableFeatures(stem.data,nfeatures = 4000)
  var_features<-stem.data@assays$RNA@var.features
  exp_mat<-FetchData(stem.data,vars=var_features)
  top_2a<-exp_mat$TOP2A
  var_sub<-c()
  corlist<-c()
  ccycle_features<-c()
  for(vf in var_features){
    exp_ge<-exp_mat[[vf]]
    cor_result<-cor(top_2a,exp_ge,method="pearson")
    corlist<-c(corlist,cor_result)
    if(cor_result<=0.15){
      var_sub<-c(var_sub,vf)
    }
    else{
      ccycle_features<-c(ccycle_features,vf)
    }
  }
  stem.data@assays$RNA@var.features<-var_sub
  s.genes <- cc.genes.updated.2019$s.genes
  g2m.genes <- cc.genes.updated.2019$g2m.genes
  stem.data <- CellCycleScoring(object = stem.data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  stem.data$CC.Difference <- stem.data$S.Score - stem.data$G2M.Score
  #stem.data<-ScaleData(stem.data,vars.to.regress="percent.mt")
  #stem.data <- ScaleData(object = stem.data, 
  #                    vars.to.regress = c("nCount_RNA",'nFeature_RNA','CC.Difference','percent.mt'),verbose = TRUE)
  stem.data<-ScaleData(stem.data)
  stem.data<-RunPCA(stem.data,npcs=100)
  ndim<-select_dim(stem.data)
  print(ndim)
  stem.data<-FindNeighbors(stem.data, reduction = "pca", dims = 1:ndim, nn.eps = 0.5)
  stem.data <- FindClusters(stem.data, resolution = 0.8, n.start = 10)
  stem.data<-RunUMAP(stem.data,reduction="pca",dims=1:ndim)
  stem.data<-RunTSNE(stem.data,reduction="pca",dims=1:ndim)
  qsave(stem.data,paste0(own_dir,fname_final))
  plan("multisession",workers=4)
  all_markers<-FindAllMarkers(stem.data,min.pct=0.25)
  write.csv(all_markers,paste0(own_dir,fname_marker),quote=F)
}

#main
merge_data(es_ips_list)
