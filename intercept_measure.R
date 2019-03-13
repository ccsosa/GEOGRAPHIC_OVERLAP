
##########################################'
require(raster)

gap_dir <- "D:/ccsosa/gap_cajanusPaper"
src_dir <- paste0(gap_dir,"/","_scripts")
source(paste0(src_dir,"/","000.zipWrite.R"))
source(paste0(src_dir,"/","000.zipRead.R"))

ENM_dir <- paste0(gap_dir,"/","Geographic_overlap")

models_dir <- paste0(gap_dir,"/","Geographic_overlap","/","ENM")
models_files <- list.files(models_dir,pattern="gz$")
models_names <- models_files
M_NAM <- gsub(pattern=("_worldclim2_5_EMN_PA_C.asc.gz"),x=as.character(models_names),replacement ="")
M_NAM <- gsub(pattern=("_samples-buffer-na.asc.gz"),x=as.character(M_NAM),replacement ="")

matrix_sp <- matrix(nrow=length(M_NAM),ncol=length(M_NAM))
second <- matrix(nrow=length(M_NAM),ncol=length(M_NAM))
name_sp <- matrix(nrow=length(M_NAM),ncol=length(M_NAM))
#matrix_sp <- as.data.frame(matrix_sp)
colnames(matrix_sp) <- M_NAM
rownames(matrix_sp) <- M_NAM

colnames(second) <- M_NAM
rownames(second) <- M_NAM

colnames(name_sp) <- M_NAM
rownames(name_sp) <- M_NAM

for(temp in 1:length(models_files)){
  
  #template <- raster(paste0(ENM_dir,"/","Helianthus_annuus_worldclim2_5_EMN_PA_C.asc"))
  template <- zipRead(models_dir,models_files[[temp]])
  #template[which(template[]==0)] <- NA
  template_cells <- sum(template[]==1,na.rm =T)
  
  cat("processing", as.character(models_files[[temp]]),"\n")
  meas <- matrix(nrow=15,ncol=1)
  sec <- matrix(nrow=15,ncol=1)
  names <- matrix(nrow=15,ncol=1)
  
  for(i in 1:length(models_files)){
    x <- zipRead(models_dir,models_files[[i]])
    #x[which(x[]==0)] <- NA
    xcells <- sum(x[]==1,na.rm =T)
    #y <- intersect(x,template)
    y=x*template
    inter_cells <- sum(y[]==1,na.rm =T)
    
    cat("processing", as.character(models_files[[temp]]),"vs",as.character(models_files[[i]]),"\n")
    
    Measure <- 2*(inter_cells/(xcells+template_cells))
    
    if(template_cells<xcells){
      n <- gsub(pattern=("_worldclim2_5_EMN_PA_C.asc.gz"),x=models_files[[temp]],replacement ="")
      n <- gsub(pattern=("_samples-buffer-na.asc.gz"),x=n,replacement ="")
      cat("minor range is ",n,"\n")
      s <- 2*(template_cells/(xcells+template_cells))
      s <- Measure/s
    }else{
      s <- 2*(xcells/(xcells+template_cells))
      s <- Measure/s
      n <- gsub(pattern=("_worldclim2_5_EMN_PA_C.asc.gz"),x=models_files[[i]],replacement ="")
      n <- gsub(pattern=("_samples-buffer-na.asc.gz"),x=n,replacement ="")
      cat("minor range is ",n,"\n")}
  
  #Measure <- (inter_cells/template_cells)
  
  # Measure=2*Measure
  
  #cat("shared cells:",inter_cells,"\n")
  #cat(as.character(models_files[[temp]])," cells are: ",template_cells,"\n")
  #cat(as.character(models_files[[i]])," cells are: ",xcells,"\n")
  cat("overlap measure is ",Measure,"\n")
  cat("second overlap measure is ",s,"\n")
  meas[[i]] <- Measure
  sec[[i]] <- s
  names[[i]] <- n
  
}

matrix_sp[,temp] <- meas
second[,temp] <- sec
name_sp[,temp] <- names

#rm(meas,Measure,x,xcells,y,inter_cells,template,template_cells)

#print(matrix_sp)
}

#matrix_sp
write.csv(matrix_sp,paste0(gap_dir,"/","Geographic_overlap.csv"))
write.csv(second,paste0(gap_dir,"/","minor_Geo_overlap.csv"))
write.csv(name_sp,paste0(gap_dir,"/","minor_names_Geo_overlap.csv"))

#heatmap(matrix_sp)

#clu <- 1-matrix_sp

#tre <- hclust(as.dist(clu))

