
#### 1.PREPARATION ####
#### 1.1 Loading packages (not all of them are for this code) #### 

libraries=c("readxl","WriteXLS","plyr","reshape2","ggplot2","tidyr","Hmisc","corrplot","dendextend","missMDA","FactoMineR","colorspace","gplots","factoextra","cluster","ggrepel")
install_this=libraries[which(lapply(libraries, require, character.only = TRUE) == F)]
install.packages(install_this) #installing what could not be loaded
install_this=libraries[which(lapply(libraries, require, character.only = TRUE) == F)]

#### 1.2 Cleaing the workspace ####
rm(list=ls(all=TRUE)) 

#### 1.3 Setting working directory ####

data_wd = getwd() # complete pathway of the directory containing the excel database
file_name="Guilayn et al. 2017 INPUT DATA.xlsx" ###name of database excel file inside data_wd
wd = getwd() # complete pathway of working directory for saving plots and output tables
setwd(wd)

##### 1.4 Creading the subsets defined in the paper (Table 1) ####

states=c("all", "raw","liq", "sld")
all=c("RW","RD","P_SF","LF","SF")
raw=c("RW","RD","P_SF")
liq="LF"
sld="SF"
liststates=list(all,raw,liq,sld)
statedescription=c("all digestate and separation fractions from wet and dry AD", 
                   "raw digestate (whole)",
                   "liquid fraction of digestate from wet and dry AD",
                   "solid fraction of digestate from wet and dry AD"
                   )
#### 1.5 Creating directories for file exportation ####

for (i in 1:length(states))
{
  temp=states[i]
  dir.create(paste(states[i]))
  dir.create(paste0(wd,"/",states[i],"/Initial_Boxplots"))
  dir.create(paste0(wd,"/",states[i],"/SHAPIRO"))
  dir.create(paste0(wd,"/",states[i],"/Corr_matrix"))
  dir.create(paste0(wd,"/",states[i],"/HCA"))
  dir.create(paste0(wd,"/",states[i],"/Optimal_n_cluster"))
  dir.create(paste0(wd,"/",states[i],"/PCA"))
  dir.create(paste0(wd,"/",states[i],"/Boxplots"))
  dir.create(paste0(wd,"/",states[i],"/Boxplots_LEGISLATION"))
}

#### 1.6 Importing databases from excel ####
#### 1.6.1. Main database with data for statistics and comparison to legistlation####
datainput = read_excel(path = file_name,
                       sheet="DATA_STATISTICS", 
                        na = "NA",
                       col_names = TRUE)

datainput=as.data.frame(datainput)
which(lapply(datainput, class) != "numeric") ###VERIFY IF NUMERIC DATA IS IMPORT AS NUMERIC

datainput.othervar=read_excel(path = file_name,
                              sheet="DATA_LEGISLATION",
                              na = "NA",
                              col_names = TRUE)

datainput.othervar=as.data.frame(datainput.othervar)
which(lapply(datainput.othervar, class) != "numeric") ###VERIFY IF NUMERIC DATA IS IMPORT AS NUMERIC

#### 1.7. Make sure every ID is different ####
print(if (length(unique(datainput$ID))==nrow(datainput)) {
  "OK" 
} else 
  {
  repetid=as.data.frame(table(datainput$ID))[(as.data.frame(table(datainput$ID))[,2]>1),]
  repetid=data.frame("ID"=repetid[,1],"Freq"=repetid[,2])
  repetid
})

#### 1.8 Defining which columns are used for the correlation matrix, HCA and PCA ####

##1st column of data with AD process input for correlation matrix
St.cor.det="RT.d" #Start of process input data. Used for correlation matrix. No NA problem. Everything from here must be numeric

##1st column of data with digestate characteristics for correlation matrix
St.charac.det="Cd" #Start of digestate charact data. Used for correlation matrix. No NA problem. Everything from here must be numeric

##1st column of numerical data for HCA and PCA
St.PCA.det="C/N" #Start of digestate charact data used for HCA and PCA. Avoid NAs as much as possible. Everything from here must be numeric

#### 1.9. Defining max. number of NAs per line of data ####
#### Lines with more than N.NA will be excluded from the PCA and HCA analysis
N.NA=1   ###limit of NA values per line (for the PCA only)

#### 1.10 Selecting datasources and exluding outliers (if necessary), as defined in the database ####
datasource=c(
  "P", #P=Probiotic
  "D",  #D=DIVA
  "S",  #S=SUEZ
  "L",  #L=LBE
  "B")  #B=BIBLIO: including probiotic

#### Excluding digestate lines according to add, usefull for outliers ####
excludeID=c() 
          
excludeIDdf=data.frame("Excluded ID"=excludeID)
datainput=subset.data.frame(datainput, !(datainput$ID %in% excludeID))
datainput=subset.data.frame(datainput, datainput$Source %in% datasource)
#exporting list of excluded digestates
write.csv2(x=excludeIDdf, file="excluded_IDs.csv",row.names = TRUE)

#### 2. Manipulating data for statistical analysis ####
#### 2.1. Creating data frames that will be used latter ####
for (i in 1:length(states))
{
  print(paste0("running i=",i," /",length(states)))
  temp=states[i]
  
  datainput.temp=subset.data.frame(datainput,datainput$State %in% liststates[[i]])
    assign(paste0("datainput.",temp),datainput.temp)
  
  
  St.cor=which(colnames(datainput.temp)==St.cor.det)    
  St.charac=which(colnames(datainput.temp)==St.charac.det)
  St.PCA=which(colnames(datainput.temp)==St.PCA.det)  
  
  dataPCA.temp=datainput.temp[,St.PCA:ncol(datainput.temp)]
  datainputNA.temp=datainput.temp[which(rowSums(is.na(dataPCA.temp))<=N.NA),]
      assign(paste0("datainputNA",".",temp), datainputNA.temp)

  
  assign(paste0("dataCORR_ch.",temp),datainput.temp[,St.charac:ncol(datainput.temp)])       
  dataPCA.temp=datainputNA.temp[,St.PCA:ncol(datainput.temp)]
  # Removing columns with over 1/3 NAs (not possible to perform statistics if greater than 2/3) 
  remove_col = as.numeric(which(apply(dataPCA.temp,2,function(x)sum(is.na(x)))>nrow(dataPCA.temp)*1/3))
  if (sum(remove_col) > 0) {
  dataPCA.temp=dataPCA.temp[,-remove_col]
  }
    assign(paste0("dataPCA.",temp),dataPCA.temp)
}

##### 2.2 Counting substition of NAs by parameter ####

for (i in 1:length(states))
{
  print(paste0("running i=",i," /",length(states)))
  temp=states[i]
dataPCA.temp=as.data.frame(get(paste0("dataPCA.",temp)))
NAtable.temp=dataPCA.temp[1:3,]
rownames(NAtable.temp)=c("NA.count","NA.percent","N.individ")

for (j in 1:ncol(dataPCA.temp))
{
  NAtable.temp[1,j]=length(which(is.na(dataPCA.temp[,j])))
  NAtable.temp[2,j]=100*as.numeric(NAtable.temp[1,j])/nrow(dataPCA.temp)
  NAtable.temp[3,j]=nrow(dataPCA.temp)
  
}
assign(paste0("NAtable.",temp),NAtable.temp)
print(NAtable.temp)
}

#### 2.3. Exporting NA count and input manipulated data to excel ####

writeinputdata=1:length(states)
for (i in 1:length(states))
{
  writeinputdata[i]=paste0("datainputNA.",states[i])
  }
WriteXLS(c("excludeIDdf",writeinputdata), ExcelFileName = "dataPCA.xls", SheetNames = c("ExcludedIDs",states), perl = "perl",
         verbose = FALSE, Encoding = c("UTF-8", "latin1", "cp1252"),
         row.names = TRUE, col.names = TRUE,
         AdjWidth = FALSE, AutoFilter = TRUE, BoldHeaderRow = TRUE,
         na = "",
         FreezeRow = 1, FreezeCol = 0,
         envir = parent.frame())

writeNAtables=1:length(states)
for (i in 1:length(states))
{
  writeNAtables[i]=paste0("NAtable.",states[i])
}

source.summary=ddply(datainputNA.all,.(Source, State),plyr::summarize, nb=length(Source))

WriteXLS(c("source.summary",writeNAtables), ExcelFileName = "NA_subst_for_PCA.xls", SheetNames = c("Source_summary",states), perl = "perl",
         verbose = FALSE, Encoding = c("UTF-8", "latin1", "cp1252"),
         row.names = TRUE, col.names = TRUE,
         AdjWidth = FALSE, AutoFilter = TRUE, BoldHeaderRow = TRUE,
         na = "",
         FreezeRow = 1, FreezeCol = 0,
         envir = parent.frame())

############################################## STOP HERE #################################################################################################
################################ Evaluate NA substitution before proceeding ###############################################################################
################################### Redefine HCA/PCA variables if necessary ######################################################################################################################

#### 3. BASIC STATISTICS ####
#### 3.1. Looking for outliers by simple boxplots ####
for (i in 1:length(states))
{
  
  print(paste0("running i=",i," /",length(states)))
  temp=states[i]
  datainputNA.temp=get(paste0("datainputNA.",temp))
  dataPCA.temp=as.data.frame(get(paste0("dataPCA.",temp)))
  pdf(paste0(wd,"/",temp,"/Initial_Boxplots/",states[i],"_all_variables_boxplot.pdf"))
  boxplot(dataPCA.temp, use.cols = TRUE)
  dev.off()
  
  pdf(paste0(wd,"/",temp,"/Initial_Boxplots/",states[i],"_all_variables2_boxplot.pdf"))
  print(ggplot(data=melt(dataPCA.temp), aes(factor(1), value)) +
  geom_boxplot() +  facet_wrap(~variable,scales="free",nrow = 2, ncol = NULL)+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(), 
        axis.title.y = element_blank()))
  dev.off()
  
  pdf(paste0(wd,"/",temp,"/Initial_Boxplots/",states[i],"_all_variables_state_boxplot.pdf"))
  print(ggplot(data=melt(data=as.data.frame(c(datainputNA.temp[c(2,4)],dataPCA.temp)),id.vars = c("ID.PCA","State")), aes(State, value)) +
          geom_boxplot(aes(fill=State)) +  facet_wrap(~variable,scales="free",nrow = 2, ncol = NULL)+
          theme(axis.text.x=element_blank(),
                axis.title.x=element_blank(), 
                axis.title.y = element_blank()))
  dev.off()
}
#### 3.2. Normality: shapiro test and qq-plots #####

for (i in 1:length(states)) {
print(paste0("running i=",i," /",length(states)))
temp=states[i]
dataPCA.temp=get(paste0("dataPCA.",temp))

shap.temp = lapply(dataPCA.temp, shapiro.test)

shapresult.temp = sapply(shap.temp, `[`, c("statistic","p.value"))
shapresult.temp=t(shapresult.temp)
assign(paste0("shapresult.",temp),shapresult.temp)

for (j in 1:ncol(dataPCA.temp)) {
  vartemp = colnames(dataPCA.temp)[j]
  vartemp = gsub("/",".",vartemp)
  pdf(paste0(wd,"/",temp,"/SHAPIRO/",
             toupper(temp),"_",vartemp,"_QQplot.pdf"))
  qqnorm(dataPCA.temp[,j], main = paste0("Normal Q-Q Plot - ",colnames(dataPCA.temp)[j]),
         xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
         plot.it = TRUE, datax = FALSE)
  qqline(dataPCA.temp[,j])
dev.off() 
pdf(paste0(wd,"/",temp,"/SHAPIRO/",
           toupper(temp),"_",vartemp,"_HISTplot.pdf"))
hist(dataPCA.temp[,j], main = paste0("Hist. - ",colnames(dataPCA.temp)[j]))
       dev.off()
}
}

#### 3.3 Correlation matrix, p-value and plots ####

for (i in 1:length(states))
{
  print(paste0("running i=",i," /",length(states)))
  temp=states[i]
  
  dataCORR_ch.temp=get(paste0("dataCORR_ch.",temp))


  matrrcor_ch.temp=rcorr(as.matrix(dataCORR_ch.temp), type=c("pearson","spearman"))

  assign(paste0("matrrcor_ch.",temp),matrrcor_ch.temp)
  assign(paste0("r2_ch.",temp),as.data.frame(matrrcor_ch.temp$r))
  assign(paste0("pv_ch.",temp),as.data.frame(matrrcor_ch.temp$P))
  assign(paste0("nb.pairs_ch",temp),as.data.frame(matrrcor_ch.temp$n))

  
  write.corr_ch.matrix.data=c(paste0("r2_ch.",temp),paste0("pv_ch.",temp),paste0("nb.pairs_ch",temp))
  # Exporting excel
  WriteXLS(write.corr_ch.matrix.data, 
           ExcelFileName = paste(wd,"/",temp,"/Corr_matrix/",toupper(temp),"_correlation_matrix_charact_data.xls",sep=""), 
           SheetNames = c("r2","p_value","n"), perl = "perl",
           verbose = FALSE, Encoding = c("UTF-8", "latin1", "cp1252"),
           row.names = TRUE, col.names = TRUE,
           AdjWidth = FALSE, AutoFilter = TRUE, BoldHeaderRow = TRUE,
           na = "",
           FreezeRow = 1, FreezeCol = 0,
           envir = parent.frame())

  #### plotting ####
  col = colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  
  pdf(paste0(wd,"/",temp,"/Corr_matrix/",toupper(temp),"_corrmatrix.dig.charac.only_pv01.pdf"))
  corrplot(matrrcor_ch.temp$r, method="color", col=col(200),
           #mar = c(0, 0,0, 0),
           type="lower", #order="hclust", hclust.method = "average", 
           addCoef.col = "black", # Ajout du coefficient de corr?lation
           number.cex= 7/ncol(matrrcor_ch.temp$r),
           tl.col="black", tl.srt=45, #Rotation des etiquettes de textes
           tl.cex = 9/ncol(matrrcor_ch.temp$r), cl.cex = 9/ncol(matrrcor_ch.temp$r),
           # Combiner avec le niveau de significativit?
           p.mat = matrrcor_ch.temp$P, sig.level = 0.01, insig = "blank", 
           # Cacher les coefficients de corr?lation sur la diagonale
           diag=FALSE,
           na.label="NA", na.label.col="grey")
  title("p-value<0.01", line=2)
  dev.off()
  
  pdf(paste0(wd,"/",temp,"/Corr_matrix/",toupper(temp),"_corrmatrix.dig.charac.only_pv05.pdf"))
  corrplot(matrrcor_ch.temp$r, method="color", col=col(200),
           type="lower", #order="hclust", hclust.method = "average", 
           addCoef.col = "black", # Ajout du coefficient de corr?lation
           number.cex= 7/ncol(matrrcor_ch.temp$r),
           tl.col="black", tl.srt=45, #Rotation des etiquettes de textes
           tl.cex = 9/ncol(matrrcor_ch.temp$r), cl.cex = 9/ncol(matrrcor_ch.temp$r),
           # Combiner avec le niveau de significativit?
           p.mat = matrrcor_ch.temp$P, sig.level = 0.05, insig = "blank", 
           # Cacher les coefficients de corr?lation sur la diagonale
           diag=FALSE,
           na.label="NA", na.label.col="grey")
  title("p-value<0.05", line=2)
  dev.off()
  
}

#### 4. Advanced statistics ####
#### 4.1. Imputing NAs ####
for (i in 1:length(states))
{
  print(paste0("running i=",i," /",length(states)))
  temp=states[i]
  datainputNA.temp=get(paste0("datainputNA.",temp))
  dataPCA.temp=get(paste0("dataPCA.",temp))
  
  # NA substitution #
  
  nb.temp=estim_ncpPCA(dataPCA.temp, scale=TRUE, method.cv = "Kfold")    #estimation of ncp (=nb$ncp)
  assign(paste0("nb.",temp),nb.temp)
  dataPCA.complete.temp=imputePCA(dataPCA.temp, ncp=nb.temp$ncp)
  if (class(dataPCA.complete.temp) != "data.frame") {
    dataPCA.complete.temp=dataPCA.complete.temp$completeObs
  }
  assign(paste0("dataPCA.complete.",temp),dataPCA.complete.temp)
}
#### 4.2. CLUSTERING ####
#### 4.2.1. HCA and heatmaps: Inital plots  ####
## FIRST HCA PLOTS --> PLOT AND ANALYZE
for (i in 1:(length(states))) 
{
  print(paste0("running i=",i," /",length(states)))
  temp=states[i]
datainputNA.temp=get(paste0("datainputNA.",temp))
dataPCA.temp=get(paste0("dataPCA.",temp))
hcmethod="ward.D2" 

d.temp=dist(as.matrix(scale(dataPCA.temp, center=TRUE, scale=TRUE)))
 assign(paste0("d.",temp),d.temp)
hc.temp=hclust(d.temp, method=hcmethod)

if (hc.temp$method=="ward.D" | hc.temp$method=="ward.D2"){
  hc.temp$height=sqrt(hc.temp$height)
}
assign(paste0("hc.",temp),hc.temp)

# STANDARD LABEL: State, Feedstock category and ID
  nbg.temp=8
  dend.temp = as.dendrogram(hc.temp)
  dend.temp= color_branches(dend.temp, k=nbg.temp) 

  # Colors of the labels: helpful to analyse the clusters
  if(length(liststates[[i]])>1) {
    labelcolor="by digestate state"
    categfactor.temp=factor(datainputNA.temp$State)
    levels(categfactor.temp)=sort(unique(datainputNA.temp$State))
    labels_colors(dend.temp) = rainbow_hcl(length(unique(datainputNA.temp$State)))[(as.numeric(categfactor.temp))[order.dendrogram(dend.temp)]] 
    }
  
  # Hanging up the dendogram
  dend.temp = hang.dendrogram(dend.temp,hang_height=0.1)
  
  # Size of the ind labels
  dend.temp = set(dend.temp, "labels_cex", (0.15+10.5/nrow(datainputNA.temp)))
  
   # Definition of the individual labels #
  dend.temp2=dend.temp
  dend.temp3=dend.temp
  labels(dend.temp)=paste0(datainputNA.temp$State[order.dendrogram(dend.temp)],
                           "-",datainputNA.temp$Categ[order.dendrogram(dend.temp)],
                           "(",datainputNA.temp$ID[order.dendrogram(dend.temp)],")")
  labels(dend.temp2)=paste0(datainputNA.temp$State[order.dendrogram(dend.temp2)],
                            "__SepM:",datainputNA.temp$Sep.meth[order.dendrogram(dend.temp2)],
                            "__RT:",datainputNA.temp$RT.d[order.dendrogram(dend.temp2)],
                            "__Tmp:",datainputNA.temp$Temp.oC[order.dendrogram(dend.temp2)],
                            "(",datainputNA.temp$ID[order.dendrogram(dend.temp2)],")")
  labels(dend.temp3)=paste0(datainputNA.temp$State[order.dendrogram(dend.temp3)],
                            "-",datainputNA.temp$Feed.description[order.dendrogram(dend.temp3)],
                            "(",datainputNA.temp$ID[order.dendrogram(dend.temp3)],")")
 
  # Ploting and saving
  
  pdf(paste0(wd,"/",temp,"/HCA/",toupper(temp),"_HCA_meth.",hc.temp$method,"_CRUDE_INITIAL.pdf"),
      width = 7, height = 8.75)
  plot(hc.temp, labels=datainputNA.temp$ID.PCA, cex=(0.15+10.5/length(datainputNA.temp$ID.PCA)))
  dev.off()
  
  pdf(paste0(wd,"/",temp,"/HCA/",toupper(temp),"_HCA_meth.",hc.temp$method,"_INITIAL.pdf"),
      width = 7, height = 8.75)
      par(mar = c(3,3,3,7))  
      plot(dend.temp, 
       main=paste0("Clustered data set for: ",paste(liststates[[i]], collapse=", "),"\n",
                   " (",statedescription[i],")","\n",
                   "method=",hc.temp$method,"\n") ,
       cex.main=1.0, cex.sub=1,
       horiz =  TRUE,  nodePar = list(cex = .007))
      legend(title="Label colors (not cluster)",
             "topleft", cex=0.5,
             legend = levels(categfactor.temp), 
             fill = rainbow_hcl(length(levels(categfactor.temp))))
    dev.off()
  
    pdf(paste0(wd,"/",temp,"/HCA/",toupper(temp),"_HCA_meth.",hc.temp$method,"_INITIAL_DETAIL.pdf"),
        width = 7, height = 8.75)
    par(mar = c(3,3,3,7))  
    plot(dend.temp2, 
       main=paste0("Clustered data set for: ",paste(liststates[[i]], collapse=", "),"\n",
                   " (",statedescription[i],")",
                   "\n","method=",hc.temp$method) ,
       cex.main=1.0, cex.sub=1,
       horiz =  TRUE,  nodePar = list(cex = .007))
    legend(title="Label colors (not cluster)",
         "topleft", cex=0.5,
         legend = levels(categfactor.temp), 
         fill = rainbow_hcl(length(levels(categfactor.temp))))
  dev.off()
  
  pdf(paste0(wd,"/",temp,"/HCA/",toupper(temp),"_HCA_meth.",hc.temp$method,"_INITIAL_feed_description.pdf"),
      width = 7, height = 8.75)
  par(mar = c(3,3,3,7))  
     plot(dend.temp3, 
       main=paste0("Clustered data set for: ",paste(liststates[[i]], collapse=", "),"\n",
                   " (",statedescription[i],")",
                   "\n","method=",hc.temp$method) ,
       cex.main=1.0, cex.sub=1,
       horiz =  TRUE,  nodePar = list(cex = .007))
     legend(title="Label colors (not cluster)",
         "topleft", cex=0.5,
         legend = levels(categfactor.temp), 
         fill = rainbow_hcl(length(levels(categfactor.temp))))
  dev.off()
  
 
  pdf(paste0(wd,"/",temp,"/HCA/",toupper(temp),"_HCA_meth.",hc.temp$method,"_HEATMAP_INITIAL.pdf"),
      width = 7, height = 8.75)
  par(mar = c(10,3,3,7))  #sets the bottom, left, top and right margins respectively of the plot region in number of lines of text.
  par(oma = c(7,0,0.5,5)) #outer margins (unity = line)
  par(xpd=NA) #allow to plot outside the plot
         heatmap.2(as.matrix(scale(dataPCA.temp,center = TRUE, scale=TRUE)),
            main=paste0("Heat map for: ",paste(liststates[[i]], collapse=", "),#"\n",
                      #" (",statedescription[i],")",
                      "\n","method=",hc.temp$method), cex.main=1.0,
            srtCol = 30, #rotate x-label
            dendrogram = "row", #dendrogram only in the lines
            Rowv = dend.temp, #dendrogram to be chosen
            Colv = "NA", # to make sure the columns are not ordered
            trace="none",          
            margins =c(5,0.1),      
            key.xlab = "stand. value",
            #denscol = "grey",
            density.info = "none",  #desactivated color by densitiy of distribution
            labRow = datainputNA.temp$ID.PCA, #labes of the lines, heatmap will already reorder
            colRow = "black",
            #cellnote = as.matrix(scale(dataPCA.temp,center = TRUE, scale=TRUE), ##legenda de valor em cada quadrado do plot
            #RowSideColors = rainbow_hcl(length(unique(datainputNA.temp$Categ)))[(as.numeric(categfactor.temp))[order.dendrogram(dend.temp)]], # to add color strips related to the dendrogram label
            col = colorRampPalette(c("red", "white", "blue"))(n = 299), #,
            #breaks=col_breaks    # enable color transition at specified limits
            RowSideColors = rainbow_hcl(length(unique(datainputNA.temp$State)))[as.numeric(factor(datainputNA.temp$State))]) # to add color strips related to the dendrogram label
          legend(x=0.3,y=-0.55, title="Dendrogram colors",
                legend = paste0(1:nbg.temp), 
                fill=rainbow_hcl(nbg.temp), cex=0.7,xpd=NA,
                ncol=if (length(nbg.temp)>5){2}else{1})
         legend(x=0.025,y=-0.55, title="First column colors",
                legend = levels(factor(datainputNA.temp$State)), 
                fill=rainbow_hcl(length(unique(datainputNA.temp$State))), cex=0.7,xpd=NA,
                ncol=if (length(nbg.temp)>5){2}else{1})
            dev.off()
            
  } 

#### 4.2.2. Tools to define number of clusters ####
for (i in 1:(length(states)))
{
temp=states[i]
dataPCA.complete.temp=get(paste0("dataPCA.complete.",temp))
scaled_data=as.matrix(scale(dataPCA.complete.temp))
hc.temp=get(paste0("hc.",temp))

#plots
pdf(paste0(wd,"/",temp,"/Optimal_n_cluster/",toupper(temp),"_HCA_n-cluster_wss.",hc.temp$hclust$method,".pdf"),
    width = 7, height = 8.75)
  plot(fviz_nbclust(scaled_data, FUN = hcut, method = "wss", k.max = 15))
dev.off()
  
pdf(paste0(wd,"/",temp,"/Optimal_n_cluster/",toupper(temp),"_HCA_n-cluster_silh.",hc.temp$hclust$method,".pdf"),
    width = 7, height = 8.75)
  plot(fviz_nbclust(scaled_data, FUN = hcut, method = "silhouette"))
dev.off()
  
gap_stat = clusGap(scaled_data, FUN = hcut, K.max = 10, B = 1000)
pdf(paste0(wd,"/",temp,"/Optimal_n_cluster/",toupper(temp),"_HCA_n-cluster_gapstat.",hc.temp$hclust$method,".pdf"),
    width = 7, height = 8.75)
 plot(fviz_gap_stat(gap_stat))
dev.off()

}

#### 4.2.3. Defining final clusters (typology) ####
############################################## STOP HERE #################################################################################################

#### 4.2.4 Manual definition of cluster number and name ####
###### CLUSTERS DEFINED FOR FERT-VALUE TYPOLOGY IN THE PAPER ####
hcagnames.all=c("") #not being used
hcagnames.description.all=c("") #not being used

# AS DEFINED ON PAPER: #
hcagnames.raw=c("1. Fibrous feedstock: Cattle slurry, silage mono/co-dig.(n=19/22)",#1    1. Lab. OFMSW/GW",
                "2. Sewage sludge, Biowaste, FAI mono/co-digestion (n=18/19)",#"2. Manure and FAI co-digestion",
                "3. OFMSW, Food Waste, FAI, Pig slurry mono/co-digestion (n=23/23)",
                "4. LAB. PERCOLATE SYSTEM: Food Waste / Green Waste co-dig. (n=4/4)",
                "5. Manure/other co-digestion (n=4/4)",
                "6. OFMSW and Biowaste mono/co-digestion (n=6/8)",
                "7. Fibrous feedstock: Cattle manure, green waste, silage (n=11/11)"
                )
hcagnames.description.raw=c() #not being used
  
hcagnames.liq=c("1. Fibrous material (n=14/16), low perf. separation (n=11/16, NA = 2/16)",
                "2. Pig slurry, FAI, Biowaste, OFMSW mono/co-digestion (n=9/10), high perf. separation (n=6/10, NA = 2/10)")
                
hcagnames.description.liq=c() #not being used

hcagnames.sld=c("1. SS mono/co-dig. (n=12/18, pig slurry mono/co-dig. (n=5/18), high perf. separation (n=13/18, NA = 4/18)",
                "2. Biowaste/Fat: high TAN (n=1/1)",
                "3. Fibrous material (n=15/17), low perf. separation (n=13/16)"
                )
                
hcagnames.description.sld=c() #not being used



#### 4.2.5 Creating cluster data frames ####
#counting the groups
for (i in 1:length(states))
{  temp=states[i]
   nbg.temp=length(get(paste0("hcagnames.",temp)))
   assign(paste0("nbg.",temp),nbg.temp)}

# Creating groups' data frames
for (i in 1:length(states))
{
  print(paste0("running i=",i," /",length(states)))
  temp=states[i]
datainputNA.temp=get(paste0("datainputNA.",temp))
hc.temp=get(paste0("hc.",temp))
nbg.temp=get(paste0("nbg.",temp))
if (nbg.temp == 0) {nbg.temp=1}
hcagnames.temp=get(paste0("hcagnames.",temp))

HCAgroups.temp=cutree(hc.temp,nbg.temp)  ###cut tree counting is maybe different from dendogram!!!
HCAgroups.temp=data.frame(datainputNA.temp$ID,HCAgroups.temp,0)
colnames(HCAgroups.temp)=c("ID","HCAgroups","group.name")
assign(paste0("HCAgroups.",temp),HCAgroups.temp)
}

# Adjusting cutree order and matching group names
for (i in 1:length(states))
{
  print(paste0("running i=",i," /",length(states)))
  temp=states[i]
  hc.temp=get(paste0("hc.",temp))
  HCAgroups.temp=get(paste0("HCAgroups.",temp))
  cutree.order.temp=unique(HCAgroups.temp$HCAgroups[hc.temp$order])
  for (k in 1:length(HCAgroups.temp$HCAgroups)){
    HCAgroups.temp$HCAgroups[k]=which(cutree.order.temp==HCAgroups.temp$HCAgroups[k])
  assign(paste0("HCAgroups.",temp),HCAgroups.temp)
  }
  
  # changing number for names  
  nbg.temp=get(paste0("nbg.",temp))
  hcagnames.temp=get(paste0("hcagnames.",temp))
  for (j in 1:nbg.temp){  
    HCAgroups.temp$group.name[HCAgroups.temp$HCAgroups==j]=hcagnames.temp[j]
  }
  assign(paste0("HCAgroups.",temp),HCAgroups.temp)
  
  ## exporting to excel ##
  HCAsummary.temp=data.frame("method"=hc.temp$method, "dist.method"=hc.temp$dist.method)
  assign(paste0("HCAsummary",temp),HCAsummary.temp)
  
  hcaheight.temp=as.data.frame(hc.temp$height)
  hcaorder.temp=as.data.frame(hc.temp$order)
  
  writehcadata=c("HCAgroups.temp","HCAsummary.temp","hcaheight.temp","hcaorder.temp")
  
  WriteXLS(writehcadata,  
           ExcelFileName = paste0(wd,"/",temp,"/HCA/",toupper(temp),"_HCA_output_data.xls"), 
           SheetNames = c("Groups","Summary","Height","Order"), perl = "perl",
           verbose = FALSE, Encoding = c("UTF-8", "latin1", "cp1252"),
           row.names = FALSE, col.names = TRUE,
           AdjWidth = FALSE, AutoFilter = TRUE, BoldHeaderRow = TRUE,
           na = "",
           FreezeRow = 1, FreezeCol = 0,
           envir = parent.frame())
}

# Creating final clusters
for (i in 1:length(states))
{
  print(paste0("running i=",i," /",length(states)))
  temp=states[i]
  dataPCA.temp=get(paste0("dataPCA.",temp))
  hc.temp=get(paste0("hc.",temp))
  nbg.temp=get(paste0("nbg.",temp))
  if (nbg.temp == 0) {nbg.temp=1}
  hcagnames.temp=get(paste0("hcagnames.",temp))
  hcagnames.description.temp=get(paste0("hcagnames.description.",temp))
  datainputNA.temp=get(paste0("datainputNA.",temp))
  
  dend.temp = as.dendrogram(hc.temp)
  dend.temp= color_branches(dend.temp, k=nbg.temp) 
  #Definition of the labels#
  labels(dend.temp)=paste0(datainputNA.temp$State[order.dendrogram(dend.temp)],"-",datainputNA.temp$Categ[order.dendrogram(dend.temp)],"(",datainputNA.temp$ID[order.dendrogram(dend.temp)],")")
  dend.temp = hang.dendrogram(dend.temp,hang_height=0.1)
  # dend.temp <- assign_values_to_leaves_nodePar(dend, 0.5, "lab.cex")  ##IF NEEDE TO REDUCE SIZE OF LABELS
  # size of the ind labels
  dend.temp = set(dend.temp, "labels_cex", (0.15+10.5/length(datainputNA.temp$ID.PCA)))
  
  
  pdf(paste0(wd,"/",temp,"/HCA/",toupper(temp),"_HCA_meth.",hc.temp$method,"_FINAL.pdf"),
      width = 7, height = 8.75)
  par(mar = c(3,3,3,7))  
      plot(dend.temp, 
         main=paste0("Clustered data set for: ",paste(liststates[[i]], collapse=", "),"\n"," (",statedescription[i],")") ,
         cex.main=1.0, #cex.sub=1,
       horiz =  TRUE,  nodePar = list(cex = 0.07))
      legend("topleft",legend = paste0(hcagnames.temp,": ",hcagnames.description.temp), fill=rainbow_hcl(nbg.temp), cex=0.55)
  dev.off()
  
 
  pdf(paste0(wd,"/",temp,"/HCA/",toupper(temp),"_HCA_HEATMAP.",hc.temp$method,"_FINAL.pdf"),
      width = 7, height = 8.75)
  par(mar = c(10,3,3,7))  #sets the bottom, left, top and right margins respectively of the plot region in number of lines of text.
  par(oma = c(7,0,0.5,5)) #outer margins (unity = line)
  par(xpd=NA) #allow to plot outside the plot
  heatmap.2(as.matrix(scale(dataPCA.temp,center = TRUE, scale=TRUE)),
            main=paste0("Heat map for: ",paste(liststates[[i]], collapse=", "),#"\n",
                        #" (",statedescription[i],")",
                        "\n","method=",hc.temp$method), cex.main=1.0,
            srtCol = 30, #rotate x-label
            dendrogram = "row", #dendrogram only in the lines
            Rowv = dend.temp, #dendrogram to be chosen
            Colv = "NA", # this to make sure the columns are not ordered
            trace="none",          
            margins =c(5,0.1),      
            key.xlab = "stand. value",
            keysize = 1.0,
            #denscol = "grey",
            density.info = "none",  #desactivated color by densitiy of distribution
            labRow = if(length(hc.temp$order)>120){FALSE} else {paste0(datainputNA.temp$State,"(",datainputNA.temp$ID,")")}, #labes of the lines
            cexRow = if(length(hc.temp$order)>120){NULL} else {1/log10(length(hc.temp$order))},
            colRow = if(length(hc.temp$order)>120){NULL} else {"black"},
            offsetRow = 0.2,
            #cellnote = as.matrix(scale(dataPCA.temp,center = TRUE, scale=TRUE), ##legenda de valor em cada quadrado do plot
            #RowSideColors = rainbow_hcl(length(unique(datainputNA.temp$Categ)))[(as.numeric(categfactor.temp))[order.dendrogram(dend.temp)]], # to add color strips related to the dendrogram label
            RowSideColors = if (length(liststates[[i]]) == 1) {rep("white",length(datainputNA.temp$State))} else {rainbow_hcl(length(unique(datainputNA.temp$State)))[as.numeric(factor(datainputNA.temp$State))]}, # to add color strips related to the dendrogram label
            col = colorRampPalette(c("red", "white", "blue"))(n = 299))
  #breaks=col_breaks    # enable color transition at specified limits
  legend(x=0.3,y=-0.53, title="Dendrogram colors",
         legend = paste0(hcagnames.temp), 
         fill=rainbow_hcl(nbg.temp), cex=0.7,xpd=NA,
         ncol=if (length(hcagnames.temp)>5){2}else{1})
  if (length(liststates[[i]]) > 1) { legend(x=0.025,y=-0.53, title="First column colors",
         legend = levels(factor(datainputNA.temp$State)), 
         fill=rainbow_hcl(length(unique(datainputNA.temp$State))), cex=0.7,xpd=NA,
         ncol=if (length(hcagnames.temp)>5){2}else{1})}
    dev.off()
  }

#### 4.3 PCA ####
#### 4.3.1. Creating PCA objects ####
for (i in 1:length(states))
{
  print(paste0("running i=",i," /",length(states)))
  temp=states[i]
  datainputNA.temp=get(paste0("datainputNA.",temp))
  dataPCA.temp=get(paste0("dataPCA.",temp))
  HCAgroups.temp=get(paste0("HCAgroups.",temp))
  dataPCA.complete.temp = get(paste0("dataPCA.complete.",temp))
  nb.temp = get(paste0("nb.",temp))

if (nb.temp$ncp==1) {
  nb.temp$ncp=2
}

#### RUNNING PCA MODEL
NAtable.temp=get(paste0("NAtable.",temp))
if (sum(NAtable.temp["NA.count",])==0) {

  digest.PCA.temp <- PCA(dataPCA.temp, scale.unit = TRUE,ncp=nb.temp$ncp,graph = FALSE)
    assign(paste0("digest.PCA.",temp),digest.PCA.temp)

} else {

  digest.PCA.temp <- PCA(dataPCA.complete.temp, scale.unit = TRUE,ncp=nb.temp$ncp,graph = FALSE)
    assign(paste0("digest.PCA.",temp),digest.PCA.temp)
}
} 

#### 4.4.2. PCA excel outputs: eigenvalues, correlations, loading matrix ####
for (i in 1:length(states))
{
  print(paste0("running i=",i," /",length(states)))
  temp=states[i]
  digest.PCA.temp=get(paste0("digest.PCA.",temp))
  nb.temp = get(paste0("nb.",temp))

  
## Excel data ##
# eigen values
eigenvalues.temp <- digest.PCA.temp$eig
  assign(paste0("eigenvalues.",temp),as.data.frame(eigenvalues.temp))

# loading matrix
loading.matrix.temp=digest.PCA.temp$var$coord
  assign(paste0("loading.matrix.",temp),as.data.frame(loading.matrix.temp))

# Correlations between PCs and input variables
PCcorr.temp=dimdesc(digest.PCA.temp, axes=1:(nb.temp$ncp))
  assign(paste0("PCcorr.",temp),PCcorr.temp)
PC1corr.temp=PCcorr.temp$Dim.1$quanti
  assign(paste0("PC1corr.",temp),as.data.frame(PC1corr.temp))
PC2corr.temp=PCcorr.temp$Dim.2$quanti
  assign(paste0("PC2corr.",temp),as.data.frame(PC2corr.temp))
if (nb.temp$ncp>2) {
PC3corr.temp=PCcorr.temp$Dim.3$quanti
assign(paste0("PC3corr.",temp),as.data.frame(PC3corr.temp))  
}
# aka values
aka.temp=digest.PCA.temp$ind$coord
  assign(paste0("aka.",temp),as.data.frame(aka.temp))
  
# EXCEL OUTPUT
write.PCA.data=c("source.summary",paste0("eigenvalues.",temp),paste0("loading.matrix.",temp),
                 paste0("PC1corr.",temp),paste0("PC2corr.",temp),if(nb.temp$ncp>2){paste0("PC3corr.",temp)})
sheetnames.PCA = c("Data.summary","Eigen_values","load_matrix","PC1_corr","PC2_corr",if(nb.temp$ncp>2){"PC3_corr"})
WriteXLS(write.PCA.data, 
         ExcelFileName = paste(wd,"/",temp,"/PCA/",toupper(temp),"_PCA_output_data.xls",sep=""), 
         SheetNames = sheetnames.PCA,
         perl = "perl",
         verbose = FALSE, Encoding = c("UTF-8", "latin1", "cp1252"),
         row.names = TRUE, col.names = TRUE,
         AdjWidth = FALSE, AutoFilter = TRUE, BoldHeaderRow = TRUE,
         na = "",
         FreezeRow = 1, FreezeCol = 0,
         envir = parent.frame())
} 

#### 4.4.3. PCA plots ####

for (i in 1:length(states))
{
  print(paste0("running i=",i," /",length(states)))
  temp=states[i]
  digest.PCA.temp=get(paste0("digest.PCA.",temp))
  datainputNA.temp=get(paste0("datainputNA.",temp))
  
# eigen values #
pdf(paste0(wd,"/",temp,"/PCA/",toupper(temp),"_PCA_eigen.pdf"))
  print(fviz_screeplot(digest.PCA.temp, ncp=10))
  dev.off()

### PCA plots 
  
# VARIABLES ONLY: contribution of variables

pdf(paste0(wd,"/",temp,"/PCA/",toupper(temp),"_PCA_var_contrib.pdf"))
  print(fviz_pca_var(digest.PCA.temp, col.var="contrib")+
  scale_color_gradient2(low="white", mid="blue", 
                        high="red", midpoint=55)+theme_bw())
  dev.off()

# INDIVIDUALS ONLY: plot by states (Raw/SF/LF,...) 
if (length(liststates[[i]])>1){
pdf(paste0(wd,"/",temp,"/PCA/",toupper(temp),"_PCA_ind_group_state.pdf"))
  print(fviz_pca_ind(digest.PCA.temp, label="none", habillage=as.factor(datainputNA.temp$State),
             addEllipses=TRUE, ellipse.level=0.95))
  dev.off()
}

#### BIPLOTS: VAR+IND 
## COS2+state 
if (length(liststates[[i]])>1) {
  pdf(paste0(wd,"/",temp,"/PCA/",toupper(temp),"_PCA_biplot_group_stateCOS2.pdf"))
  print(fviz_pca_biplot(digest.PCA.temp, habillage = as.factor(datainputNA.temp$State), addEllipses = FALSE,
                  col.var = "red", alpha.var ="cos2",
                  label = "var") +
    scale_color_brewer(palette="Dark2"))
    dev.off()
       
  pdf(paste0(wd,"/",temp,"/PCA/",toupper(temp),"_PCAbiplot.group_state.pdf"))
  print(fviz_pca_biplot(digest.PCA.temp, habillage = as.factor(datainputNA.temp$State), label="var", addEllipses = FALSE))
    dev.off()
} else 
    {
  pdf(paste0(wd,"/",temp,"/PCA/",toupper(temp),"_PCA_biplot_COS2.pdf"))
  print(fviz_pca_biplot(digest.PCA.temp, addEllipses = FALSE,
                  col.var = "red", alpha.var ="cos2",
                  label = "var") +
    scale_color_brewer(palette="Dark2"))
    dev.off()
  
  pdf(paste0(wd,"/",temp,"/PCA/",toupper(temp),"_PCAbiplot.pdf"))
    print(fviz_pca_biplot(digest.PCA.temp, label="var", addEllipses = TRUE))
  dev.off()
}


if (length(datainputNA.temp$Sep.meth)!=sum(is.na(datainputNA.temp$Sep.meth))) {
pdf(paste0(wd,"/",temp,"/PCA/",toupper(temp),"_PCAbiplot_label_sep.method.pdf"))
plot(fviz_pca_biplot(digest.PCA.temp, habillage = as.factor(datainputNA.temp$State), label="var", addEllipses = FALSE) 
  + annotate(geom="text", x=digest.PCA.temp$ind$coord[,1], y=digest.PCA.temp$ind$coord[,2], 
           label=datainputNA.temp$Sep.meth, color="black")
  + theme(legend.title=element_text(colour="black", size=12, face="bold"),
            legend.position="bottom",legend.direction="vertical")
 )
dev.off()
}

## BIplot label=ID
pdf(paste0(wd,"/",temp,"/PCA/",toupper(temp),"_PCAbiplot_label_ID.pdf"))
print(fviz_pca_biplot(digest.PCA.temp, addEllipses = FALSE, col.var = "contrib", label = "var") +
  scale_color_gradient2(low="white", mid="blue", high="red", midpoint=55)+
  geom_text_repel(aes(x=digest.PCA.temp$ind$coord[,1], 
                      y=digest.PCA.temp$ind$coord[,2],label=datainputNA.temp$ID)))
dev.off()

# Biplot label = feed categ  #
pdf(paste0(wd,"/",temp,"/PCA/",toupper(temp),"_PCAbiplot_label_feed.categ.pdf"))
print(fviz_pca_biplot(digest.PCA.temp, addEllipses = FALSE, col.var = "contrib", label = "var") +
  scale_color_gradient2(low="white", mid="blue", high="red", midpoint=55)+
  geom_text_repel(aes(x=digest.PCA.temp$ind$coord[,1], 
                      y=digest.PCA.temp$ind$coord[,2],
                      label=datainputNA.temp$Feed.description)))
dev.off()

# Plot by HCA groups 
  HCAgroups.temp=get(paste0("HCAgroups.",temp))
  nbg.temp=get(paste0("nbg.",temp))
  shapevaluesb=c(15,16,17,18,25,7,8,9,10,11,12,13,4,3,14,0,1,2,3,4,5,6,21,22,23,24,25)
# WITH ELLIPSE  
  if (nbg.temp>1) {
    pdf(paste0(wd,"/",temp,"/PCA/",toupper(temp),"_PCAbiplot12.HCAgroups.pdf"))
    plot(fviz_pca_biplot(digest.PCA.temp, 
                         habillage = as.factor(HCAgroups.temp$group.name), addEllipses = TRUE, label="var") ########## if there is group names
         + scale_shape_manual(values=shapevaluesb[1:length(unique(HCAgroups.temp$group.name))])
         + guides(fill=guide_legend(ncol=2))
         + theme(legend.title=element_text(colour="black", size=12, face="bold"),
                 legend.position="bottom",legend.direction="vertical")
    )
    dev.off()
    if (digest.PCA.temp$call$ncp > 2) {
    pdf(paste0(wd,"/",temp,"/PCA/",toupper(temp),"_PCAbiplot13.HCAgroups.pdf"))
    plot(fviz_pca_biplot(digest.PCA.temp, axes = c(1,3),
                         habillage = as.factor(HCAgroups.temp$group.name), addEllipses = TRUE, label="var") ########## if there is group names
         + scale_shape_manual(values=shapevaluesb[1:length(unique(HCAgroups.temp$group.name))])
         + guides(fill=guide_legend(ncol=2))
         + theme(legend.title=element_text(colour="black", size=12, face="bold"),
                 legend.position="bottom",legend.direction="vertical")
    )
    dev.off()
    }
  }
# WITHOUT ELLIPSE
if (nbg.temp>1) {
  pdf(paste0(wd,"/",temp,"/PCA/",toupper(temp),"_PCAbiplot12.HCAgroups NO ELLIPSE.pdf"))
  plot(fviz_pca_biplot(digest.PCA.temp, 
                        habillage = as.factor(HCAgroups.temp$group.name), addEllipses = FALSE, label="var") ########## if there is group names
        +  scale_shape_manual(values=shapevaluesb[1:length(unique(HCAgroups.temp$group.name))])
        + scale_fill_manual(name="HCA groups")
       + guides(fill=guide_legend(ncol=2))
        + theme(legend.title=element_text(colour="black", size=12, face="bold"),
        legend.position="bottom",legend.direction="vertical")
               )
  if (digest.PCA.temp$call$ncp > 2) {
    dev.off()
    pdf(paste0(wd,"/",temp,"/PCA/",toupper(temp),"_PCAbiplot13.HCAgroups NO ELLIPSE.pdf"))
    plot(fviz_pca_biplot(digest.PCA.temp, axes = c(1,3),
                         habillage = as.factor(HCAgroups.temp$group.name), addEllipses = FALSE, label="var") ########## if there is group names
         + scale_shape_manual(values=shapevaluesb[1:length(unique(HCAgroups.temp$group.name))])
         + scale_fill_manual(name="HCA groups")
         + guides(fill=guide_legend(ncol=2))
         + theme(legend.title=element_text(colour="black", size=12, face="bold"),
                 legend.position="bottom",legend.direction="vertical")
    )
    dev.off()
  }
    }

} 

#### 5. Boxplots with absolute values ####
#### 5.1. Selecting retained HCA groups
#### 5.1.1. As used for the fertilizing value typology ####
Final.groups=list(raw=c(1,2,3,5,6,7),
                  liq=c(1,2),
                  sld=c(1,3))

#### 5.2.Ploting for the original HCA/PCA variables (but not center-scaled) ####
for (i in 1:length(Final.groups)){
  print(paste0("running i=",i," /",length(Final.groups)))
  temp=names(Final.groups)[i]
  HCAgroups.temp=get(paste0("HCAgroups.",temp))
  datainputNA.temp=get(paste0("datainputNA.",temp))
  dataPCA.temp=get(paste0("dataPCA.",temp))
  datainputBP.temp=data.frame(HCAgroups.temp,group=0,dataPCA.temp)
  datainputBP.temp$group=paste0(datainputBP.temp$HCAgroups)
  datainputBP.temp=subset.data.frame(datainputBP.temp, datainputBP.temp$HCAgroups %in% Final.groups[[i]])
  datainputBP_fert_potential.temp = ddply(datainputBP.temp,
                                          .(ID,HCAgroups,group.name,group),
                                          summarise, 
                                          Fert_pot = 0.1*sum(TN,TP*2.2915,TK*1.2047))
  datainputBP_fert_potential.temp=na.omit(datainputBP_fert_potential.temp)
  assign(paste0("datainputBP.",temp),datainputBP.temp)
  datainputBP_fert_potential.temp = data.frame(datainputBP_fert_potential.temp,state=temp)
  assign(paste0("datainputBP_fert_potential.",temp),datainputBP_fert_potential.temp)
    dimelt.temp=melt(datainputBP.temp, id.vars=c("ID","HCAgroups","group.name","group"))

### PLOT
  
  plot=ggplot(data=dimelt.temp, aes(x=group, y=value))+geom_boxplot(aes(fill=group.name))
  plot_fertpot=ggplot(data=datainputBP_fert_potential.temp, aes(x=group, y=Fert_pot))+geom_boxplot(aes(fill=group.name))
  
  give.n <- function(x){
    return(c(y = 1.02*as.numeric(quantile(x,prob=0.75)), label = length(x)))}
  mean.n <- function(x){
    return(c(y = mean(x), label = round(mean(x),2)))} 
  
  pdf(paste0(wd,"/",temp,"/Boxplots/",toupper(temp),"_singleplot.pdf"), width = 7, height = 4)
  plot(plot
       + facet_wrap(~variable ,scales="free",nrow = 1, ncol = NULL) #,labeller=labeller(.cols=labelsVAR))
       + labs(title="Boxplots for HCA groups",x="HCA groups",y="Values in different scale")
       + theme(plot.title = element_text(size = rel(1.5), colour = "black"),
               legend.title=element_blank(),
               legend.position="bottom",
               legend.direction="horizontal")
       + guides(fill = guide_legend(nrow = 2))
       + stat_summary(fun.data = give.n, geom = "text", aes(shape="Nb. of ind."), fun.y = median,show.legend=FALSE, size=3 )
       + stat_summary(fun.data = mean.n, aes(shape="Mean"), colour = "black", geom="point")
       + scale_shape_manual("", values=c("Mean"="x","Nb. of ind."="N"))
  )
  dev.off()
  
  pdf(paste0(wd,"/",temp,"/Boxplots/",toupper(temp),"_FERT_POTENTIAL_singleplot.pdf"), width = 7, height = 4)
  plot(plot_fertpot
       + labs(title="Boxplots for HCA groups",x="HCA groups",y="N + P2O5 + K2O (%DM)")
       + theme(plot.title = element_text(size = rel(1.5), colour = "black"),
               legend.title=element_blank(),
               legend.position="bottom",
               legend.direction="horizontal")
       + geom_hline(yintercept = 7,color="red")
       + guides(fill = guide_legend(nrow = 2))
       + stat_summary(fun.data = give.n, geom = "text", aes(shape="Nb. of ind."), fun.y = median,show.legend=FALSE, size=3 )
       + stat_summary(fun.data = mean.n, aes(shape="Mean"), colour = "black", geom="point")
       + scale_shape_manual("", values=c("Mean"="x","Nb. of ind."="N"))
  )
  dev.off()
  
}
#### Not in paper: PLOTTING FERTILIZER POTENTIAL (DM NUTRIENT CONTENT AS MINERAL EQUIVALENTS) FOR EVERY STATE IN THE SAME PLOT ####
data_fert_potent_all = rbind(datainputBP_fert_potential.raw,datainputBP_fert_potential.sld,datainputBP_fert_potential.liq)
plot_fertpot_all=ggplot(data=data_fert_potent_all, aes(x=group, y=Fert_pot))+geom_boxplot(aes(fill=group.name))

pdf(paste0(wd,"/FERT_POTENTIAL_singleplot.pdf"), width = 7, height = 4)
plot(plot_fertpot_all
     + facet_wrap(~state ,scales="free",nrow = 1, ncol = NULL) #,labeller=labeller(.cols=labelsVAR))
     + labs(title="Boxplots for HCA groups",x="HCA groups",y="N + P2O5 + K2O (%DM)")
     + theme(plot.title = element_text(size = rel(1.5), colour = "black"),
             legend.title=element_blank(),
             legend.position="bottom",
             legend.direction="horizontal")
     + geom_hline(yintercept = 7,color = "red")
     + guides(fill = guide_legend(nrow = 2))
     + stat_summary(fun.data = give.n, geom = "text", aes(shape="Nb. of ind."), fun.y = median,show.legend=FALSE, size=3 )
     + stat_summary(fun.data = mean.n, aes(shape="Mean"), colour = "black", geom="point")
     + scale_shape_manual("", values=c("Mean"="x","Nb. of ind."="N"))
)
dev.off()

pdf(paste0(wd,"/FERT_POTENTIAL_singleplot_scalefixed.pdf"), width = 7, height = 4)
plot(plot_fertpot_all
     + facet_wrap(~state ,scales="free_x",nrow = 1, ncol = NULL) #,labeller=labeller(.cols=labelsVAR))
     + labs(title="Boxplots for HCA groups",x="HCA groups",y="N + P2O5 + K2O (%DM)")
     + theme(plot.title = element_text(size = rel(1.5), colour = "black"),
             legend.title=element_blank(),
             legend.position="bottom",
             legend.direction="horizontal")
     + geom_hline(yintercept = 7,color = "red")
     + guides(fill = guide_legend(nrow = 2))
     + stat_summary(fun.data = give.n, geom = "text", aes(shape="Nb. of ind."), fun.y = median,show.legend=FALSE, size=3 )
     + stat_summary(fun.data = mean.n, aes(shape="Mean"), colour = "black", geom="point")
     + scale_shape_manual("", values=c("Mean"="x","Nb. of ind."="N"))
)
dev.off()



#### 5.3. Ploting for new variables within "other var." input data ####
#### OBS. COLUMNS MANUALLY ADJUSTED ####
for (i in 1:length(Final.groups)){
  print(paste0("running i=",i," /",length(Final.groups)))
  temp=names(Final.groups)[i]
  HCAgroups.temp=get(paste0("HCAgroups.",temp))
  datainput.othervar.temp=datainput.othervar[which(datainput.othervar$ID %in% HCAgroups.temp$ID),]
  dataVAR.temp=datainput.othervar.temp[,9:ncol(datainput.othervar)]
  variables=c("C/N",
              "C/Norg",
              "DM (%)",
              "VS (%FM)",
              "VS (%DM)",
              "TAN (g.kgFM-1)",
              "TAN/TN (%)",
              "Norg (g.kgFM-1)",
              "TN (g.kgFM-1)",
              "P (g.kgFM-1)",
              "K (g.kgFM-1)",
              "S (g.kgFM-1)",
              "Ca (g.kgFM-1)",
              "Mg (g.kgFM-1)",
              "P2O5 (g.kgFM-1)",
              "K2O (g.kgFM-1)",
              "SO3 (g.kgFM-1)",
              "CaO (g.kgFM-1)",
              "MgO (g.kgFM-1)",
              "TN+P2O5+K2O",
              "Max(TN,P2O5,K2O)",
              "As (g.kgDM-1)",
              "Cd (g.kgDM-1)",
              "Cr (g.kgDM-1)",
              "Cu (g.kgDM-1)",
              "Hg (g.kgDM-1)",
              "Ni (g.kgDM-1)",
              "Pb (g.kgDM-1)",
              "Se (g.kgDM-1)",
              "Zn (g.kgDM-1)",
              "Res. Biogas (NL.gVS-1)"
  )
  colnames(dataVAR.temp)=variables
  datainput.othervar.temp=data.frame(HCAgroups.temp,group=0,dataVAR.temp)
  datainput.othervar.temp$group=paste0("",datainput.othervar.temp$HCAgroups) # tinha feito isso pra add um G antes do numero do grupo
  datainput.othervar.temp=subset.data.frame(datainput.othervar.temp, datainput.othervar.temp$HCAgroups %in% Final.groups[[i]])
  
  
  databasic.temp=datainput.othervar.temp[,c(1:4,5:13,19,20,22:25,35)]
  heavymetals.temp=datainput.othervar.temp[,c(1:4,26:34)]
  print(colnames(databasic.temp))
  colnames(databasic.temp)[5:20]=variables[-4+c(5:13,19,20,22:25,35)]
  print(colnames(databasic.temp))
  print(colnames(heavymetals.temp))
  colnames(heavymetals.temp)[5:13]=variables[-4+c(26:34)]
  print(colnames(heavymetals.temp))
  
  # Removing variablies with only NA for all HCA groups 
  remove1 = as.numeric(which(lapply(databasic.temp[,5:ncol(databasic.temp)],mean,na.rm = T) == "NaN"))
  remove2 = as.numeric(which(lapply(heavymetals.temp[,5:ncol(heavymetals.temp)],mean,na.rm = T) == "NaN"))
  if (sum(remove1) != 0) {
    databasic.temp=databasic.temp[,-4-remove1]
  }
  if (sum(remove2) != 0) {
    heavymetals.temp=heavymetals.temp[,-4-remove2]
  }
  
  
  dimelt.databasic.temp=melt(databasic.temp, id.vars=c("ID","HCAgroups","group.name","group"))
  dimelt.heavymetals.temp=melt(heavymetals.temp, id.vars=c("ID","HCAgroups","group.name","group"))
  
  # PLOT
  
  give.n <- function(x){
    return(c(y = 1.02*as.numeric(quantile(x,prob=0.75)), label = length(x)))}
  mean.n <- function(x){
    return(c(y = mean(x), label = round(mean(x),2)))} 
  plot2=ggplot(data=dimelt.databasic.temp, aes(x=group, y=value))+geom_boxplot(aes(fill=group.name))
  plot4=ggplot(data=dimelt.heavymetals.temp, aes(x=group, y=value))+geom_boxplot(aes(fill=group.name))
  
  pdf(paste0(wd,"/",temp,"/Boxplots_LEGISLATION/",toupper(temp),"_legislation_basics_singleplot.pdf"))
  plot(plot2
       + facet_wrap(~variable ,scales="free",nrow = NULL, ncol = NULL) #,labeller=labeller(.cols=labelsVAR))
       + labs(title="Boxplots for HCA groups",x="HCA groups",y="Values in different scale")
       + theme(plot.title = element_text(size = rel(1.5), colour = "black"),
               legend.title=element_blank(),
               legend.position="bottom",
               legend.direction="horizontal")
       + guides(fill = guide_legend(nrow = 2))
       + stat_summary(fun.data = give.n, geom = "text", aes(shape="Nb. of ind."), fun.y = median,show.legend=FALSE, size=3 )
       + stat_summary(fun.data = mean.n, aes(shape="Mean"), colour = "black", geom="point")
       + scale_shape_manual("", values=c("Mean"="x","Nb. of ind."="N"))
  )
  dev.off()
  
  pdf(paste0(wd,"/",temp,"/Boxplots_LEGISLATION/",toupper(temp),"_heavy_metals_singleplot.pdf"))
  plot(plot4
       + facet_wrap(~variable ,scales="free",nrow = NULL, ncol = NULL) #,labeller=labeller(.cols=labelsVAR))
       + labs(title="Boxplots for HCA groups",x="HCA groups",y="Values in different scale")
       + theme(plot.title = element_text(size = rel(1.5), colour = "black"),
               legend.title=element_blank(),
               legend.position="bottom",
               legend.direction="horizontal")
       + guides(fill = guide_legend(nrow = 2))
       + stat_summary(fun.data = give.n, geom = "text", aes(shape="Nb. of ind."), fun.y = median,show.legend=FALSE, size=3 )
       + stat_summary(fun.data = mean.n, aes(shape="Mean"), colour = "black", geom="point")
       + scale_shape_manual("", values=c("Mean"="x","Nb. of ind."="N"))
  )
  dev.off()
  
}

