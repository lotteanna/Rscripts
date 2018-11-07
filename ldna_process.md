Evaluating various inputs for ldna
===

Script by Lotte van Boheemen 2018

```{r}
library(dplyr)
library(plyr)
library(tidyr)
library(gdata)
library(igraph)
#devtools::install_github("petrikemppainen/LDna", ref = 'v.64')
library(LDna)
```

Run analyses
```{r}
####EDIT###

#infiles
LDmat = N1
ldfile = LDnaRaw(LDmat)
tmpfile = "na50_aw_phi"

#range of values to test
edge = seq(1,20,by=1)
phis = seq(1,6, by=1)

#Set criteria for post-test summary filtering
Merge = .2
SNPS = 10


###DO NOT EDIT#####
Dir <- paste0("LD/", tmpfile,"/")
sums = "summs"
lis = "lists"

if (file.exists(Dir)){
   print("dir exists")
} else {
    dir.create(file.path(Dir))
}

setwd(Dir)

if (file.exists(sums)){
   print("summs exists")
} else {
    dir.create(file.path(sums))
}
if (file.exists(lis)){
   print("lists exists")
} else {
    dir.create(file.path(lis))
}

pdf(paste0("ld_",tmpfile,".pdf"))
tryCatch(for(i in edge){
  tryCatch(for(j in phis){
    tmp = extractClusters(ldfile, min.edges = i, plot.tree = TRUE, extract=TRUE, rm.COCs = T, phi=j)
    save(tmp, file=paste0(tmpfile,"_clus_",i,"_",j))
    sink(paste0("lists/",tmpfile,"_clus_",i,"_",j,".txt")); print(tmp$clusters);sink() 
    }, error = function(e) {an.error.occured <<- TRUE})
}, error = function(e) {an.error.occured <<- TRUE})
dev.off()


#Loop through all cluster files in phi
files <- list.files( pattern="*clus*", full.names=FALSE, recursive=FALSE)
tryCatch(for(i in 1:length(files)){
  a = files[i]
  load(a) # see above, this is now called tmp
  out <-  summaryLDna(ldfile, tmp, LDmat) 
  
  out$clustersize = "NA" #make new empty column where clustersize will be written to
  
  #get the size of the cluster
  for(j in 1:length(tmp$clusters)){
  tmp2 = as.data.frame(tmp$clusters[j])
  V1 = colnames(tmp2) 
  tmp3 = tmp2 %>% tidyr::separate(V1, c("contig", "location"), "__")
  
  if(tmp3$contig[[1]]==tmp3$contig[[nrow(tmp3)]]){
    clus_size = as.numeric(tmp3$location[[nrow(tmp3)]]) -  as.numeric(tmp3$location[[1]])
  } else {
    clus_size = print("n/a")
  }
    out[j,8] =clus_size
  }
  
  write.table(out, paste0(sums,"/sum_" ,a, ".txt"))
}, warning = function(w) invokeRestart("muffleWarning"), error = function(e) {an.error.occured <<- TRUE})


##Filter the summary files 
summaries <- list.files(path=paste0(sums,"/"), pattern="sum", full.names=FALSE, recursive=FALSE)
results = list()
for( i in 1:length(summaries)){
  tst = read.table(paste0(sums,"/",summaries[i]), row.names = 1)
  tst$cluster = rownames(tst)
  tst$file = summaries[i]
  tst2 = tst %>% 
    filter(!Type=="COC") %>%
    filter(Merge.at >= Merge) %>%
    filter(nLoci >= SNPS)
  results[[i]] = tst2
  fullsumm =  do.call(rbind, results)
}
write.table(fullsumm, paste0("../sum_",tmpfile,"M",Merge,"_S",SNPS,"10.txt"), row.names = FALSE)

```
