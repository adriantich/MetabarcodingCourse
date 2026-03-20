# Load library vegan
library(vegan)
library(MASS)

# Set the working directory
setwd("~/Downloads")

# Load database
db <- read.table("ULOY.final_dataset_curated.csv",head=T,sep=";",stringsAsFactors = F)
# Check names
names(db)
# Define sample columns (excluding blanks)
sample_columns <- 18:121
# Create abundance matrix
abundances <- db[,sample_columns]
# Add MOTU id as row names
rownames(abundances) <- paste(db$id,db$scientific_name)

# Check total reads
total_reads <- colSums(abundances)
total_reads[order(total_reads)]
# We will need to remove sample 07_20_S1_A_O, but before...

# Remove the initial X from sample names and generate a vector
sample_names <- gsub("X","",colnames(abundances))

# Load the metadata table
metadata <- read.table("ULOY_metadata.csv", sep=";", header=T,stringsAsFactors = F)

# Calculate metadata vectors for date and location
# This step is to ensure that the order of the dates and location is the same than in the samples
date <- NULL 
for (i in 1:length(sample_names)) date[i] <- metadata$date[metadata$original_samples==sample_names[i]]
location <- NULL  
for (i in 1:length(sample_names)) location[i] <- metadata$position[metadata$original_samples==sample_names[i]]

# Remove the failed sample and define metadata vectors as factors
# We have to remove the sample in all the vectors
abundances <- abundances[,sample_names !="07_20_S1_A_O"]
date <- as.factor(date[sample_names !="07_20_S1_A_O"])
location <- as.factor(location[sample_names !="07_20_S1_A_O"])

# Plot abundances per sample, color by date
barplot(colSums(abundances),col=date,las=2,
        ylab = "total reads",xlab="",cex.names = 0.6)

# Define function for calculating relative abundances per sample
normalize_samples <- function(db,sample_cols=1:ncol(db)){
  sumcols <- colSums(db[,sample_cols])
  for (i in sample_cols) db[,i] <- db[,i]/sumcols[i+1-sample_cols[1]]
  return(db)
}

# Calculate normalized database
abundances_relative <- normalize_samples(abundances)
# Check calculation
colSums(abundances_relative)

# Define function for drawing PCA
plot_PCA <- function(db,sample_cols=1:ncol(db),main){
  abundances <- db[,sample_cols]
  pca <-rda(t(abundances))
  pdf(paste(main,".pdf",sep=""),width=12,height=12)
  par(fig=c(0,1,0,1))
  scl=2
  scrs <- scores(pca, display = c("sites", "species"), scaling = scl)
  xlim <- with(scrs, range(species[,1], sites[,1]))*1.3
  xlim[1] <- xlim[1]-0.1
  xlim[2] <- xlim[2]+0.1
  ylim <- with(scrs, range(species[,2], sites[,2]))*1.1
  ylim[1] <- ylim[1]-0.1
  ylim[2] <- ylim[2]+0.1
  plot.window(xlim = xlim, ylim = ylim, asp = 1)
  plot(0,0, type = "n",scaling=scl, choices = F,xlim=xlim,ylim=ylim,las=1,
     xlab=paste("PC1 (",round(pca$CA$eig[1]/sum(pca$CA$eig)*100,2),"%)",sep=""),
     ylab=paste("PC2 (",round(pca$CA$eig[2]/sum(pca$CA$eig)*100,2),"%)",sep=""))
  expand.txt=1.2
  text(expand.txt*scrs$species[,1],expand.txt*scrs$species[,2],
     labels = rownames(db),col = "black", font=3,cex = .8)
  expand=.9
  arrows(0,0,expand*scrs$species[,1],expand*scrs$species[,2],length = 0.08,col="grey50")
  colors <- date
  color_names <- 1:7
  shape <- location
  shape_names <- c(21,21,21,24,25)
  for (j in 1:ncol(abundances)){
    points(scrs$sites[j,1],scrs$sites[j,2], col = "black",pch = shape_names[location[j]],
           bg = color_names[date[j]],cex=3)
    text(scrs$sites[j,1],scrs$sites[j,2], labels=colnames(abundances)[j],cex=.6)
  }
  dev.off()
}

# Plot the PCAs
plot_PCA(abundances,main="PCA ULOY absolute abundances")
plot_PCA(abundances_relative,main="PCA ULOY relative abundances")

################################################################
# Rarefaction of the abundance matrix
# We will rarefy to 225 reads (the minimum of the total reads for the remaining samples)

# This is a loop to rarefy every column in the abundance matrix to the rarefy_size
rarefy_size <- 225
abundances_rarefacted <- abundances
for (i in 1:ncol(abundances)) {
  extended_vector <- NULL
  for (motu in rownames(abundances)) extended_vector <- c(extended_vector,rep(motu,abundances[motu,i]))
  rarefacted_sample <- sample(extended_vector,size=rarefy_size,replace = T)
  for (motu in 1:nrow(abundances)) abundances_rarefacted[motu,i] <- sum(rarefacted_sample==rownames(abundances)[motu])
}
# Check the total reads for the rarefacted matrix
colSums(abundances_rarefacted)

# Plot abundances per sample, color by date
barplot(colSums(abundances_rarefacted),col=date,las=2,
        ylab = "total reads",xlab="",cex.names = 0.6)

# Plot the PCA
plot_PCA(abundances_rarefacted,main="PCA ULOY rarefacted abundances")


#######################################################

# Define function for drawing MDS

plot_MDS <- function(db,sample_cols=1:ncol(db),main,fourth_root_transform=F,method = "bray",rotate_X=F,rotate_Y=F){
  abundances <-db[,sample_cols]
  if (fourth_root_transform) distances <- vegdist(t(sqrt(sqrt(abundances))),method = method) else
  {distances <- vegdist(t(abundances),method = method)}
  mds <- isoMDS(distances,maxit = 300)
  if (rotate_X) mds$points[,1] <- -mds$points[,1]
  if (rotate_Y) mds$points[,2] <- -mds$points[,2]
  draw_mds <- function(code,codesgroup,tone,pch,cex,label){
    points(mds$points[code,], pch=pch,col="black",bg=tone,cex=cex)
    ordiellipse(mds,codes,kind="sd",conf=0.75,col=tone,show.groups=codesgroup)
    ordispider(mds,codes,label=F,col=tone,show.groups=codesgroup)
    text(mds$points[code,],labels = label,cex=.4)
  }
  codes <- as.factor(paste(date,substr(as.character(location),1,5)))
  levels(codes)
  cages_07_20 <- codes == levels(codes)[1]
  north_07_20 <- codes == levels(codes)[2]
  south_07_20 <- codes == levels(codes)[3]
  cages_09_03 <- codes == levels(codes)[4]
  north_09_03 <- codes == levels(codes)[5]
  south_09_03 <- codes == levels(codes)[6]
  cages_09_12 <- codes == levels(codes)[7]
  north_09_12 <- codes == levels(codes)[8]
  south_09_12 <- codes == levels(codes)[9]
  cages_10_01 <- codes == levels(codes)[10]
  north_10_01 <- codes == levels(codes)[11]
  south_10_01 <- codes == levels(codes)[12]
  cages_10_15 <- codes == levels(codes)[13]
  north_10_15 <- codes == levels(codes)[14]
  south_10_15 <- codes == levels(codes)[15]
  cages_10_25 <- codes == levels(codes)[16]
  north_10_25 <- codes == levels(codes)[17]
  south_10_25 <- codes == levels(codes)[18]
  cages_11_07 <- codes == levels(codes)[19]
  north_11_07 <- codes == levels(codes)[20]
  south_11_07 <- codes == levels(codes)[21]
  
  pdf(paste(main,".pdf",sep=""),width=12,height=12)
  par(fig=c(0,1,0,1))
  plot(mds$points, type="n",xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=c(min(mds$points[,1])-0.1,max(mds$points[,1])+0.1),
     ylim=c(min(mds$points[,2])-0.1,max(mds$points[,2])+0.1),
     main=c(paste(main," ",nrow(abundances)," MOTUs Metabarcoding COI",sep=""),paste("stress = ",round(mds$stress,2)," %",sep="")))
  draw_mds(cages_07_20 ,levels(codes)[1],tone="grey50",21,2,"07-20")
  draw_mds(north_07_20 ,levels(codes)[2],tone="grey50",24,2,"07-20")
  draw_mds(south_07_20 ,levels(codes)[3],tone="grey50",25,2,"07-20")
  draw_mds(cages_09_03 ,levels(codes)[4],tone=2,21,2,"09-03")
  draw_mds(north_09_03 ,levels(codes)[5],tone=2,24,2,"09-03")
  draw_mds(south_09_03 ,levels(codes)[6],tone=2,25,2,"09-03")
  draw_mds(cages_09_12 ,levels(codes)[7],tone=3,21,2,"09-12")
  draw_mds(north_09_12 ,levels(codes)[8],tone=3,24,2,"09-12")
  draw_mds(south_09_12 ,levels(codes)[9],tone=3,25,2,"09-12")
  draw_mds(cages_10_01 ,levels(codes)[10],tone=4,21,2,"10-01")
  draw_mds(north_10_01 ,levels(codes)[11],tone=4,24,2,"10-01")
  draw_mds(south_10_01 ,levels(codes)[12],tone=4,25,2,"10-01")
  draw_mds(cages_10_15 ,levels(codes)[13],tone=5,21,2,"10_15")
  draw_mds(north_10_15 ,levels(codes)[14],tone=5,24,2,"10_15")
  draw_mds(south_10_15 ,levels(codes)[15],tone=5,25,2,"10_15")
  draw_mds(cages_10_25 ,levels(codes)[16],tone=6,21,2,"10_25")
  draw_mds(north_10_25 ,levels(codes)[17],tone=6,24,2,"10_25")
  draw_mds(south_10_25 ,levels(codes)[18],tone=6,25,2,"10_25")
  draw_mds(cages_11_07 ,levels(codes)[19],tone=7,21,2,"11_07")
  draw_mds(north_11_07 ,levels(codes)[20],tone=7,24,2,"11_07")
  draw_mds(south_11_07 ,levels(codes)[21],tone=7,25,2,"11_07")
  legend("bottomright",legend=c("salmon cages","North","South"),pch=c(21,24,25))
  legend("bottomleft",legend=c("2017-07-20","2017-09-03","2017-09-12","2017-10-01",
                               "2017-10-15","2017-10-25","2017-11-07"),pch=22,pt.bg=1:7)
  # Add envfit
  julian_dates <- julian(as.POSIXlt(gsub("_","/",as.character(date))))
  env <- data.frame(as.numeric(julian_dates),location)
  en <- envfit(mds,env)
  plot(en)
  dev.off()
}

# Plot the MDS
plot_MDS(abundances,main="MDS ULOY absolute no transformation",fourth_root_transform=F)
plot_MDS(abundances,main="MDS ULOY absolute fourth root",fourth_root_transform=T)
plot_MDS(abundances_relative,main="MDS ULOY relative no transformation",fourth_root_transform=F)
plot_MDS(abundances_relative,main="MDS ULOY relative fourth root",fourth_root_transform=T)

plot_MDS(abundances_rarefacted,main="MDS ULOY rarefacted fourth root N=225 reads",fourth_root_transform=T)

# Generate a presence/absence matrix
presence_absence <- abundances_relative
presence_absence[presence_absence>0] <- 1
plot_MDS(presence_absence,main="MDS ULOY Presence-Absence",fourth_root_transform=F,method="jaccard")

