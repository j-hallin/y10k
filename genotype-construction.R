##################################################################
## This script reads in the genotypes of the 172 segregants and ##
## creates all possible diploid hybrids                         ##
##################################################################

######
## Function for creating the crosses
####
crossCreate <- function(data.list, strain) {
    names <- gsub('X', '', names(data))
    cross.names <- paste(names[87:172], names[strain], sep='-')
    
    crossName <- list()
    for (i in 1:(length(data)/2)) {
        crossName[i] <- paste(data.list[[strain]], data.list[[86+i]], collapse=' ', sep='')
    }
    
    ##split the strings so that each marker is separated
    ##have to do it like this in order for the list to work correctly
    final <- list()
    for (j in 1:length(crossName)) {
        final[j] <- strsplit(crossName[[j]], ' ')
    }
    #create the data.frame with cross names
    
    final.2 <- data.frame(matrix(unlist(final), nrow=length(data[,1]), byrow=F))
    #final <- data.frame(final)
    names(final.2) <- cross.names
    return(final.2)
}

## Read in the data
data <- read.csv('~/Dropbox/myResearch/strains/2way_4way_segregants/originalGenotypes/2way/F12segregants/genotypes_mpileupSegregantsWithParents_jun7-2way-SegSites_andersbe.numbered.txt', sep='\t', header=T)
nam.pos <- data[,1:3]
data <- data[,-c(1:3)]


## Making a list of the columns
## Each individual is its own vector in the list
data.list <- list()
for (i in 1:length(data[1,])) {
  data.list[i] <- list(data[,i])
}

## Creating the crosses
crosses <- lapply(as.list(1:86), crossCreate, data.list = data.list)

hybrid.genotypes <- do.call(cbind, crosses)

##Since nr 45 was contaminated I removed it from the experiments,
##I will also have to remove it here
rm45 <- grep('45$', names(hybrid.genotypes))
hybrid.genotypes <- hybrid.genotypes[-rm45]

##And lastly I put back the marker names and everything
hybrid.genotypes <- hybrid.genotypes[,with(hybrid.genotypes, order(names(hybrid.genotypes)))]
hybrid.genotypes <- cbind(nam.pos, hybrid.genotypes)
hybrid.genotypes$X <- gsub(pattern = 'chr', replacement = '', hybrid.genotypes$X)


