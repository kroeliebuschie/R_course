#!/usr/bin/Rscript
# student name: Mehdi Nellen
# student nr: 3701263

# style guide used: https://google-styleguide.googlecode.com/svn/trunk/Rguide.xml

options(stringsAsFactors = FALSE)

# Packages for plotting
require(ggplot2)
require(grid)
require(gridBase)

# URLs with online data
spt7.url <- "http://bioinformatics.holstegelab.nl/courses/intro_R_2015_jan/spt7.txt"
SGD.url  <- "http://bioinformatics.holstegelab.nl/courses/intro_R_2015_jan/SGD.txt"



############################
######    PLOTTING    ######
############################

# Read data
spt7.table <- read.delim(spt7.url, row.names = 1)

# calculate expression differences and make logical vectors
# which can be used to make subsets
trans.diff <- spt7.table$spt7.del - spt7.table$wt
ind.up     <- trans.diff >=  1
ind.down   <- trans.diff <= -1

#### PLOT
# Most of the code here is for combining the ggplot
# with the base plot.

# create pdf
pdf("spt7_scatterd_plot.pdf", width = 20, height = 10)
# settings to combine Base and grid plots
plot.new() 
gl <- grid.layout(1,2)
vp.1 <- viewport(layout.pos.col = 1, layout.pos.row = 1) 
vp.2 <- viewport(layout.pos.col = 2, layout.pos.row = 1) 
pushViewport(viewport(layout=gl))
pushViewport(vp.1)

### Plotting with R's graphics package
par(new = TRUE, fig = gridFIG(), pch = 16)
plot(wt~spt7.del, data = spt7.table,
     xlab = "Transcriptional profile of SPT7 deletion", 
     ylab = "Transcriptional profile of wild-type",
     main = "Using R's standard plotting package")
points(spt7.table[ ind.up  , ], col = "red")
points(spt7.table[ ind.down, ], col = "green3")
points(spt7.table[ "SPT7"  , ], col = "blue")
with(spt7.table[ "SPT7"  , ], text(wt~spt7.del, labels = "SPT7", pos = 2))
legend("bottomright", pch = 16, 
       c("upregulated", "unchanged", "downregulated", "SPT7"), 
       col = c("red", "black", "green3", "blue"))
popViewport() #finished first viewport

### Plotting with ggplot
# Make an extra column with conditions for the plotting
spt7.table$transcription <- "unchanged"
spt7.table$lab           <- ""
spt7.table[ ind.up  , "transcription" ] <- "upregulated"
spt7.table[ ind.down, "transcription" ] <- "downregulated"
spt7.table[ "SPT7"  , "transcription" ] <- "SPT7"
spt7.table[ "SPT7"  , "lab" ]           <- "SPT7"


# Make the plot
pushViewport(vp.2) #start with a new view port
p <- ggplot(spt7.table, aes(x = spt7.del, y = wt, color=transcription)) + 
  geom_point(shape = 16, size = 3) +   
  scale_color_manual(values = c("green3", "blue", "black", "red")) +
  geom_text(aes(label = lab),hjust = 1, vjust = 1, show_guide = FALSE) +
  ggtitle("Using the ggplot2 package") +
  theme(plot.title   = element_text(lineheight = .8, face="bold", size = rel(2)),
        axis.title.x = element_text(lineheight = .8, size = rel(2)),
        axis.title.y = element_text(lineheight = .8, size = rel(2)))

print(p, newpage = FALSE) # This will add it to the existing plot
popViewport(1)
dev.off()


############################
####   COUNTING GENES   ####
############################

## This assignment has been done in two ways. The first method is 
## short and elegant and applies to this data, the second method 
## is applicable to other data sets because of the use of functions.
## The second method also deals with the genes who don't have a 
## chromosome.

# Read data
SGD.table <- read.delim(SGD.url, row.names = 1)

### METHOD 1 ###

# Count genes for every chromosome
n.chr <- tapply(rownames(SGD.table), SGD.table$chr, length)

# In case the chromosomes are factors (which you might get when 
# importing them this way) the previous result can also be 
# achieved like this:
n.chr <- table(factor(SGD.table$chr))

# Concatenate the result
cat(paste("Chromosome", rownames(n.chr), "contains", 
          n.chr, "gene(s)", collapse = "\n"))

## This will also print the amount of genes with no chromosome

### METHOD 2 ###
# In this method for loops are used to show they also work.
# Duplicate gene names are removed (they are not present 
# in the data). It is also possible to remove NA's and
# empty names.

# Function to count genes in chromosomes
nGenes <- function(chr, genes, na.rm = FALSE, empty.rm = FALSE) {
  # Count the amount of non duplicated genes on each chromosome
  #
  # Args: 
  #   chr     : character vector specifying chromosome names. 
  #   genes   : character vector with length and order of 
  #             chr specifying gene names.
  #   na.rm   : a logical value indicating whether NA values .
  #             should be stripped before the computation proceeds.
  #   empty.rm: a logical value indicating whether empty  values 
  #             should be stripped before the computation proceeds.
  #
  # Returns:
  #   A list with number of genes as value and chromosome name as 
  #   name. The list has class nGenesChrom.
  if(length(chr) != length(genes)) {
    stop("chr and genes differ in length")
  }
  
  # Remove NA's if desired
  if(na.rm) {
    to.rm <- !(is.na(chr) | is.na(genes))
    chr   <- chr[to.rm]
    genes <- genes[to.rm]
  }
  
  # remove empty chromosome names
  if(empty.rm) {
    to.rm <- !(chr == "" | genes == "")
    chr   <- chr[to.rm]
    genes <- genes[to.rm]
  } else if(any(chr == "", genes == "")) {
    warning("There are empty strings in your data.", call. = FALSE)
  }
  
  # Make an empty list
  chr.list <- list()
  
  #fill the empty list in a for loop
  for(chrom in unique(chr)) {
    n.gen <- sum(!duplicated(genes[chr == chrom]))
    chr.list[[chrom]] <- n.gen
  }
  chr.list <- structure(chr.list, class = "nGenesChrom")
  return(chr.list)
}

print.nGenesChrom <- function(list) {
  # Outputs the list, concatenating the representations. 
  #
  # Args: 
  #   list: a list created by nGenes()
  #
  # Returns:
  #   Prints results of nGenes.
  
  # prevent errors with empty names
  names(list)[names(list) == ""] <- "NoChrName"
  
  cat("Results:\n")
  for (name in names(list)) {
    # The if statement discriminates 1 so it does not get the plural word of "gene"
    if(list[[name]] == 1) {
      cat("\tChromosome", name, "conatins", list[[name]], "gene", "\n")
    } else {
      cat("\tChromosome", name, "conatins", list[[name]], "genes", "\n")
    }
  }
}

# run the functions on the data
result.list <- nGenes(SGD.table$chromosome, rownames(SGD.table), 
                   na.rm = TRUE, empty.rm = TRUE)

# result.list is now a list of class "nGenesChrom" 
# output can be printed using print or nothing as bellow
# but under the hood it uses cat to return some output
print(result.list)
result.list

# The unspecified chromosome got removed (you can avoid this by setting empty.rm = False)
# this will give a warning



############################
####   REPLACING NA's   ####
############################

## The following function can be made shorter (and
## faster) by removing the if/else statement & error
## handling.

NAreplacR <- function(x) {
  # Replaces NA's in a vector with the mean
  # of all real values in that vector.
  #
  # Args: 
  #   x: a numeric vector
  #
  # Returns:
  #   A vector with all NA's replaced with the mean
  if(!is.numeric(x) || !is.vector(x) ){
    # Error handling
    stop("Argument should be a numeric vector")
  } 
  x[is.na(x)] <- mean(x, na.rm = TRUE)           # replace NA's with the mean
  return(x)
}
