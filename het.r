## This is a script to filter individuals with high F coefficients > m+-3*s
## het format has the following columns:
# FID	Family ID
# IID	Within-family ID
# O(HOM)	Observed number of homozygotes
# E(HOM)	Expected number of homozygotes
# N(NM)	Number of non-missing autosomal genotypes
# F	Method-of-moments F coefficient estimate

args <- commandArgs(TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

dat <- read.table(args[1], header=F)	# Read in the EUR.het file
m <- mean(dat$V6)				# Calculate the mean  
s <- sd(dat$V6)					# Calculate the SD
valid <- subset(dat, V6 <= m+3*s & V6 >= m-3*s)	# Get any samples with F coefficient within 3 SD of the population mean
write.table(valid[,c(1,2)], args[2], quote=F, row.names=F)	# print FID and IID for valid samples
q()
