# merge_salmon_iso.R, revised to base R on 6/14/2022 and updated 6/17/2022 
args <- commandArgs(TRUE)

read.salmon <- function(x) {
	message(sprintf("Importing %s...", x))
	y <- read.table(file = x, header = TRUE, sep = "\t", fill = TRUE)
	n <- apply(y, 1, function(x) length(which(is.na(x))))
	if (sum(n) > 0) {
		message(sprintf("** Omitting potential incomplete record(s): %s **", paste(which(n > 0) - 1, collapse = " ")))
		y <- as.data.frame(na.omit(y))
	}
	return(y)
}                

x1 <- read.csv(args[1])
x1 <- x1[,colnames(x1) %in% c("sample", args[2])]
colnames(x1)[match(args[2], colnames(x1))] <- "chain_salmon_iso"
x1$chain_salmon_iso <- trimws(x1$chain_salmon_iso)

y1.list <- lapply(x1$chain_salmon_iso, read.salmon)
y1.Name <- unique(unlist(lapply(y1.list, function(x) x$Name)))
y1 <- matrix(0, nrow = length(y1.Name), ncol = length(x1$sample))
for (i in 1:length(y1.list)) {
	y1[,i] <- y1.list[[i]][,args[3]][match(y1.Name, y1.list[[i]]$Name)]
}

x.out <- data.frame(ID = y1.Name, y1)
colnames(x.out) <- c("ID", x1$sample)
write.table(x.out, file = args[4], sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
