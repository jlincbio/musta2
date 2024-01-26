args <- commandArgs(trailing = TRUE)
listFL <- unlist(strsplit(args[2], split = ","))
#listFL <- sort(list.files(path = getwd(), full.names = TRUE, recursive = TRUE, pattern = "fl_count.txt"))
fl <- vector(mode = "list", length = length(listFL))
names(fl) <- gsub("\\/fl_count.*$", "", gsub("^\\/", "", gsub(getwd(), "", listFL)))
for (i in 1:length(fl)) {
	fl[[i]] <- read.table(file = listFL[[i]], sep = "\t", header = TRUE)
}
pbid <- gtools::mixedsort(unique(unlist(sapply(fl, function(x) x$pbid))))
flm <- matrix(NA, nrow = length(pbid), ncol = length(fl))
for (i in 1:length(fl)) {
	flm[,i] <- fl[[i]]$count_fl[match(pbid, fl[[i]]$pbid)]
}
flm <- apply(flm, 2, function(x) ifelse(is.na(x), yes = 0, no = x))
flx <- data.frame(pbid = pbid, flm)
colnames(flx) <- c("id", names(fl))
write.csv(flx, file = args[1], row.names = FALSE, quote = FALSE)
message(sprintf("# Merged FLcounts file written to %s.\n", args[1]))
