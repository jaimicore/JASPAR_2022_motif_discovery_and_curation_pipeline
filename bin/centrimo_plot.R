args = commandArgs(trailingOnly=TRUE);

infile = args[1];
outfile = args[2];

sites = read.delim(infile, header=FALSE);
pdf(outfile);
h = hist(((sites$V8 + sites$V9) / 2) - 250, breaks=501, plot=F);
plot(h$mids, h$counts, type='l', main='', 
     xlab="Distance to peak centre", ylab="Number of motif occurrences");
null = dev.off();
message("; Exported centrimo plot PDF format")


outfile.png <- outfile
outfile.png <- gsub(outfile.png, pattern = "\\.pdf$", replacement = "")
outfile.png <- gsub(outfile.png, pattern = "\\.", replacement = "_")
outfile.png <- paste0(outfile.png, ".png")
png(outfile.png);
h = hist(((sites$V8 + sites$V9) / 2) - 250, breaks=501, plot=F);
plot(h$mids, h$counts, type='l', main='', 
     xlab="Distance to peak centre", ylab="Number of motif occurrences");
null = dev.off();
message("; Exported centrimo plot png format")