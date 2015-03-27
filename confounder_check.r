# written by Georg Zeller EMBL Heidelberg 2012-2015
# version 0.1.0

### parse commandline arguments
suppressMessages(library('optparse'))
# define arguments
option_list = list(
  make_option(c('-s', '--srcdir'), type='character', help='Source directory of this and other utility scripts'),
  make_option('--metadata_in', type='character', help='Input file containing meta-data'),
  make_option('--label_in', type='character', help='Input file containing labels'),
  make_option('--plot', type='character', help='Output pdf file which will contain resulting plots')
)
# parse arguments
opt = parse_args(OptionParser(option_list=option_list))
source.dir = opt$srcdir
fn.in.meta = opt$metadata_in
fn.in.label = opt$label_in
fn.plot = opt$plot
cat('source.dir =', source.dir, '\n')
cat('fn.in.meta =', fn.in.meta, '\n')
cat('fn.in.label =', fn.in.label, '\n')
cat('fn.plot =', fn.plot, '\n')
cat('\n')
if (substr(source.dir, nchar(source.dir), nchar(source.dir)) != '/') {
  source.dir = paste(source.dir, '/', sep='')
}

start.time = proc.time()[1]


### imports
source(paste(source.dir, 'utils.r', sep=''))


### read label and meta- data
# labels (assuming the label file has 1 column)
label = read.table(file=fn.in.label, sep='\t', header=FALSE, row.names=1, check.names=FALSE, quote='', comment.char='#')
n = rownames(label)
label = as.numeric(label[,1])
names(label) = n
con = file(fn.in.label, 'rt') 
label.header = readLines(con, 1)
close(con)
label.info = parse.label.header(label.header)
stopifnot(label.info$type == 'BINARY')
PL = max(label.info$class.descr)
NL = min(label.info$class.descr)
# meta-data
meta.data = read.table(file=fn.in.meta, sep='\t', header=TRUE, row.names=1, check.names=FALSE, quote='', comment.char='#')
stopifnot(all(names(label) %in% rownames(meta.data)) && all(rownames(meta.data) %in% names(label)))
m = match(names(label), rownames(meta.data))
meta.data = meta.data[m,]
stopifnot(all(names(label) == rownames(meta.data)))

### generate plots with the results of confounder analysis
pdf(fn.plot)

n.idx = label==NL
n.lab = gsub('[_.-]', ' ', names(label.info$class.descr)[label.info$class.descr==NL])
p.idx = label==PL
p.lab = gsub('[_.-]', ' ', names(label.info$class.descr)[label.info$class.descr==PL])

for (m in 1:dim(meta.data)[2]) {
  mname = gsub('[_.-]', ' ', colnames(meta.data)[m])
  mname = paste(toupper(substring(mname, 1, 1)), substring(mname, 2), sep="")
  cat('checking', mname, 'as a potential confounder...\n')

  mvar = as.numeric(meta.data[,m])
  u.val = unique(mvar)
  u.val = u.val[!is.na(u.val)]
  if (length(u.val) == 1) {
    cat('  skipped because all subjects have the same value\n')
  } else if (length(u.val) <= 5) {
    cat('  using a bar plot\n')
    par(mar=c(6.1,4.1,4.1,4.1))
    ct = matrix(NA, nrow=2, ncol=length(u.val))
    c = rainbow(length(u.val), s=0.4, v=0.6, alpha=1, start=0.0)
    for (i in 1:length(u.val)) {
      ct[1,i] = sum(mvar[n.idx] == u.val[i], na.rm=TRUE)
      ct[2,i] = sum(mvar[p.idx] == u.val[i], na.rm=TRUE)
    }
    freq = t(ct)
    for (i in 1:dim(freq)[2]) {
      freq[,i] = freq[,i] / sum(freq[,i])
    }
    barplot(freq, ylim=c(0,1), main=mname, names.arg=c(n.lab, p.lab), col=c)
    p.val = fisher.test(ct)$p.value
    mtext(paste('Fisher test p-value:', format(p.val, digits=4)), side=1, line=3, at=1, adj=0)
  } else {
    cat('  using a Q-Q plot\n')
    par(mar=c(5.1,4.1,4.1,4.1))
    ax.int = c(min(mvar, na.rm=TRUE), max(mvar, na.rm=TRUE))
    qqplot(mvar[n.idx], mvar[p.idx], xlim=ax.int, ylim=ax.int, pch=16, cex=0.6,
           xlab=n.lab, ylab=p.lab, main=paste('Q-Q plot for', mname))
    abline(0, 1, lty=3)
    p.val = wilcox.test(mvar[n.idx], mvar[p.idx], exact=FALSE)$p.value
    text(ax.int[1]+0.9*(ax.int[2]-ax.int[1]), ax.int[1]+0.1*(ax.int[2]-ax.int[1]),
         paste('MWW test p-value:', format(p.val, digits=4)), pos=2)
  }
}


# close pdf device
tmp = dev.off()

cat('\nSuccessfully analyzed meta-data for potential confounding in ', proc.time()[1] - start.time, ' seconds\n', sep='')
