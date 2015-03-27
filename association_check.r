# written by Georg Zeller EMBL Heidelberg 2012-2015
# version 0.1.0

### parse commandline arguments
suppressMessages(library('optparse'))
# define arguments
option_list = list(
  make_option(c('-s', '--srcdir'), type='character', help='Source directory of this and other utility scripts'),
  make_option('--label_in', type='character', help='Input file containing labels'),
  make_option('--feat_in', type='character', help='Input file containing features'),
  make_option('--plot', type='character', help='Output pdf file which will contain resulting plots'),
  make_option('--col_scheme', type='character', default='RdYlBu', help='Color scheme'),
  make_option('--alpha', type='double', default=0.05, help='Significance level: only features with p-values < alpha will be reported'),
  make_option('--min_fc', type='double', default=0, help='Fold-change cutoff: only features with absolute log-10 fold change > min_fc will be reported'),
  make_option('--mult_test', type='character', default='fdr', help='Method to correct for multiple testing (one of \"fdr\", \"holm\", \"bonferroni\", \"BHY\", or \"none\")'),
  make_option('--detect_limit', type='double', default=10^-8, help='Lower detection limit for feature values (for log-transform and plots)'),
  make_option('--max_show', type='integer', default=50, help='Maximum number of significant features to be shown in result plots')
)
# parse arguments
opt = parse_args(OptionParser(option_list=option_list))
source.dir = opt$srcdir
fn.in.label = opt$label_in
fn.in.feat = opt$feat_in
fn.plot = opt$plot
color.scheme = opt$col_scheme
alpha = opt$alpha
min.fc = opt$min_fc
mult.corr = opt$mult_test
detect.lim = opt$detect_limit
max.show = opt$max_show
cat('source.dir =', source.dir, '\n')
cat('fn.in.label =', fn.in.label, '\n')
cat('fn.in.feat =', fn.in.feat, '\n')
cat('fn.plot =', fn.plot, '\n')
cat('color.scheme =', color.scheme, '\n')
cat('alpha =', alpha, '\n')
cat('min.fc =', min.fc, '\n')
cat('mult.corr =', mult.corr, '\n')
cat('detect.lim =', detect.lim, '\n')
cat('max.show =', max.show, '\n')
cat('\n')
if (substr(source.dir, nchar(source.dir), nchar(source.dir)) != '/') {
  source.dir = paste(source.dir, '/', sep='')
}

start.time = proc.time()[1]

### further parameters (not exposed)
sort.by = 'pv' # 'fc' or 'pv' (fold change or p-value respectively)

### imports
suppressMessages(library('colorRamps'))
suppressMessages(library('RColorBrewer'))
source(paste(source.dir, 'utils.r', sep=''))
source(paste(source.dir, 'eval_utils.r', sep=''))


### read label, feature and meta-data
# features
feat = as.matrix(read.table(file=fn.in.feat, sep='\t', header=TRUE, row.names=1, stringsAsFactors=FALSE, check.names=FALSE, quote='', comment.char='#'))
# labels (assuming the label file has 1 column)
label = read.table(file=fn.in.label, sep='\t', header=FALSE, row.names=1, check.names=FALSE, quote='', comment.char='#')
n = rownames(label)
stopifnot(n == colnames(feat))
label = as.numeric(label[,1])
names(label) = n
con = file(fn.in.label, 'rt') 
label.header = readLines(con, 1)
close(con)
# parse label description
label.info = parse.label.header(label.header)
stopifnot(label.info$type == 'BINARY')
PL = max(label.info$class.descr)
NL = min(label.info$class.descr)

### some color pre-processing
if (color.scheme == 'matlab') {
  color.scheme = matlab.like(100)
} else {
  # TODO check for valid param!
  color.scheme = rev(colorRampPalette(brewer.pal(11, color.scheme))(100))
}
col.p = color.scheme[length(color.scheme)-4]
col.n = color.scheme[1+4]


n.idx = label==NL
n.lab = gsub('[_.-]', ' ', names(label.info$class.descr)[label.info$class.descr==NL])
p.idx = label==PL
p.lab = gsub('[_.-]', ' ', names(label.info$class.descr)[label.info$class.descr==PL])

# calculate abundance fold change and test for significant associations
# between features and labes using Wilcoxon test
p.val = vector('numeric', nrow(feat))
fc = vector('numeric', nrow(feat))
for (i in 1:nrow(feat)) {
  fc[i] = median(log10(feat[i,p.idx] + detect.lim)) - median(log10(feat[i,n.idx] + detect.lim))
  p.val[i] = wilcox.test(feat[i,n.idx], feat[i,p.idx], exact=FALSE)$p.value
}
if (mult.corr == 'none') {
  p.adj = p.val
} else if (tolower(mult.corr) == 'bonferroni') {
  p.adj = p.adjust(p.val, method='bonferroni')
} else if (tolower(mult.corr) == 'holm') {
  p.adj = p.adjust(p.val, method='holm')
} else if (tolower(mult.corr) == 'fdr' ) {
  p.adj = p.adjust(p.val, method='fdr')
} else if (tolower(mult.corr) == 'bhy') {
  p.adj = p.adjust(p.val, method='BY')
} else {
  stop('Unknown multiple testing correction method:', mult.corr)
}
cat('Found', sum(p.adj < alpha, na.rm=TRUE), 'significant associations at a significance level <', alpha, '\n')
idx = which(p.adj < alpha)
if (min.fc > 0) {
  idx = which(p.adj < alpha & abs(fc) > min.fc)
  cat('Found', length(idx), 'significant associations with absolute log10 fold change >', min.fc, '\n')
}

if (length(idx) > 0) {
  if (sort.by == 'fc') {
    idx = idx[order(fc[idx], decreasing=FALSE)]
  } else if (sort.by == 'pv') {
    idx = idx[order(p.adj[idx], decreasing=TRUE)]
  } else {
    cat('Unknown sorting option:', sort.by, 'order by p-value...\n')
    idx = idx[order(p.adj[idx], decreasing=TRUE)]
  }
  for (i in idx) {
    cat(sprintf('%-40s', rownames(feat)[i]), 'p-value:', format(p.adj[i], digits=4), '\n')
  }
  # truncated the list for the following plots
  if (length(idx) > max.show) {
    idx = idx[(length(idx)-max.show+1):length(idx)]
    cat('Truncating the list of significant associations to the top', max.show, '\n')  
  }

  # compute single-feature AUCs
  cat('\nCalculating the area under the ROC curve for each significantly associated feature\n')
  aucs = vector('numeric', nrow(feat))
  for (i in idx) {
    f = feat[i,]
    ev = eval.classifier(f, label)
    aucs[i] = calc.auroc(ev)
    if (aucs[i] < 0.5) {
      # convert for negative associations
      aucs[i] = 1-aucs[i]
    }
  }
  for (i in idx) {
    cat(sprintf('%-40s', rownames(feat)[i]), aucs[i], '\n')
  }


  ### generate plots with significant associations between features and labels
  pdf(fn.plot, paper='special', height=8.27, width=11.69) # format: A4 landscape
  lmat = cbind(1,2,3)
  layout(lmat, widths=c(0.6,0.2,0.2))

  # range plot of data with p-values
  x = log10(as.matrix(feat[idx, p.idx, drop=FALSE]) + detect.lim)
  y = log10(as.matrix(feat[idx, n.idx, drop=FALSE]) + detect.lim)

  par(mar=c(5.1, 25.1, 4.1, 5.1))
  col = c(paste(col.n, '77', sep=''), paste(col.p, '77', sep=''), 'gray')

  plot.data.range(x, y, rownames(feat)[idx], x.col=col[2], y.col=col[1],
                  x.suff=paste(' (', p.lab, ')', sep=''), y.suff=paste(' (', n.lab, ')', sep=''))  

  p.val.annot = formatC(p.adj[idx], format='E', digits=2)
  for (i in 1:length(p.val.annot)) {
    mtext(p.val.annot[i], 4, line=0.5, at=i, las=1, cex=min(0.7, 1-(length(idx)/100)))
  }
  mtext('Adj. p-value', 4, line=0.5, at=length(idx)+1, cex=0.7, font=2, las=1)
  if (sum(p.adj < alpha, na.rm=TRUE) <= max.show) {
    title(main='Differentially abundant features', xlab='Abundance (log10-scale)')
  } else {
    title(main=paste('Differentially abundant features\ntruncated to the top', max.show),
          xlab='Abundance (log10-scale)')
  }

  # plot fold changes
  par(mar=c(5.1, 2.1, 4.1, 2.1))
  bcol = ifelse(fc[idx] > 0, col[2], col[1])
  mn = floor(min(fc[idx]))
  mx = ceiling(max(fc[idx]))
  mx = max(abs(mn), abs(mx))
  if (!is.finite(mx)) {
    mx = 10
  }
  mn = -mx
  plot(NULL, xlab='', ylab='', , xaxs='i', yaxs='i', axes=FALSE,
       xlim=c(mn, mx), ylim=c(0.2, length(idx)+0.2), type='n')
  barplot(fc[idx], horiz=TRUE, width=0.6, space=2/3, col=bcol, axes=FALSE, add=TRUE)
  ticks = mn:mx
  for (v in ticks) {
    abline(v=v, lty=3, col='lightgrey')
  }
  tick.labels = formatC(10^ticks, format='E', digits=0)
  axis(side=1, at=ticks, labels=tick.labels, cex.axis=0.7)
  title(main='Fold change', xlab='FC (log10-scale)')

  # plot single-feature AUCs
  par(mar=c(5.1, 1.1, 4.1, 3.1))
  plot(NULL, xlab='', ylab='', , xaxs='i', yaxs='i', axes=FALSE,
       xlim=c(0.5,1), ylim=c(0.5, length(idx)+0.5), type='n')
  ticks = seq(0.5, 1.0, length.out=6)
  for (v in ticks) {
    abline(v=v, lty=3, col='lightgrey')
  }
  for (y in 1:length(idx)) {
    i = idx[y]
    points(aucs[i], y, pch=18, col=bcol[y])
    points(aucs[i], y, pch=5, col='black', cex=0.9)
  }
  axis(side=1, at=ticks, cex.axis=0.7)
  title(main='Feature AUCs', xlab='AU-ROC')

  # close pdf device
  tmp = dev.off()
}

cat('\nSuccessfully analyzed statistically significant associations between individual features and labels in ', proc.time()[1] - start.time, ' seconds\n', sep='')