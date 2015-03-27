# written by Georg Zeller EMBL Heidelberg 2012-2015
# version 0.1.0

### default parameters
# None

suppressMessages(library('optparse'))
# define arguments
option_list = list(
  make_option(c('-s', '--srcdir'), type='character', help='Source directory of this and other utility scripts'),
  make_option('--label', type='character', help='Input file containing labels'),
  make_option('--pred', type='character', help='Input file containing the trained classification model(s)'),
  make_option('--plot', type='character', help='Output file for plotting')
)

# parse arguments
opt = parse_args(OptionParser(option_list=option_list))
source.dir = opt$srcdir
fn.label = opt$label
fn.pred = opt$pred
fn.plot = opt$plot
cat('source.dir =', source.dir, '\n')
cat('fn.label =', fn.label, '\n')
cat('fn.pred =', fn.pred, '\n')
cat('fn.plot =', fn.plot, '\n')
cat('\n')
if (substr(source.dir, nchar(source.dir), nchar(source.dir)) != '/') {
  source.dir = paste(source.dir, '/', sep='')
}

start.time = proc.time()[1]

### imports
suppressMessages(library('pROC'))
source(paste(source.dir, 'utils.r', sep=''))
source(paste(source.dir, 'eval_utils.r', sep=''))

### read labels and predictions
label = read.table(file=fn.label, sep='\t', header=FALSE, row.names=1, check.names=FALSE, quote='', comment.char='#')
n = rownames(label)
label = as.numeric(label[,1])
names(label) = n
stopifnot(length(unique(label)) == 2)
con = file(fn.label, 'rt') 
label.header = readLines(con, 1)
close(con)
# parse label description
label.info = parse.label.header(label.header)
stopifnot(label.info$type == 'BINARY')
PL = max(label.info$class.descr)
NL = min(label.info$class.descr)
# predictions
# TODO compare header to label
pred = read.table(file=fn.pred, sep='\t', header=FALSE, row.names=1, check.names=FALSE, quote='')
if (dim(pred)[2] > 1) {
  pred = read.table(file=fn.pred, sep='\t', header=TRUE, row.names=1, check.names=FALSE, quote='')
}
pred = as.matrix(pred)
#cat(dim(pred), '\n')

### make sure that label and prediction are in the same order 
stopifnot(all(names(label) %in% rownames(pred)) && all(rownames(pred) %in% names(label)))
m = match(names(label), rownames(pred))
#cat(m, '\n')
pred = pred[m,,drop=FALSE]
stopifnot(all(names(label) == rownames(pred)))

### generate evaluation plots with ROC and precision-recall curves
pdf(fn.plot)

# ROC curve
auroc = 0
plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab='False positive rate', ylab='True positive rate', type='n')
title('ROC curve for LASSO model')
abline(a=0, b=1, lty=3)
if (dim(pred)[2] > 1) {
  aucs = vector('numeric', dim(pred)[2])
  for (c in 1:dim(pred)[2]) {
    roc.c = roc(response=label, predictor=pred[,c], ci=FALSE)
    lines(1-roc.c$specificities, roc.c$sensitivities, col=gray(runif(1,0.2,0.8)))
    aucs[c] = roc.c$auc
    cat('AU-ROC (resampled run ', c, '): ', format(aucs[c], digits=3), '\n', sep='') 
  }
  l.vec = rep(label, dim(pred)[2])
} else {
  l.vec = label
}
# average data for plotting one mean prediction curve
summ.stat = 'mean'
roc.summ = roc(response=label, predictor=apply(pred, 1, summ.stat),
                  ci=TRUE, of="se", sp=seq(0, 1, 0.05))
lines(1-roc.summ$specificities, roc.summ$sensitivities, col='black', lwd=2)
auroc = roc.summ$auc
# plot CI
x = as.numeric(rownames(roc.summ$ci))
yl = roc.summ$ci[,1]
yu = roc.summ$ci[,3]
polygon(1-c(x, rev(x)), c(yl, rev(yu)), col='#88888844', border=NA)


if (dim(pred)[2] > 1) {
  cat('Mean-pred. AU-ROC:', format(auroc, digits=3), '\n')
  cat('Averaged AU-ROC: ', format(mean(aucs), digits=3), ' (sd=', format(sd(aucs), digits=4), ')\n', sep='')
  text(0.7, 0.1, paste('Mean-prediction AUC:', format(auroc, digits=3)))
} else {
  cat('AU-ROC:', format(auroc, digits=3), '\n')
  text(0.7, 0.1, paste('AUC:', format(auroc, digits=3)))
}


# precision recall curve
plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab='Recall', ylab='Precision', type='n')
title('Precision-recall curve for LASSO model')
abline(h=mean(label==PL), lty=3)

if (dim(pred)[2] > 1) {
  aucs = vector('numeric', dim(pred)[2])
  for (c in 1:dim(pred)[2]) {
    ev = eval.classifier(pred[,c], label)
    pr = get.pr(ev)
    lines(pr$x, pr$y, col=gray(runif(1,0.2,0.8)))
    aucs[c] = calc.aupr(ev)
    cat('AU-PRC (resampled run ', c, '): ', format(aucs[c], digits=3), '\n', sep='') 
  }
  ev = eval.classifier(apply(pred, 1, summ.stat), label)
} else {
  ev = eval.classifier(pred, label)
}
pr = get.pr(ev)
lines(pr$x, pr$y, col='black', lwd=2)
aupr = calc.aupr(ev)
if (dim(pred)[2] > 1) {
  cat('Mean-pred. AU-PRC:', format(aupr, digits=3), '\n')
  cat('Averaged AU-PRC: ', format(mean(aucs), digits=3), ' (sd=', format(sd(aucs), digits=4), ')\n', sep='')
  text(0.7, 0.1, paste('Mean-prediction AUC:', format(aupr, digits=3)))
} else {
  cat('AU-PRC:', format(aupr, digits=3), '\n')
  text(0.7, 0.1, paste('AUC:', format(aupr, digits=3)))
}

# close pdf device
tmp = dev.off()

cat('\nSuccessfully evaluated predictions in ', proc.time()[1] - start.time,
    ' seconds\n', sep='')
