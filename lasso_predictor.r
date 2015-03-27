# written by Georg Zeller EMBL Heidelberg 2012-2015
# version 0.1.0

### default parameters
fn.test.sample = NULL
fn.test.label = NULL

### parse commandline arguments
suppressMessages(library('optparse'))
# define arguments
option_list = list(
  make_option(c('-s', '--srcdir'), type='character', help='Source directory of this and other utility scripts'),
  make_option('--label_in', type='character', help='Input file containing labels'),
  make_option('--feat_in', type='character', help='Input file containing features'),
  make_option('--test_sets', type='character', help='Input file specifying which examples to use for testing'),
  make_option('--model', type='character', help='Input file containing the trained models'),
  make_option('--pred', type='character', help='Output file to which predictions will be written')
)
# parse arguments
opt = parse_args(OptionParser(option_list=option_list))
source.dir = opt$srcdir
fn.test.label = opt$label_in
fn.test.feat = opt$feat_in
fn.model = opt$model
fn.pred = opt$pred
fn.test.sample = opt$test_sets
cat('source.dir =', source.dir, '\n')
cat('fn.test.label =', fn.test.label, '\n')
cat('fn.test.feat =', fn.test.feat, '\n')
cat('fn.test.sample =', fn.test.sample, '\n')
cat('fn.model =', fn.model, '\n')
cat('fn.pred =', fn.pred, '\n')
cat('\n')
if (substr(source.dir, nchar(source.dir), nchar(source.dir)) != '/') {
  source.dir = paste(source.dir, '/', sep='')
}
# optional parameters will be reset to NULL if specified as 'NULL', 'NONE' or 'UNKNOWN'
if (is.null(fn.test.sample) || toupper(fn.test.sample)=='NULL' || toupper(fn.test.sample)=='NONE' || toupper(fn.test.sample)=='UNKNOWN') {
  fn.test.sample = NULL
  cat('fn.test.sample not specified: applying LASSO model(s) on whole data set\n')
}
if (is.null(fn.test.label) || toupper(fn.test.label)=='NULL' || toupper(fn.test.label)=='NONE' || toupper(fn.test.label)=='UNKNOWN') {
  fn.test.label = NULL
  cat('fn.test.label not specified: skipping evaluation\n')
}

start.time = proc.time()[1]

### imports
source(paste(source.dir, 'utils.r', sep=''))
source(paste(source.dir, 'eval_utils.r', sep=''))
source(paste(source.dir, 'lasso.r', sep=''))

### read test data and the trained LASSO model(s)
# features
feat = read.table(file=fn.test.feat, sep='\t', header=TRUE, row.names=1, stringsAsFactors=FALSE, check.names=FALSE, quote='', comment.char='#')
# transpose feature matrix as a convenience preprocessing for LASSO
feat = t(feat)

# LASSO model (coefficients)
model = NULL
model$W = read.table(file=fn.model, sep='\t', header=TRUE, row.names=1, stringsAsFactors=FALSE, check.names=FALSE, quote='')
num.runs = dim(model$W)[2]
stopifnot(dim(model$W)[1]-1 == dim(feat)[2])
# parse model header
con = file(fn.model, 'r')
model.header = readLines(con, 1)
close(con)
model.header = parse.model.header(model.header)
label.header = model.header$label.header
# parse label description
label.info = parse.label.header(label.header)
stopifnot(label.info$type == 'BINARY')
PL = max(label.info$class.descr)
NL = min(label.info$class.descr)

### subselect test examples as specified in fn.test.sample (if given)
fold.name = list()
fold.exm.idx = list()
if (!is.null(fn.test.sample)) {
  con = file(fn.test.sample, 'r') 
  input = readLines(con)
  m.idx = 0
  for (i in 1:length(input)) { 
    l = input[[i]]
    if (substr(l, 1, 1) != '#') {
      m.idx = m.idx + 1
      s = unlist(strsplit(l, '\t'))
      fold.name[[m.idx]] = substr(s[1], 2, nchar(s[1]))
      fold.exm.idx[[m.idx]] = which(rownames(feat) %in% as.vector(s[2:length(s)]))
#      cat(fold.name[[m.idx]], 'contains', length(fold.exm.idx[[m.idx]]), 'test examples\n')
#      cat(fold.exm.idx[[m.idx]], '\n\n')
    }
  }
  close(con)
  stopifnot(length(fold.name) == num.runs)
  stopifnot(length(fold.exm.idx) == num.runs)
  stopifnot(all(paste('M', unlist(fold.name), sep='_') == colnames(model$W)))
} else {
  # apply each LASSO model on whole data set when only single test set is given
  for (r in 1:num.runs) {
    fold.name[[r]] = paste('whole data set predicted by model', r)
    fold.exm.idx[[r]] = rownames(feat)
  }
}
fold.name = unlist(fold.name)
cat('\nPreparing to make predictions with', num.runs, 'LASSO model(s)...\n')

### apply one LASSO model per test sample (i.e. CV fold)
# predictions are made on a concatenation of test examples from all test samples
pred = NULL
fold.pred.idx = list()
for (r in 1:num.runs) {
  cat('Applying ', colnames(model$W)[r], ' on ', fold.name[r], ' (', r, ' of ', num.runs, ')...\n', sep='')
  # subselect appropriate model
  m = model
  m$W = m$W[,r]

  # subselect test examples
  test.feat = feat[fold.exm.idx[[r]],,drop=FALSE]

  # apply LASSO model
  p = lasso.predict(test.feat, m)
  pred = c(pred, p)
  fold.pred.idx[[r]] = (length(pred)-length(p)+1):length(pred)
}
cat('\nTotal number of predictions made:', length(pred), '\n')

if (!is.null(fn.test.label)) {
  ### if test labels are given do some evaluation as well
  label = read.table(file=fn.test.label, sep='\t', header=FALSE, row.names=1, check.names=FALSE, quote='', comment.char='#')
  n = rownames(label)
  stopifnot(n == rownames(feat))
  label = as.numeric(label[,1])
  names(label) = n
  con = file(fn.test.label, 'rt') 
  label.header = readLines(con, 1)
  close(con)
  # parse label description
  label.info = parse.label.header(label.header)
  stopifnot(label.info$type == 'BINARY')
  PL = max(label.info$class.descr)
  NL = min(label.info$class.descr)
  # get the appropriate labels for all test sets
  test.label = NULL
  aucs = vector('numeric', num.runs)
  for (r in 1:num.runs) {
    l = label[fold.exm.idx[[r]]]
    test.label = c(test.label, l)
    p.idx = (length(test.label)-length(l)+1):length(test.label)
    # accuracy of individual test sets
    if (length(unique(test.label[p.idx])) == 2) {
      ev = eval.classifier(pred[p.idx], as.vector(test.label[p.idx]))
      aucs[r] = calc.auroc(ev)
    }
  }
  stopifnot(length(test.label) == length(pred))
  stopifnot(names(test.label) == names(pred))

  # in case of cross-validation there should be exactly one prediction per labeled example,
  # so we reorder them according to the order of label
  if (length(label) == length(pred) && all(names(label) %in% names(pred)) && all(names(pred) %in% names(label))) {
    m = match(names(label), names(pred))
    pred = pred[m]
    test.label = test.label[m]
    stopifnot(all(names(label) == names(pred)))
  }
  
  # test accuracy of combined test set
  c.auc = NA
  if (length(unique(test.label)) == 2) {
    ev = eval.classifier(pred, as.vector(test.label))
    c.auc = calc.auroc(ev)
  }
  cat('Combined test AUC = ', format(c.auc, digits=3),
      ' (m=', format(mean(aucs, na.rm=TRUE), digits=3),
      ', s.d.=', format(sd(aucs, na.rm=TRUE), digits=3), ')\n', sep='')
}

### reformat predictions in case models were trained in repeated cross-validation
if (length(unique(names(pred))) < length(pred)) {
  ref.names = NULL
  if (any(substr(fold.name,1,14) == 'whole data set')) {
    r.idx = as.numeric(sapply(strsplit(fold.name, 'predicted by model '), '[[', 2))
    runs = sort(unique(r.idx))
    stopifnot(all(runs == 1:length(runs)))
    if (!is.null(fn.test.label)) {
      ref.names = names(label)
    } else {
      ref.names = unique(names(pred))
    }
  } else {
    r.idx = as.numeric(sapply(strsplit(fold.name, 'rep'), '[[', 2))
    runs = sort(unique(r.idx))
    stopifnot(all(runs == 1:length(runs)))
    if (!is.null(fn.test.label)) {
      ref.names = names(label)
    } else {
      ref.names = names(pred)[unlist(fold.pred.idx[r.idx==1])]
    }
  }
#  cat(ref.names, '\n\n')
#  cat(names(label), '\n\n')
#  cat(names(pred), '\n\n')

  pred.mat = matrix(data=NA, nrow=length(ref.names), ncol=length(runs))
  rownames(pred.mat) = ref.names
  if (any(substr(fold.name,1,14) == 'whole data set')) {
    colnames(pred.mat) = paste('Model', runs, sep='') 
  } else {
    colnames(pred.mat) = paste('CV_rep', runs, sep='') 
  }
  
  for (r in runs) {
    idx = which(r.idx == r)
    p = unlist(fold.pred.idx[idx])
    m = match(names(pred)[p], ref.names)
#    cat(sort(m), '\n\n')
#    cat(length(m), '\n\n')
#    cat(length(label), '\n\n')

    stopifnot(all(sort(m) == 1:length(label)))
    pred.mat[m,r] = pred[p]
    stopifnot(all(names(pred)[p] == rownames(pred.mat)[m]))
  }
  c = cor(pred.mat, method='spearman')
  cat('\nCorrelation between predictions from repeated CV:\n')
  cat('Min: ', min(c), ', Median: ', median(c), ', Mean: ', mean(c), '\n', sep='')
}

### save prediction
pred.header = paste('#Predictions for ', PL, ':', names(label.info$class.descr)[label.info$class.descr==PL],
  ' [', label.header, ']', sep='')
write(pred.header, file=fn.pred, append=FALSE)
if (length(unique(names(pred))) < length(pred)) {
  suppressWarnings(write.table(pred.mat, file=fn.pred, quote=FALSE, sep='\t', row.names=TRUE, col.names=NA, append=TRUE))
} else {
  write.table(pred, file=fn.pred, quote=FALSE, sep='\t', row.names=TRUE, col.names=FALSE, append=TRUE)
}
cat('\nSaved all predictions\n')

cat('\nSuccessfully applied LASSO model in ', proc.time()[1] - start.time,
    ' seconds\n', sep='')
