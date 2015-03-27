# written by Georg Zeller EMBL Heidelberg 2012-2015
# version 0.1.0

### default parameters
fn.train.sample = NULL
num.folds = 5
cv.stratified = TRUE
modsel.crit = 'auroc'
min.nonzero.coeff = 1

# parameters that cannot be configured externally # TODO
C.vec = 10^seq(-2,2,length=5+4+8)   # values for C (regularization strength) to be tested in model selection
r.seed = 2013                       # initiliazing the pseudo-random number generator to get reproducible results
DEBUG.CHECKS = FALSE                # performs additional checks (asserting that custom log-reg prediction code is correct)


### parse commandline arguments
suppressMessages(library('optparse'))
# define arguments
option_list = list(
  make_option(c('-s', '--srcdir'), type='character', help='Source directory of this and other utility scripts'),
  make_option('--label_in', type='character', help='Input file containing labels'),
  make_option('--feat_in', type='character', help='Input file containing features'),
  make_option('--train_sets', type='character', help='Input file specifying which examples to use for training'),
  make_option('--model', type='character', help='Output file to which the trained models will be written'),
  make_option('--num_folds', type='integer', default=10, help='Number of cross-validation folds for model selection (i.e. subsets, needs to be >= 2)'),
  make_option('--stratify', type='logical', default=TRUE, help='Should cross-validation for model selection be stratified such that an approx. equal proportion of positive examples are contained in each subset (only for binary labels)?'),
  make_option('--sel_criterion', type='character', default='auroc', help='Evaluation criterion for model selection (options: \'acc\', \'auroc\', \'auprc\', \'auroc.2\')'),
  make_option('--min_nonzero_coeff', type='integer', default=1, help='Minimum number of non-zero coefficients required for a model to be considered in model selection')
)
# parse arguments
opt = parse_args(OptionParser(option_list=option_list))
source.dir = opt$srcdir
fn.train.label = opt$label_in
fn.train.feat = opt$feat_in
fn.model = opt$model
fn.train.sample = opt$train_sets
num.folds = opt$num_folds
stratify = opt$stratify
modsel.crit = opt$sel_criterion
min.nonzero.coeff = opt$min_nonzero_coeff
cat('source.dir =', source.dir, '\n')
cat('fn.train.label =', fn.train.label, '\n')
cat('fn.train.feat =', fn.train.feat, '\n')
cat('fn.train.sample =', fn.train.sample, '\n')
cat('fn.model =', fn.model, '\n')
cat('num.folds =', num.folds, '\n')
cat('stratify =', stratify, '\n')
cat('modsel.crit =', modsel.crit, '\n')
cat('min.nonzero.coeff =', min.nonzero.coeff, '\n')
cat('\n')
if (substr(source.dir, nchar(source.dir), nchar(source.dir)) != '/') {
  source.dir = paste(source.dir, '/', sep='')
}
# optional parameters will be reset to NULL if specified as 'NULL', 'NONE' or 'UNKNOWN'
if (toupper(fn.train.sample)=='NULL' || toupper(fn.train.sample)=='NONE' || toupper(fn.train.sample)=='UNKNOWN') {
  fn.train.sample = NULL
  cat('fn.train.sample not specified: using whole data set for training\n')
}

start.time = proc.time()[1]
set.seed(r.seed)

### imports
suppressMessages(library('pROC'))
source(paste(source.dir, 'lasso.r', sep=''))
source(paste(source.dir, 'utils.r', sep=''))
source(paste(source.dir, 'eval_utils.r', sep=''))


### read training data
# features
feat = read.table(file=fn.train.feat, sep='\t', header=TRUE, row.names=1, stringsAsFactors=FALSE, check.names=FALSE, quote='', comment.char='#')
# labels (assuming the label file has 1 column)
label = read.table(file=fn.train.label, sep='\t', header=FALSE, row.names=1, check.names=FALSE, quote='', comment.char='#')
n = rownames(label)
stopifnot(n == colnames(feat))
label = as.numeric(label[,1])
names(label) = n
con = file(fn.train.label, 'rt') 
label.header = readLines(con, 1)
close(con)
# parse label description
label.info = parse.label.header(label.header)
stopifnot(label.info$type == 'BINARY')
PL = max(label.info$class.descr)
NL = min(label.info$class.descr)
# transpose feature matrix as a convenience preprocessing for LASSO
feat = t(feat)


### subselect training examples as specified in fn.train.sample (if given)
num.runs = 1
fold.name = list()
fold.exm.idx = list()
if (!is.null(fn.train.sample)) {
  num.runs = 0
  con = file(fn.train.sample, 'r') 
  input = readLines(con)
  for (i in 1:length(input)) { 
    l = input[[i]]
    if (substr(l, 1, 1) != '#') {
      num.runs = num.runs + 1
      s = unlist(strsplit(l, '\t'))
      fold.name[[num.runs]] = substr(s[1], 2, nchar(s[1]))
      fold.exm.idx[[num.runs]] = which(names(label) %in% as.vector(s[2:length(s)]))
      cat(fold.name[[num.runs]], 'contains', length(fold.exm.idx[[num.runs]]), 'training examples\n')
#      cat(fold.exm.idx[[num.runs]], '\n\n')
#    } else {
#      cat('Ignoring commented line:', l, '\n\n')
    }
  }
  close(con)
} else {
  # train on whole data set
  fold.name[[1]] = 'whole data set'
  fold.exm.idx[[1]] = names(label)
}
fold.name = unlist(fold.name)
stopifnot(length(fold.name) == num.runs)
stopifnot(length(fold.exm.idx) == num.runs)
cat('\nPreparing to train LASSO models on', num.runs, 'training set samples...\n\n')

### train one LASSO model per training sample (i.e. CV fold)
W.mat = matrix(data=NA, nrow=dim(feat)[2]+1, ncol=num.runs)
rownames(W.mat) = c(colnames(feat), 'Bias')
colnames(W.mat) = paste('M', fold.name, sep='_')
for (r in 1:num.runs) {
  cat('Training on ', fold.name[r], ' (', r, ' of ', num.runs, ')...\n', sep='')

  ### subselect examples for training
  train.label = label[fold.exm.idx[[r]]]
  train.feat = feat[fold.exm.idx[[r]],]
  stopifnot(dim(train.feat)[1] == length(train.label))
  stopifnot(all(rownames(train.feat) == names(train.label)))
  # reorder training examples so that class order is the same for all models
  exm.order = sort(train.label, index.return=TRUE)$ix
  train.label = train.label[exm.order]
  train.feat = train.feat[exm.order,]  

  ### assign training data to internal folds for model selection
  foldid = assign.fold(train.label, num.folds, cv.stratified)
  
  ### internal cross-validation for model selection
  opt.C = lasso.model.selection.cv(train.feat, train.label, foldid, C.vec, PL=PL, NL=NL,
                                   eval.crit=modsel.crit, min.nonzero=min.nonzero.coeff)
  cat('  optimal C=', opt.C, ' (', which(opt.C==C.vec), ' of ', length(C.vec), ')\n', sep='')

  ### retrain with best parameter setting
  model = lasso.train(train.feat, train.label, opt.C)
  m = model$m
  # model header string
  mh = paste('#', m$Tag, ' (', m$TypeDetail, '): [BINARY:',
    m$ClassNames[1], '=', names(label.info$class.descr)[label.info$class.descr==m$ClassNames[1]], ';',
    m$ClassNames[2], '=', names(label.info$class.descr)[label.info$class.descr==m$ClassNames[2]], ']', sep='')

  if (DEBUG.CHECKS) {
    ### some extra checks for debugging only
    suppressMessages(library('LiblineaR'))
    m.check = LiblineaR(train.feat, train.label, type=6, cost=opt.C, bias=TRUE, epsilon=0.00001)
    # the threshold for this check needs to be reasonable with respect to epsilon used for LASSO optimization
    # (3 orders of magnitudes larger seems to be alright, if in doubt recompute m.check several times and compare them)
    cat(max(abs(m$W - m.check$W)), '\n')
    stopifnot(all(abs(m$W - m.check$W) < 0.1))

    # check prediction on all examples
    m.check1 = m.check
    m.check1$W = t(m.check$W)
    p1 = lasso.predict(feat, m.check1)
    p2 = as.vector(predict(m.check, feat, proba=TRUE)$probabilities[,as.character(PL)])
    cat(max(abs(p1 - p2)), '\n')
    stopifnot(all(abs(p1 - p2) < 0.0001))
  }

  ### collect model parameters (feature weights)
  if (r==1) {
    model.header = mh
  } else {
    stopifnot(model.header == mh)
  }
  stopifnot(all(names(model$W) == rownames(W.mat)))
  W.mat[,r] = model$W

  stopifnot(!all(model$W == 0))  

  cat('\n')
}

### save models
write(model.header, file=fn.model, append=FALSE)
suppressWarnings(write.table(W.mat, file=fn.model, quote=FALSE, sep='\t', row.names=TRUE, col.names=NA, append=TRUE))
cat('Saved all trained models.\n')

cat('\nSuccessfully built ', num.runs, ' LASSO models in ', proc.time()[1] - start.time,
    ' seconds\n', sep='')
