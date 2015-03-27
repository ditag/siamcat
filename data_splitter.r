# written by Georg Zeller EMBL Heidelberg 2012-2015
# version 0.1.0

# parameters that cannot be specified in the interface
r.seed = 2013  # TODO expose

### parse commandline arguments
suppressMessages(library('optparse'))
# define arguments
option_list = list(
  make_option(c('-s', '--srcdir'), type='character', help='Source directory of this and other utility scripts'),
  make_option('--label_in', type='character', help='Input file containing labels'),
  make_option('--train_sets', type='character', help='Output file containing training sets'),
  make_option('--test_sets', type='character', help='Output file containing test sets'),
  make_option('--num_folds', type='integer', default=10, help='Number of cross-validation folds (i.e. subsets, needs to be >= 2)'),
  make_option('--resample', type='integer', default=0, help='Resampling rounds (values <= 1 deactivate resampling)'),
  make_option('--stratify', type='logical', default=TRUE, help='Should cross-validation be stratified such that an approx. equal proportion of positive examples are contained in each subset (only for binary labels)?')
)
# parse arguments
opt = parse_args(OptionParser(option_list=option_list))
source.dir = opt$srcdir
fn.in.label = opt$label_in
fn.train.folds = opt$train_sets
fn.test.folds = opt$test_sets
num.folds = opt$num_folds
num.resample = opt$resample
stratify = opt$stratify
cat('source.dir =', source.dir, '\n')
cat('fn.in.label =', fn.in.label, '\n')
cat('fn.train.folds =', fn.train.folds, '\n')
cat('fn.test.folds =', fn.test.folds, '\n')
cat('num.folds =', num.folds, '\n')
cat('num.resample =', num.resample, '\n')
cat('stratify =', stratify, '\n')
cat('\n')
if (substr(source.dir, nchar(source.dir), nchar(source.dir)) != '/') {
  source.dir = paste(source.dir, '/', sep='')
}
source(paste(source.dir, 'utils.r', sep=''))

start.time = proc.time()[1]
set.seed(r.seed)

### read label data
# (assuming the label file has 1 column)
label = read.table(file=fn.in.label, sep='\t', header=FALSE, row.names=1, check.names=FALSE, quote='', comment.char='#')
n = rownames(label)
label = as.numeric(label[,1])
names(label) = n
con = file(fn.in.label, 'rt') 
label.header = readLines(con, 1)
close(con)
exm.ids = names(label)
# parse label description
label.info = parse.label.header(label.header)
classes = label.info$class.descr

### check arguments
if (num.resample < 1) {
  cat('\nresetting num.resample = 1 (', num.resample, ' is an invalid number of resampling rounds)\n', sep='')
  num.resample = 1
}
if (num.folds < 2) {
  cat('\nresetting num.folds = 2 (', num.folds, ' is an invalid number of folds)\n', sep='')
  num.folds = 2
}
if (num.folds >= length(label)) {
  cat('Performing leave-one-out (LOO) cross-validation\n')
  if (stratify == TRUE) {
    cat('Stratification is not possible with LOO cross-validation\n')
    stratify = FALSE
  }
  num.folds = length(label)
}

### generate files with example partitions, one line per fold
write('#Cross-validation training folds', file=fn.train.folds, append=FALSE)
write('#Cross-validation test folds', file=fn.test.folds, append=FALSE)
# in case of repeated CV, resample subsets several times
for (r in 1:num.resample) {
  ### assign data to crossvalidation folds
  foldid = rep(0, length(label))
  if (stratify) {
    for (c in classes) {
      urn = rep(1:num.folds, ceiling(sum(label==c) / num.folds))
      s = sample(urn, sum(label==c))
      foldid[label==c] = s
    }
    # try to swap some foldids to obtain equally sized folds
    pop = vector('numeric', num.folds)
    for (f in 1:num.folds) {
      pop[f] = sum(foldid==f)
    }
    while (any(pop < max(pop)-1)) {
      cat('.')
      too.big = which.max(pop)
      too.small = which.min(pop)
      idx = sample(which(foldid==too.big), 1)
      foldid[idx] = too.small
      pop[too.big] = pop[too.big] - 1
      pop[too.small] = pop[too.small] + 1
    }
  } else {
    perm = sample(1:length(label), length(label)) / length(label)
    for (f in num.folds:1) {
      foldid[perm <= f/num.folds] = f
    }
  }
  stopifnot(length(label) == length(foldid))
  stopifnot(length(unique(foldid)) == num.folds)
#  cat(foldid, '\n')

  for (f in 1:num.folds) {
    # select test examples
    test.idx = which(foldid == f)
    train.idx = which(foldid != f)
    stopifnot(length(intersect(train.idx, test.idx)) == 0)
    
    cat('Fold ', f, ' contains ', sum(foldid==f), ' examples\n', sep='')

    # append training and test examples, one line per fold
    fold.name = paste('>cv_fold', ifelse(num.resample>1, paste(f, '_rep', r, sep=''), as.character(f)), sep='')
    write.table(t(c(fold.name, exm.ids[train.idx])), file=fn.train.folds, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE, append=TRUE)
    write.table(t(c(fold.name, exm.ids[test.idx])), file=fn.test.folds, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE, append=TRUE)
  }
}
  if (stratify) {
    cat('Successfully created data split for ', num.folds, '-fold stratified cross-validation', sep='')
  } else {
    cat('Successfully created data split for ', num.folds, '-fold cross-validation', sep='')
  }
if (num.resample > 1) {
  cat(' with ', num.resample, ' times repeated resampling\n', sep='')
} else {
  cat('\n')
}
