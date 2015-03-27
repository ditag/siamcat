# written by Georg Zeller EMBL Heidelberg 2012-2015
# version 0.1.0

### parse commandline arguments
suppressMessages(library('optparse'))
# define arguments
option_list = list(
  make_option('--feat_in', type='character', help='Input file containing features'),
  make_option('--feat_out', type='character', help='Output file to which features after selection are written'),
  make_option('--param_out', type='character', help='Output file to which normalization parameters will be written'),
  make_option('--method', type='character', default='log.std', help='Normalization method (either \"log.std\", \"rlog2.std\" or \"log.unit\")'),
  make_option('--log_n0', type='double', default=10^-8, help='Pseudocount that is added before log-transformation'),
  make_option('--sd_min_quantile', type='double', default=0.1, help='Quantile of the distribution of standard deviation of all feature that will be added to the denominator during standardization of each feature in order to avoid underestimation (only for metod==\"log.std\")'),
  make_option('--vector_norm', type='integer', default=2, help='Vector norm to use (either 1 or 2, only for method==\"log.unit\")'),
  make_option('--norm_feature', type='logical', default=FALSE, help='Normalize by feature (only for method==\"log.unit\")?'),
  make_option('--norm_sample', type='logical', default=TRUE, help='Normalize by sample (after feature normalization, only for method==\"log.unit\")?'),
  make_option('--norm_global', type='logical', default=FALSE, help='Normalize by global rescaling (only if both norm_feature and norm_sample are FALSE and only for method==\"log.unit\")?')
)
# parse arguments
opt = parse_args(OptionParser(option_list=option_list))
fn.in.feat = opt$feat_in
fn.out.feat = opt$feat_out
fn.out.param = opt$param_out
norm.method = opt$method
log.n0 = opt$log_n0
sd.min.q = opt$sd_min_quantile
n.p = opt$vector_norm
n.feature = opt$norm_feature
n.sample = opt$norm_sample
n.global = opt$norm_global
cat('fn.in.feat =', fn.in.feat, '\n')
cat('fn.out.feat =', fn.out.feat, '\n')
cat('fn.out.param =', fn.out.param, '\n')
cat('norm.method =', norm.method, '\n')
cat('log.n0 =', log.n0, '\n')
cat('sd.min.q =', sd.min.q, '\n')
cat('n.p =', n.p, '\n')
cat('n.sample =', n.sample, '\n')
cat('n.feature =', n.feature, '\n')
cat('n.global =', n.global, '\n')
cat('\n')

start.time = proc.time()[1]

### read feature data
feat = read.table(file=fn.in.feat, sep='\t', header=TRUE, row.names=1, stringsAsFactors=FALSE, check.names=FALSE, quote='', comment.char='#')

### remove features with missing values
# TODO there may be better ways of dealing with NA features
num.orig.feat = nrow(feat)
keep.idx = rowSums(is.na(feat) == 0)
if (any(!keep.idx)) {
  feat = feat[keep.idx,]
  cat('Removed ', nrow(feat)-num.orig.feat, ' features with missing values (retaining ', nrow(feat),' )\n', sep='')
}

### keep track of normalization parameters
par = list()
par$norm.method = norm.method
par$log.n0 = log.n0
par$n.p = n.p
par$n.sample = n.sample
par$n.feature = n.feature
par$n.global = n.global
par$retained.feat = rownames(feat)

### apply normalization
if (norm.method == 'rank.unit') {
  for (c in 1:ncol(feat)) {
    feat[,c] = rank(feat[,c], ties.method='average')
  }
  stopifnot(!any(is.na(feat)))
  for (c in 1:ncol(feat)) {
    feat[,c] = feat[,c] / sqrt(sum(feat[,c]^2))
  }
} else if (norm.method == 'rlog2.std') {
  suppressMessages(library('DESeq2'))
  feat = round(feat)
  cn = colnames(feat)
  rn = rownames(feat)
  d = dim(feat)
  col.data = data.frame(condition=rep(NA, ncol(feat)),
			type=rep('paired-end', ncol(feat)))  
  dds = DESeqDataSetFromMatrix(countData=feat, colData=col.data, design=~1)
  rld = rlogTransformation(dds, blind=TRUE)
  feat = assay(rld)
  stopifnot(all(d == dim(feat)))
  colnames(feat) = cn
  rownames(feat) = rn

  m = apply(feat, 1, mean)
  s = apply(feat, 1, sd)
  q = quantile(s, sd.min.q, names=FALSE)
  stopifnot(q > 0)
  for (r in 1:nrow(feat)) {
    feat[r,] = (feat[r,] - m[r]) / (s[r] + q)
  }
  par$feat.mean = m
  par$feat.adj.sd = s + q
  stopifnot(!any(is.na(feat)))
} else if (norm.method == 'log.std') {
  feat = log10(feat + log.n0)
  m = apply(feat, 1, mean)
  s = apply(feat, 1, sd)
  q = quantile(s, sd.min.q, names=FALSE)
  stopifnot(q > 0)
  # TODO needs an apply-style rewrite!
  for (r in 1:nrow(feat)) {
    feat[r,] = (feat[r,] - m[r]) / (s[r] + q)
  }
  par$feat.mean = m
  par$feat.adj.sd = s + q
  stopifnot(!any(is.na(feat)))
} else if (norm.method == 'log.unit') {
  cat('Feature sparsity before normalization: ', 100*mean(feat==0), '%\n', sep='')
  feat = log10(feat + log.n0)
  if (n.p == 1) {
    if (n.feature) {
      feat.norm.denom = vector('numeric', nrow(feat))
      # TODO needs an apply-style rewrite!
      for (r in 1:nrow(feat)) {
	feat.norm.denom[r] = sum(feat[r,])
        feat[r,] = feat[r,] / feat.norm.denom[r]
      }
      par$feat.norm.denom = feat.norm.denom
    }
    if (n.sample) {
      for (c in 1:ncol(feat)) {
        feat[,c] = feat[,c] / sum(feat[,c])
      }
    }
  } else if (n.p == 2) {
    if (n.feature) {
      feat.norm.denom = vector('numeric', nrow(feat))
      for (r in 1:nrow(feat)) {
	feat.norm.denom[r] = sqrt(sum(feat[r,]^2))
        feat[r,] = feat[r,] / feat.norm.denom[r]
      }
      par$feat.norm.denom = feat.norm.denom
    }
    if (n.sample) {
      for (c in 1:ncol(feat)) {
        feat[,c] = feat[,c] / sqrt(sum(feat[,c]^2))
      }
    }
  } else {
    stop('unknown norm!')
  }
  if (!n.feature && !n.sample && n.global) {
    global.norm.denom = max(feat)
    feat = feat / global.norm.denom
    par$global.norm.denom = global.norm.denom
  }
  cat('Feature sparsity after normalization: ', 100*mean(feat==0), '%\n', sep='')
  stopifnot(!any(is.na(feat)))
} else {
  cat('\nunrecognized norm.method, exiting!\n')
  quit(save='no', status = 1)
}

### write output
write.table(feat, file=fn.out.feat, quote=FALSE, sep='\t', row.names=TRUE, col.names=NA)

### write parameters
write('#normalization parameters', fn.out.param, append=FALSE, sep='')
for (p in 1:length(par)) {
  write(paste('###par', names(par)[p], mode(par[[p]]), sep=':'), fn.out.param, append=TRUE, sep='')
  write(par[[p]], fn.out.param, ncolumns=1, append=TRUE, sep='')
}

cat('\nSuccessfully normalized feature data in ', proc.time()[1] - start.time, ' seconds\n', sep='')
