# written by Georg Zeller EMBL Heidelberg 2012-2015
# version 0.1.0

### parse commandline arguments
suppressMessages(library('optparse'))
# define arguments
option_list = list(
  make_option('--feat_in', type='character', help='Input file containing features'),
  make_option('--feat_out', type='character', help='Output file to which features after selection are written'),
  make_option('--method', type='character', default='abundance', help='Filtering method (one of \"abundance\", \"cum.abundance\", or \"prevalence\")'),
  make_option('--cutoff', type='double', default=0.001, help='abundance / prevalence cutoff applied for filtering'),
  make_option('--recomp_prop', type='logical', default=FALSE, help='Should relative abundances be be recomputed?'),
  make_option('--rm_unmapped', type='logical', default=TRUE, help='Should the abundance of unmapped reads be removed?')
)
# parse arguments
opt = parse_args(OptionParser(option_list=option_list))
source.dir = opt$srcdir
fn.in.feat = opt$feat_in
fn.out.feat = opt$feat_out
filter.method = opt$method
abundance.cutoff = opt$cutoff
recomp.prop = opt$recomp_prop
rm.unmapped = opt$rm_unmapped
cat('fn.in.feat =', fn.in.feat, '\n')
cat('fn.out.feat =', fn.out.feat, '\n')
cat('filter.method =', filter.method, '\n')
cat('abundance.cutoff =', abundance.cutoff, '\n')
cat('recomp.prop =', recomp.prop, '\n')
cat('rm.unmapped =', rm.unmapped, '\n')
cat('\n')

start.time = proc.time()[1]


#### read feature data and preprocess
feat = as.matrix(read.table(file=fn.in.feat, header=TRUE, row.names=1, 
                            sep='\t', stringsAsFactors=FALSE, check.names=FALSE, quote='', comment.char='#'))
  
if (recomp.prop) {
  # recompute relative abundance values (proportions)
  ra.feat = prop.table(feat, 2)
} else {
  ra.feat = feat
}

### apply filters
if (filter.method == 'abundance') {
  # remove features whose abundance is never above the threshold value (e.g. 0.5%) in any of the samples
  f.max = apply(ra.feat, 1, max)
  f.idx = which(f.max >= abundance.cutoff)
} else if (filter.method == 'cum.abundance') {
  # remove features with very low abundance in all samples i.e. ones that are never among the most abundant
  # entities that collectively make up (1-cutoff) of the reads in any sample
  f.idx = vector('numeric', 0)
  # sort features per sample and apply cumsum to identify how many collectively have weight K
  for (s in 1:ncol(ra.feat)) {
    srt = sort(ra.feat[,s], index.return=TRUE)
    cs = cumsum(srt$x)
    m = max(which(cs < abundance.cutoff))
    f.idx = union(f.idx, srt$ix[-(1:m)])
  }
  # an index of those features that collectively make up more than 1-K of the read mass in any sample
  f.idx = sort(f.idx)
} else if (filter.method == 'prevalence') {
  # remove features with low prevalence across samples
  # i.e. ones that are 0 (undetected) in more than (1-cutoff) proportion of samples
  f.idx = which(rowSums(ra.feat > 0) / ncol(ra.feat) > abundance.cutoff)
} else {
  cat('\nunrecognized filter.method, exiting!\n')
  quit(save='no', status = 1)
}

cat('Removed ', nrow(feat)-length(f.idx), ' features whose abundance did not exceed ', abundance.cutoff,
    ' in any sample (retaining ', length(f.idx), ')\n', sep='')
feat = feat[f.idx,]

### postprocessing and output generation
if (rm.unmapped) {
  # remove 'unmapped' feature
  unm.idx = rownames(feat) == 'UNMAPPED' | rownames(feat) == 'unmapped' | rownames(feat) == '-1' | rownames(feat) == 'UNCLASSIFIED' | rownames(feat) == 'unclassified' 
  if (any(unm.idx)) {
    feat = feat[!unm.idx,]
    cat('Removed ', sum(unm.idx), ' features corresponding to UNMAPPED reads',
    ' (retaining ', nrow(feat), ')\n', sep='')
  }
}

# write filtered feature table
write.table(feat, file=fn.out.feat, quote=FALSE, sep='\t', row.names=TRUE, col.names=NA)

cat('\nSuccessfully filtered feature data in ', proc.time()[1] - start.time, ' seconds\n', sep='')
