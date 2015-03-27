# written by Georg Zeller EMBL Heidelberg 2012-2015
# version 0.1.0

### parse commandline arguments
suppressMessages(library('optparse'))
# define arguments
option_list = list(
  make_option('--feat_in', type='character', help='Input file containing features'),
  make_option('--metadata_in', type='character', help='Input file containing metadata'),
  make_option('--feat_out', type='character', help='Output file to which features after selection are written'),
  make_option('--pred_names', type='character', help='names (comma-separated list) of the metavariables to be added to the feature matrix as predictors')
)
# parse arguments
opt = parse_args(OptionParser(option_list=option_list))
fn.in.feat = opt$feat_in
fn.in.meta = opt$metadata_in
fn.out.feat = opt$feat_out
pred.names = opt$pred_names
cat('fn.in.feat =', fn.in.feat, '\n')
cat('fn.in.meta =', fn.in.meta, '\n')
cat('fn.out.feat =', fn.out.feat, '\n')
cat('pred.names =', pred.names, '\n')
cat('\n')

start.time = proc.time()[1]


#### read features and metadata
feat = as.matrix(read.table(file=fn.in.feat, header=TRUE, row.names=1, 
                            sep='\t', stringsAsFactors=FALSE, check.names=FALSE, quote='', comment.char='#'))

meta.data = read.table(file=fn.in.meta, sep='\t', header=TRUE, row.names=1, check.names=FALSE, quote='', comment.char='#')
stopifnot(all(colnames(feat) %in% rownames(meta.data)) && all(rownames(meta.data) %in% colnames(feat)))
m = match(colnames(feat), rownames(meta.data))
meta.data = meta.data[m,]
stopifnot(all(colnames(feat) == rownames(meta.data)))


### add metadata as predictors to teh feature matrix
cnt = 0
if (pred.names != '' && pred.names != 'NULL') {
  pred.names = strsplit(pred.names, ',', fixed=TRUE)[[1]]
  for (p in pred.names) {
    idx = which(colnames(meta.data) == p)
    stopifnot(length(idx) == 1)
    cat('adding ', p, '\n', sep='')
    m = meta.data[,idx]
    if (!all(is.finite(m))) {
      na.cnt = sum(!is.finite(m))
      cat('filling in', na.cnt, 'missing values by mean imputation\n')
      mn = mean(m, na.rm=TRUE)
      m[!is.finite(m)] = mn
    }
    feat = rbind(feat, m)
    rownames(feat)[nrow(feat)] = paste('META-', toupper(p), sep='')
    cnt = cnt + 1
  }
}
cat('added', cnt, 'meta-variables as predictors to the feature matrix\n')

### write combined feature table
write.table(feat, file=fn.out.feat, quote=FALSE, sep='\t', row.names=TRUE, col.names=NA)

cat('\nSuccessfully added metadata to features in ', proc.time()[1] - start.time, ' seconds\n', sep='')
