# written by Georg Zeller EMBL Heidelberg 2012-2015
# version 0.1.0

### parse commandline arguments
suppressMessages(library('optparse'))
# define arguments
option_list = list(
  make_option(c('-s', '--srcdir'), type='character', help='Source directory of this and other utility scripts'),
  make_option('--metadata_in', type='character', help='Input file containing meta-data'),
  make_option('--metadata_out', type='character', help='Output file to which meta-data after selection is written'),
  make_option('--label_in', type='character', help='Input file containing labels'),
  make_option('--label_out', type='character', help='Output file to which labels after selection are written'),
  make_option('--feat_in', type='character', help='Input file containing features'),
  make_option('--feat_out', type='character', help='Output file to which features after selection are written'),
  make_option('--filter_var', type='character', help='Meta-data variable on which filtering is based'),
  make_option('--allowed_range', type='character', default=NULL, help='Range of allowed values (closed interval)'),
  make_option('--allowed_set', type='character', default=NULL, help='Set of allowed values (comma-separated list)')
)
# parse arguments
opt = parse_args(OptionParser(option_list=option_list))
source.dir = opt$srcdir
fn.in.meta = opt$metadata_in
fn.out.meta = opt$metadata_out
fn.in.label = opt$label_in
fn.out.label = opt$label_out
fn.in.feat = opt$feat_in
fn.out.feat = opt$feat_out
filter = opt$filter_var
allowed.range = opt$allowed_range
# parse interval
if(!is.null(allowed.range)) {
  allowed.range = gsub('\\[|\\]','', allowed.range)
  allowed.range = as.numeric(unlist(strsplit(allowed.range,',')))
  stopifnot(length(allowed.range) == 2)
  stopifnot(allowed.range[1] <= allowed.range[2])
}
allowed.set = opt$allowed_set
if (!is.null(allowed.set)) {
  # parse set
  allowed.set = gsub('\\{|\\}','', allowed.set)
  allowed.set = as.numeric(unlist(strsplit(allowed.set,',')))
  allowed.set = sort(unique(allowed.set))
}
cat('source.dir =', source.dir, '\n')
cat('fn.in.meta =', fn.in.meta, '\n')
cat('fn.out.meta =', fn.out.meta, '\n')
cat('fn.in.label =', fn.in.label, '\n')
cat('fn.out.label =', fn.out.label, '\n')
cat('fn.in.feat =', fn.in.feat, '\n')
cat('fn.out.feat =', fn.out.feat, '\n')
cat('filter =', filter, '\n')
if (!xor(is.null(allowed.range), is.null(allowed.set))) {
  cat('Neither allowed.range nor allowed.set (or both at the same time) have been provided, exiting!\n')
  quit(save='no', status = 1)
} else {
  if (!is.null(allowed.range)) {
  cat('allowed.range = [', paste(allowed.range, collaps=','), ']\n', sep='')
  } else {
    cat('allowed.set = {', paste(allowed.set, collapse=','), '}\n', sep='')
  }
}
cat('\n')
if (substr(source.dir, nchar(source.dir), nchar(source.dir)) != '/') {
  source.dir = paste(source.dir, '/', sep='')
}

start.time = proc.time()[1]


### imports
source(paste(source.dir, 'utils.r', sep=''))


### read label, feature and meta- data
# features
feat = read.table(file=fn.in.feat, sep='\t', header=TRUE, row.names=1, stringsAsFactors=FALSE, check.names=FALSE, quote='', comment.char='#')
# labels (assuming the label file has 1 column)
label = read.table(file=fn.in.label, sep='\t', header=FALSE, row.names=1, check.names=FALSE, quote='', comment.char='#')
n = rownames(label)
stopifnot(n == colnames(feat))
label = as.numeric(label[,1])
names(label) = n
con = file(fn.in.label, 'rt') 
label.header = readLines(con, 1)
close(con)
# meta-data
meta.data = read.table(file=fn.in.meta, sep='\t', header=TRUE, row.names=1, check.names=FALSE, quote='', comment.char='#')
stopifnot(all(names(label) %in% rownames(meta.data)) && all(rownames(meta.data) %in% names(label)))
m = match(names(label), rownames(meta.data))
meta.data = meta.data[m,]
stopifnot(all(names(label) == rownames(meta.data)))


### select samples fulfilling the filter criteria
# (i.e. samples having certain metadata values)
stopifnot(filter %in% colnames(meta.data))
filter.var = meta.data[,filter]
if (!is.null(allowed.range)) {
  s.idx = !is.na(filter.var) & filter.var >= allowed.range[1] & filter.var <= allowed.range[2]
  cat('Removed ', sum(!s.idx), ' samples with ', filter,
      ' not in [', paste(allowed.range, collapse=', '), '] (retaining ', sum(s.idx), ')\n', sep='')
} else {
  s.idx = !is.na(filter.var) & filter.var %in% allowed.set
  cat('Removed ', sum(!s.idx), ' samples with ', filter,
      ' not in {', paste(allowed.set, collapse=', '), '} (retaining ', sum(s.idx), ')\n', sep='')
}
feat = feat[,s.idx]
label = label[s.idx]
meta.data = meta.data[s.idx,]


### write label, feature and meta-data with selected sample set 
# labels
write(label.header, file=fn.out.label, append=FALSE)
write.table(label, file=fn.out.label, quote=FALSE, sep='\t', row.names=TRUE, col.names=FALSE, append=TRUE)
# features
write.table(feat, file=fn.out.feat, quote=FALSE, sep='\t', row.names=TRUE, col.names=NA)
# meta-data
write.table(meta.data, file=fn.out.meta, quote=FALSE, sep='\t', row.names=TRUE, col.names=NA)


cat('\nSuccessfully selected samples in ', proc.time()[1] - start.time, ' seconds\n', sep='')
