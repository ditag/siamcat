# written by Georg Zeller EMBL Heidelberg 2012-2015
# version 0.1.0

##### auxiliary function to trim whitespace from string
# returns string without leading or trailing whitespace
trim = function (x) {
  gsub("^\\s+|\\s+$", "", x)
}

##### function to parse the header of a label file
### label.header - string in the format: #<TYPE>:<L1>=<class1>;<L2>=<class2>[;<L3>=<class3>]
###   where <TYPE> is a string specifying the type of label variable such as
###   BINARY (for binary classification), CATEGORICAL (for multi-class classification), or CONTINUOUS (for regression)
###   <L1> is a short numeric label for the first class with description <class1> (similarly for the other classes)
parse.label.header = function(label.header) {
  s = strsplit(label.header, ':')[[1]]
  type = trim(s[1])
  if (substr(type, 1, 1) == '#')
    type = trim(substr(type, 2, nchar(type)))
  class.descr = unlist(strsplit(strsplit(trim(s[2]), ';')[[1]], '='))
  l = class.descr[seq(2,length(class.descr),2)]
  class.descr = as.numeric(class.descr[seq(1,length(class.descr)-1,2)])
  names(class.descr) = l

  label.info = list()
  label.info$type = type
  label.info$class.descr = class.descr
  return(label.info)
}

##### function to parse the header of a model file
### TODO documentation
parse.model.header = function(model.header) {
  s = strsplit(model.header, ':')[[1]]
  type = trim(s[1])
  if (substr(type, 1, 1) == '#')
    type = trim(substr(type, 2, nchar(type)))
  label.header = trim(paste(s[2:length(s)], collapse=':'))
  if (substr(label.header, 1, 1) == '[') {
    stopifnot(substr(label.header, nchar(label.header), nchar(label.header)) == ']')
    label.header = substr(label.header, 2, nchar(label.header)-1)
  }
  p = grep('\\(.*\\)', type)
  properties = NULL
  if (length(p) > 0) {
    stopifnot(length(p) == 1)
    stopifnot(substr(type, nchar(type), nchar(type)) == ')')
    properties = substr(type, p+1, nchar(type)-1)
    type = trim(substr(type, 1, p-1))
  }

  model.info = list()
  model.info$type = type
  model.info$properties = properties
  model.info$label.header = label.header
  return(model.info)
}

### TODO docu!
plot.data.range = function(x, y=NULL, x.col='black', y.col='black', labels=NULL, x.suff=NULL, y.suff=NULL) {
  if (is.null(y)) {
    p.m = min(x, na.rm=TRUE)
  } else {
    stopifnot(dim(x)[1] == dim(y)[1])
    p.m = min(c(min(x, na.rm=TRUE), min(y, na.rm=TRUE)))
  }
  plot(rep(p.m, dim(x)[1]), 1:dim(x)[1],
       xlab='', ylab='', yaxs='i', axes=FALSE,
       xlim=c(p.m, 0), ylim=c(0.5, dim(x)[1]+0.5), frame.plot=FALSE, type='n')
  for (v in seq(p.m,-1,1)) {
    abline(v=v, lty=3, col='lightgrey')
  }

  tck = floor(p.m):0
  axis(1, tck, formatC(10^tck, format='E', digits=0), las=1, cex.axis=0.7)

  x.q = apply(x, 1, function (x) quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm=TRUE, names=FALSE))
  if (is.null(y)) {
    # inter-quartile range
    rect(x.q[2,], (1:dim(x)[1])-0.2, x.q[4,], (1:dim(x)[1])+0.2)
    # 90% interval
    segments(x.q[1,], 1:dim(x)[1], x.q[2,], 1:dim(x)[1])
    segments(x.q[4,], 1:dim(x)[1], x.q[5,], 1:dim(x)[1])
    segments(x.q[1,], y0=(1:dim(x)[1])-0.15, y1=(1:dim(x)[1])+0.15)
    segments(x.q[5,], y0=(1:dim(x)[1])-0.15, y1=(1:dim(x)[1])+0.15)
    # median
    segments(x.q[3,], y0=(1:dim(x)[1])-0.2, y1=(1:dim(x)[1])+0.2, lwd=2)
    # scatter plot on top
    for (i in 1:dim(x)[1]) {
      if (nchar(x.col) > 7) {
        # adjust alpha channel by reducing transparency
	a = substr(x.col,nchar(x.col)-1, nchar(x.col))
	a = 1 - (1 - as.numeric(paste('0x', a, sep=''))/255)/2
	x.col = gsub('..$', toupper(as.hexmode(round(a*255))), x.col)
      }
      points(x[i,], rep(i, dim(x)[2])+rnorm(ncol(x),sd=0.05), pch=16, cex=0.6, col=x.col)
    }
  } else {
    y.q = apply(y, 1, function (x) quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm=TRUE, names=FALSE))
    # inter-quartile range
    rect(x.q[2,], 1:dim(x)[1], x.q[4,], (1:dim(x)[1])+0.3, col=x.col)
    rect(y.q[2,], (1:dim(y)[1])-0.3, y.q[4,], 1:dim(y)[1], col=y.col)
    # 90% interval
    segments(x.q[1,], 1:dim(x)[1], x.q[5,], 1:dim(x)[1])#, col=x.col)
    segments(y.q[1,], 1:dim(x)[1], y.q[5,], 1:dim(x)[1])#, col=x.col)
    segments(x.q[1,], y0=1:dim(x)[1], y1=(1:dim(x)[1])+0.2)
    segments(y.q[1,], y0=(1:dim(x)[1])-0.2, y1=1:dim(x)[1])
    segments(x.q[5,], y0=1:dim(x)[1], y1=(1:dim(x)[1])+0.2)
    segments(y.q[5,], y0=(1:dim(x)[1])-0.2, y1=1:dim(x)[1])
    # median
    segments(x.q[3,], y0=1:dim(x)[1], y1=(1:dim(x)[1])+0.3, lwd=3)#, col=x.col)
    segments(y.q[3,], y0=(1:dim(x)[1])-0.3, y1=1:dim(x)[1], lwd=3)#, col=y.col)
    # scatter plot on top
    for (i in 1:dim(x)[1]) {
      if (nchar(x.col) > 7) {
        # adjust alpha channel by reducing transparency
	a = substr(x.col,nchar(x.col)-1, nchar(x.col))
	a = 1 - (1 - as.numeric(paste('0x', a, sep=''))/255)/2
	x.col = gsub('..$', toupper(as.hexmode(round(a*255))), x.col)
      }
      if (nchar(y.col) > 7) {
        # adjust alpha channel by reducing transparency
	a = substr(y.col,nchar(y.col)-1, nchar(y.col))
	a = 1 - (1 - as.numeric(paste('0x', a, sep=''))/255)/2
	y.col = gsub('..$', toupper(as.hexmode(round(a*255))), y.col)
      }
      points(x[i,], rep(i+0.15, ncol(x))+rnorm(ncol(x),sd=0.03), pch=16, cex=0.6, col=x.col)
      points(y[i,], rep(i-0.15, ncol(y))+rnorm(ncol(y),sd=0.03), pch=16, cex=0.6, col=y.col)
    }
  }

  if (!is.null(labels)) {
    stopifnot(length(labels) == dim(c)[1])
    if (!is.null(y) && !is.null(x.suff) && !is.null(y.suff)) {
      for (i in 1:dim(x)[1]) {
        mtext(paste(labels[i], x.suff), 2, line=1, at=i+0.2, las=1, cex=min(0.7, 1-(nrow(x)/70)))
        mtext(y.suff, 2, line=1, at=i-0.2, las=1, cex=min(0.7, 1-(nrow(x)/70)))
      }
    } else {
      for (i in 1:dim(x)[1]) {
        mtext(labels[i], 2, line=1, at=i, las=1, cex=min(0.7, 1-(nrow(x)/50)))
      }
    }
  }
}