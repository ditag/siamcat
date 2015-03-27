# written by Georg Zeller EMBL Heidelberg 2012-2015
# version 0.1.0

# evaluates the predictive performance of a classifier on a labeled data sets
# returns a list with vectors containing TP, FP, TN, FN for each threshold value on the predictions
# (where TP = true positives, FP = false positives, TN = true negatives, FN = false negatives)
eval.classifier = function(predictions, labels) {
  stopifnot(dim(labels) == NULL)
  stopifnot(length(unique(labels)) == 2)
  stopifnot(all(is.finite(predictions)))
  PL = max(labels)
  NL = min(labels)
  # calculate thresholds, one between each subsequent pair of sorted prediction values
  # this is ignorant to whether predictions is in matrix or vector format (see below)
  t = predictions
  dim(t) = NULL
  t = sort(unique(t))
  t = rev(c(min(t)-1, (t[-1]+t[-length(t)])/2, max(t)+1))
  if (is.null(dim(predictions))) {
    # assuming that a single model was applied to predict the data set
    stopifnot(length(labels) == length(predictions))
    # actual evaluations per threshold value
    tp = vector('numeric', length(t))
    fp = vector('numeric', length(t))
    tn = vector('numeric', length(t))
    fn = vector('numeric', length(t))
    for (i in 1:length(t)) {
      tp[i] = sum(labels==PL & predictions>t[i])
      fp[i] = sum(labels==NL & predictions>t[i])
      tn[i] = sum(labels==NL & predictions<t[i])
      fn[i] = sum(labels==PL & predictions<t[i])
    }
  } else {
    # assuming that several models were applied to predict the same data and predictions of each model occupy one column
    stopifnot(length(labels) == nrow(predictions))
    tp = matrix(0, nrow=length(t), ncol=ncol(predictions))
    fp = matrix(0, nrow=length(t), ncol=ncol(predictions))
    tn = matrix(0, nrow=length(t), ncol=ncol(predictions))
    fn = matrix(0, nrow=length(t), ncol=ncol(predictions))
    for (c in 1:ncol(predictions)) {
      for (r in 1:length(t)) {
        tp[r,c] = sum(labels==PL & predictions[,c] > t[r])
        fp[r,c] = sum(labels==NL & predictions[,c] > t[r])
        tn[r,c] = sum(labels==NL & predictions[,c] < t[r])
        fn[r,c] = sum(labels==PL & predictions[,c] < t[r])
      }
    }
  }
  return(list(tp=tp, tn=tn, fp=fp, fn=fn, thresholds=t))
}

# returns a vector of x and y values for plotting a receiver operating characteristic curve
get.roc = function(eval) {
  if (is.null(dim(eval$tp))) {
    stopifnot(!is.null(eval$fp))
    stopifnot(!is.null(eval$tn))
    stopifnot(!is.null(eval$fn))
    fpr = eval$fp / (eval$tn + eval$fp)
    tpr = eval$tp / (eval$tp + eval$fn)  
  } else {
    stopifnot(ncol(eval$tp) == ncol(eval$fp))
    stopifnot(ncol(eval$tp) == ncol(eval$tn))
    stopifnot(ncol(eval$tp) == ncol(eval$fn))
    fpr = matrix(NA, nrow=nrow(eval$tp), ncol=ncol(eval$tp))
    tpr = matrix(NA, nrow=nrow(eval$tp), ncol=ncol(eval$tp))
    for (c in 1:ncol(eval$tp)) {
      fpr[,c] = eval$fp[,c] / (eval$tn[,c] + eval$fp[,c])
      tpr[,c] = eval$tp[,c] / (eval$tp[,c] + eval$fn[,c])
    }
  }
  return(list(x=fpr, y=tpr))
}

# plots a receiver operating characteristic curve
plot.roc.curve = function(eval) {
  roc = get.roc(eval)
  plot(roc$x, roc$y,
       xlim=c(0,1), ylim=c(0,1), type='l', # type='o' pch=16, cex=0.4,
       xlab='False positive rate', ylab='True positive rate')
}

# returns a vector of x and y values for plotting a precision-recall curve
get.pr = function(eval) {
  tpr = eval$tp / (eval$tp + eval$fn)
  ppv = eval$tp / (eval$tp + eval$fp)
  # at thresholds where the classifier makes no positive predictions at all,
  # we (somewhat arbitrarily) set its precision to 1
  ppv[is.na(ppv)] = 1.0
  return(list(x=tpr, y=ppv))
}

# plots a precision-recall curve
plot.pr.curve = function(eval) {
  pr = get.pr(eval)
  plot(pr$x, pr$y,
       xlim=c(0,1), ylim=c(0,1), type='l', # type='o' pch=16, cex=0.4,
       xlab='Recall (TPR)', ylab='Precision (PPV)')
}


# calculates the area under a curve using a trapezoid approximation
area.trapez = function(x, y) {
  if (x[1] > x[length(x)]) {
    x = rev(x)
    y = rev(y)
  }
  xd = x[-1] - x[-length(x)]
  ym = 0.5 * (y[-1] + y[-length(y)])
  return(xd %*% ym)
}

# calculates the area under the ROC curve (over the interval [0, max.fpr], if specified)
calc.auroc = function(eval, max.fpr=1.0) {
  roc = get.roc(eval)
  idx = roc$x <= max.fpr
  return(area.trapez(roc$x[idx], roc$y[idx]))
}

# calculates the area under the precision-recall curve (over the interval [0, max.tpr], if specified)
calc.aupr = function(eval, max.tpr=1.0) {
  pr = get.pr(eval)
  idx = pr$x <= max.tpr
  return(area.trapez(pr$x[idx], pr$y[idx]))
}

# returns the positions of local minima in a given (ordered) vector
# local minima are defined as having ext adjacent smaller values in both directions
local.min = function(y, ext) {
  m = list()
  for (i in 1:length(y)) {
    env = c((i-ext):(i-1), (i+1):(i+ext))
    env = env[env > 0 & env <= length(y)]
    if (y[i] > max(y[env])) {
      m = c(m,i)
    }
  }
  return(unlist(m))
}
