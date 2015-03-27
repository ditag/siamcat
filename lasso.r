# written by Georg Zeller EMBL Heidelberg 2012-2015
# version 0.1.0

require('LiblineaR')

liblinear.type = 6
class.weights = c(5, 1)
model.tag = 'LASSO'
eps = 1e-8

##### function to partition training set into cross-validation folds for model selection
assign.fold = function(label, num.folds, stratified) {
  perm = sample(1:length(label), length(label)) / length(label)
  foldid = rep(0, length(label))
  for (f in num.folds:1) {
    foldid[perm <= f/num.folds] = f
  }
  # stratify positive examples
  if (stratified) {
    classes = unique(label)
    for (c in classes) {
      urn = rep(1:num.folds, ceiling(sum(label==c) / num.folds))
      s = sample(urn, sum(label==c))
      foldid[label==c] = s
    }
  }
  stopifnot(length(label) == length(foldid))
  return(foldid)
}

##### function to draw a stratified sample from the label vector
sample.strat = function(label, frac) {
  classes = unique(label)
  s = NULL
  for (c in classes) {
    idx = which(label==c)
    s = c(s, sample(which(label==c), round(frac * sum(label==c))))
  }
  return(s)
}

##### function to train a LASSO model for a single given C
lasso.train = function(x, y, opt.C) {
  # call liblinear training
  names(class.weights) = as.character(sort(unique(y)))
  m = LiblineaR(x, y, type=liblinear.type, cost=opt.C, wi=class.weights, bias=TRUE, epsilon=eps)

  # coefficient vector
  W = m$W
  names(W) = c(colnames(x), 'Bias')

  m$Tag = model.tag
  model = list(W=W, m=m)
  return(model)
}

##### function to apply a trained LASSO model for prediction
lasso.predict = function(x, model) {
  # TODO assert(that classes in model header have correct order!
  stopifnot(dim(model$W)[2] == 1)
  b = model$W[length(model$W)]
  c = as.vector(model$W[-length(model$W)])
  p = vector('numeric', dim(x)[1])
  names(p) = rownames(x)
#  cat('x', dim(x), '\n')
#  cat('c', length(c), '\n')
  for (i in 1:dim(x)[1]) {
    le = b + (as.vector(x[i,]) %*% c)
    p[i] = 1.0 / (1.0 + exp(le))
    # taking negative exponent gives probablities for class being = 2
    # this is equivalent to taking 1 - P("class==1") 
  }
  return(p)
}


##### function to perform model selection in internal cross-validation
lasso.model.selection.cv = function(x, y, foldid, C.vec, PL=NULL, NL=NULL, eval.crit='acc', min.nonzero=1) {
  starttime = proc.time()
  folds = sort(unique(foldid))
  num.folds = length(folds)
  classes = sort(unique(y))
  names(class.weights) = as.character(sort(unique(y)))

  # possible values for eval.crit:
  # 'acc' - accuracy
  # 'auroc' - area under the ROC curve
  # 'auroc.1' - area under the ROC curve on the FPR interval [0, 0.1] (other decimal values are also possible)
  # 'auprc' - area under the precision-recall curve
  # 'auprc.1' - area under the precision-recall curve on the TPR interval [0, 0.1] (other decimal values possible)
 
  if (is.null(PL)) {
    stopifnot(!is.null(NL))
    cat('No label information provided: accuracy will be used as model selection criterion.\n')
    eval.crit = 'acc' 
  }

#  cat('Running model selection grid search in ', num.folds, '-fold cross-validation\n', sep='')
  # performance vector
  perf = vector('numeric', length(C.vec))
  nonzero.coeff = matrix(0, num.folds, length(C.vec))
  cnt = 1
  for (c in C.vec) {
    cat('  C=', c, '\n', sep='')
    pred = vector('numeric', length(y))
    for (f in 1:num.folds) {
      test.idx = which(foldid == folds[f])
      train.idx = which(foldid != folds[f])
      stopifnot(length(intersect(train.idx, test.idx)) == 0)
      trainx = x[train.idx,]
      trainy = y[train.idx]
      # sort labels to make sure that all models contain classes in the same order
      idx = sort(trainy, index.return=TRUE)$ix
      trainy = trainy[idx]
      trainx = trainx[idx,]
      stopifnot(all(unique(trainy) == classes))
      testx  = x[test.idx,]
      
      model = LiblineaR(trainx, trainy, type=liblinear.type, cost=c, wi=class.weights, bias=TRUE, epsilon=eps)
#      cat(model$W, '\n')
#      cat(dim(trainx), '\n')
#      cat('nonzero:', sum(model$W[1:(ncol(trainx)-1)] != 0), '\n')
      nonzero.coeff[f,cnt] = sum(model$W[1:(ncol(trainx)-1)] != 0)
      if (eval.crit=='acc') {
        p = predict(model, testx, proba=FALSE)
        pred[test.idx] = as.vector(p$predictions)
      } else {
        p = predict(model, testx, proba=TRUE)
        pred[test.idx] = as.vector(p$probabilities[,as.character(PL)])
        stopifnot(length(test.idx) == length(p$probabilities[,as.character(PL)]))
      }
    }
    et = proc.time() - starttime
    et = et['elapsed']

    if (eval.crit=='acc') {
      # TODO doesn't work for continous predictions (like posterior prob)!
      acc = mean(pred == y)
      perf[cnt] = acc
      cat('  (ACCURACY=')
    } else if (substr(eval.crit, 1, nchar('auprc')) == 'auprc') {
      s = strsplit(eval.crit, '.', fixed=TRUE)
      max.tpr = ifelse(length(s[[1]]) > 1, as.numeric(s[[1]][2])/10, 1.0)
      ev = eval.classifier(pred, y)
      aupr = calc.aupr(ev, max.tpr)
      perf[cnt] = aupr
      cat('  (AU-PRC[0,', max.tpr, ']=', sep='')
    } else if (substr(eval.crit, 1, nchar('auroc')) == 'auroc') {
      s = strsplit(eval.crit, '.', fixed=TRUE)
      max.fpr = ifelse(length(s[[1]]) > 1, as.numeric(s[[1]][2])/10, 1.0)
      #ev = eval.classifier(pred, y)
      #auroc = calc.auroc(ev, max.fpr)
      #perf[cnt] = auroc
      ev.roc = roc(response=y, predictor=pred, sp=seq(0, max.fpr, 0.025), ci=FALSE)
      perf[cnt] = ev.roc$auc
      cat('  (AU-ROC[0,', max.fpr, ']=', sep='')
    } else {
      stop('unrecognized eval.crit argument:', eval.crit)
    }
    cat(format(perf[cnt], digits=3), ' [', paste(nonzero.coeff[,cnt], collapse=','), '], elapsed time: ', et, ')\n', sep='')
    cnt = cnt+1
  }
  # check whether these models have sufficiently many nonzero coefficients 
  suff.nonzero = apply(nonzero.coeff, 2, min) > min.nonzero
  stopifnot(sum(suff.nonzero) > 0) # otherwise ALL models will be discarded
  # TODO better error handling!
  #cat('removing',  'models with insufficient non-zero coefficients...\n')
  #cat(apply(nonzero.coeff, 2, min), '\n')
  nonzero.coeff = nonzero.coeff[,suff.nonzero]
  C.vec = C.vec[suff.nonzero]
  perf = perf[suff.nonzero]
  opt.idx = which.max(perf)
  opt.C = C.vec[opt.idx]
  et = proc.time() - starttime
  et = et['elapsed']
  if (eval.crit=='acc') {
    cat('  ACCURACY=',)
  } else if (substr(eval.crit, 1, nchar('auprc')) == 'auprc') {
    s = strsplit(eval.crit, '.', fixed=TRUE)
    max.tpr = ifelse(length(s[[1]]) > 1, as.numeric(s[[1]][2])/10, 1.0)
    if (max.tpr < 1.0)
      cat('  AU-PRC[0,', max.tpr, ']=', sep='')
    else
      cat('  AU-PRC=')
  } else if (substr(eval.crit, 1, nchar('auroc')) == 'auroc') {
    s = strsplit(eval.crit, '.', fixed=TRUE)
    max.fpr = ifelse(length(s[[1]]) > 1, as.numeric(s[[1]][2])/10, 1.0)
    if (max.fpr < 1.0)
      cat('  AU-ROC[0,', max.fpr, ']=', sep='')
    else
      cat('  AU-ROC=')
  } else {
    stop('unrecognized eval.crit argument:', eval.crit)
  }
  cat(max(perf), ' (optimal C=', opt.C, ' [', paste(nonzero.coeff[,opt.idx], collapse=','), '], elapsed time: ', et, ')\n', sep='')
  return(opt.C)
}
##### end function
