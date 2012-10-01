stargazer <-
function(..., title="", style="default", covariate.labels=NULL, dep.var.labels=NULL, decimal.mark=NULL, digit.separate=NULL, digit.separator=NULL, digits=NULL, digits.extra=NULL, initial.zero=NULL, intercept.top=NULL, model.names=NULL, model.numbers=NULL, notes=NULL, notes.align=NULL, notes.label=NULL, omit=NULL, omit.labels=NULL, omit.yes.no=c("Yes","No"), ord.intercepts=FALSE, star.char=NULL, star.cutoffs=NULL, nobs=TRUE, mean.sd=TRUE, min.max=TRUE, median=FALSE, iqr=FALSE) {
  return(.stargazer.wrap(..., title=title, style=style, covariate.labels=covariate.labels, 
                         dep.var.labels=dep.var.labels, decimal.mark=decimal.mark, 
                         digit.separate=digit.separate, digit.separator=digit.separator, 
                         digits=digits, digits.extra=digits.extra, initial.zero=initial.zero, 
                         intercept.top=intercept.top, model.names=model.names, model.numbers=model.numbers,
                         notes=notes, notes.align=notes.align, notes.label=notes.label, 
                         omit=omit, omit.labels=omit.labels, omit.yes.no=omit.yes.no, 
                         ord.intercepts=ord.intercepts, star.char=star.char, 
                         star.cutoffs=star.cutoffs, nobs=nobs, mean.sd=mean.sd, 
                         min.max=min.max, median=median, iqr=iqr))
}

