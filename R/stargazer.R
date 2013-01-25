stargazer <-
function(..., title="", style="default", summary=TRUE, covariate.labels=NULL, dep.var.labels=NULL, align=FALSE, coef=NULL, se=NULL, t=NULL, p=NULL,  t.auto=TRUE, p.auto=TRUE, decimal.mark=NULL, digit.separate=NULL, digit.separator=NULL, digits=NULL, digits.extra=NULL, initial.zero=NULL, intercept.top=NULL, label="", model.names=NULL, model.numbers=NULL, notes=NULL, notes.align=NULL, notes.label=NULL, omit=NULL, omit.labels=NULL, omit.stat=NULL, omit.yes.no=c("Yes","No"), ord.intercepts=FALSE, star.char=NULL, star.cutoffs=NULL, zero.component=FALSE, summary.logical=TRUE, nobs=TRUE, mean.sd=TRUE, min.max=TRUE, median=FALSE, iqr=FALSE) {
  return(.stargazer.wrap(..., title=title, style=style, summary=summary, covariate.labels=covariate.labels, 
                         dep.var.labels=dep.var.labels, align=align, coef=coef, se=se, t=t, p=p, t.auto=t.auto, p.auto=p.auto, decimal.mark=decimal.mark,
                         digit.separate=digit.separate, digit.separator=digit.separator, 
                         digits=digits, digits.extra=digits.extra, initial.zero=initial.zero, 
                         intercept.top=intercept.top, label = label,
                         model.names=model.names, model.numbers=model.numbers,
                         notes=notes, notes.align=notes.align, notes.label=notes.label, 
                         omit=omit, omit.labels=omit.labels, omit.stat=omit.stat, omit.yes.no=omit.yes.no, 
                         ord.intercepts=ord.intercepts, star.char=star.char, 
                         star.cutoffs=star.cutoffs, zero.component=zero.component,
                         summary.logical = summary.logical,
                         nobs=nobs, mean.sd=mean.sd, 
                         min.max=min.max, median=median, iqr=iqr))
}