stargazer <-
function(..., type = "latex", title="", style="default", summary=TRUE, out=NULL, column.labels=NULL, column.separate = NULL, covariate.labels=NULL, dep.var.caption=NULL, dep.var.labels=NULL, dep.var.labels.include=TRUE, align=FALSE, coef=NULL, se=NULL, t=NULL, p=NULL,  t.auto=TRUE, p.auto=TRUE, ci=FALSE, ci.custom=NULL, ci.level=0.95, ci.separator=NULL, apply.coef=NULL, apply.se=NULL, apply.t=NULL, apply.p=NULL, apply.ci=NULL, column.sep.width = "5pt", decimal.mark=NULL, digit.separate=NULL, digit.separator=NULL, digits=NULL, digits.extra=NULL, float=TRUE, float.env="table", font.size=NULL, header=TRUE, initial.zero=NULL, intercept.bottom=TRUE, intercept.top=FALSE, keep=NULL, keep.stat=NULL, label="", model.names=NULL, model.numbers=NULL, no.space=NULL, notes=NULL, notes.align=NULL, notes.append=TRUE, notes.label=NULL, omit=NULL, omit.labels=NULL, omit.stat=NULL, omit.yes.no=c("Yes","No"), order=NULL, ord.intercepts=FALSE, perl=FALSE, rq.se = "nid", selection.equation=FALSE, single.row=FALSE, star.char=NULL, star.cutoffs=NULL, suppress.errors=FALSE, table.placement = "!htbp", zero.component=FALSE, summary.logical=TRUE, nobs=TRUE, mean.sd=TRUE, min.max=TRUE, median=FALSE, iqr=FALSE) {
  save.warn.option <- getOption("warn") 
  options(warn=-1)
  return(.stargazer.wrap(..., type=type, title=title, style=style, summary=summary, out=out,
                         column.labels=column.labels, column.separate = column.separate,
                         covariate.labels=covariate.labels, dep.var.caption = dep.var.caption,
                         dep.var.labels=dep.var.labels, dep.var.labels.include=dep.var.labels.include,
                         align=align, coef=coef, se=se, t=t, p=p, t.auto=t.auto, p.auto=p.auto, 
                         ci=ci, ci.custom=ci.custom, ci.level=ci.level, ci.separator = ci.separator,
                         apply.coef=apply.coef, apply.se=apply.se, apply.t=apply.t, 
                         apply.p=apply.p, apply.ci=apply.ci, 
                         column.sep.width=column.sep.width, decimal.mark=decimal.mark,
                         digit.separate=digit.separate, digit.separator=digit.separator, 
                         digits=digits, digits.extra=digits.extra, 
                         float=float, float.env=float.env,
                         font.size=font.size, header=header, 
                         initial.zero=initial.zero, 
                         intercept.bottom=intercept.bottom, intercept.top=intercept.top, 
                         keep = keep, keep.stat = keep.stat,
                         label = label,
                         model.names=model.names, model.numbers=model.numbers,
                         no.space=no.space, notes=notes, notes.align=notes.align, 
                         notes.append=notes.append, notes.label=notes.label, 
                         omit=omit, omit.labels=omit.labels, omit.stat=omit.stat, omit.yes.no=omit.yes.no, 
                         order=order, ord.intercepts=ord.intercepts, perl=perl, rq.se=rq.se, selection.equation=selection.equation,
                         single.row=single.row, star.char=star.char, 
                         star.cutoffs=star.cutoffs, suppress.errors=suppress.errors,
                         table.placement = table.placement, zero.component=zero.component,
                         summary.logical = summary.logical,
                         nobs=nobs, mean.sd=mean.sd, 
                         min.max=min.max, median=median, iqr=iqr, warn=save.warn.option))
}
