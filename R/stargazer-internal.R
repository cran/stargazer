.onAttach <- 
function(libname, pkgname) {
  packageStartupMessage("\nPlease cite as: \n")
  packageStartupMessage(" Hlavac, Marek (2013). stargazer: LaTeX code for well-formatted regression and summary statistics tables.")
  packageStartupMessage(" R package version 3.0. http://CRAN.R-project.org/package=stargazer \n")
}

.stargazer.wrap <-
  function(...,title, style, summary, covariate.labels, dep.var.labels, align, coef, se, t, p, t.auto, p.auto, decimal.mark, digit.separate, digit.separator, digits, digits.extra, initial.zero, intercept.top, label, model.names, model.numbers, notes, notes.align, notes.label, omit, omit.labels, omit.stat, omit.yes.no, ord.intercepts, star.char, star.cutoffs, zero.component, summary.logical, nobs, mean.sd, min.max, median, iqr) {
    
  .add.model <-
  function(object.name, user.coef=NULL, user.se=NULL, user.t=NULL, user.p=NULL, auto.t=TRUE, auto.p=TRUE) {
    
    .summary.object <<- summary(object.name)

  	.global.models <<- append(.global.models, .model.identify(object.name))
  	
  	.global.dependent.variables <<- append(.global.dependent.variables, .dependent.variable(object.name))
  	.global.dependent.variables.written <<- append(.global.dependent.variables.written, .dependent.variable.written(object.name))

  	.global.N <<- append(.global.N, .number.observations(object.name))
  	.global.LL <<- append(.global.LL, .log.likelihood(object.name))
  	.global.R2 <<- append(.global.R2, .r.squared(object.name))
  	.global.max.R2 <<- append(.global.max.R2, .max.r.squared(object.name))
  	.global.adj.R2 <<- append(.global.adj.R2, .adj.r.squared(object.name))
  	.global.AIC <<- append(.global.AIC, .AIC(object.name))
    .global.BIC <<- append(.global.BIC, .BIC(object.name))
  	.global.scale <<- append(.global.scale, .get.scale(object.name))
    .global.UBRE <<- append(.global.UBRE, .gcv.UBRE(object.name))
  	.global.sigma2 <<- append(.global.sigma2, .get.sigma2(object.name))
  	
    .global.theta <<- cbind(.global.theta, .get.theta(object.name))
  	.global.SER <<- cbind(.global.SER, .SER(object.name))
  	.global.F.stat <<- cbind(.global.F.stat, .F.stat(object.name))
  	.global.chi.stat <<- cbind(.global.chi.stat, .chi.stat(object.name))
  	.global.wald.stat <<- cbind(.global.wald.stat, .wald.stat(object.name))
  	.global.lr.stat <<- cbind(.global.lr.stat, .lr.stat(object.name))
  	.global.logrank.stat <<- cbind(.global.logrank.stat, .logrank.stat(object.name))
  	.global.null.deviance <<- cbind(.global.null.deviance, .null.deviance(object.name))
  	.global.residual.deviance <<- cbind(.global.residual.deviance, .residual.deviance(object.name))

  	max.length <- length(.global.coefficient.variables)+length(.coefficient.variables(object.name))

  	# add RHS variables and coefficients
  	coef.var <- .coefficient.variables(object.name)

  	temp.gcv <- rep(NA,each=1,times=max.length)

  	temp.gcv[1:length(.global.coefficient.variables)] <- .global.coefficient.variables

  	how.many.gcv <- length(.global.coefficient.variables)

  	position <- 0
  	for (i in seq(1:length(coef.var))) {

  		found <- FALSE
  		for (j in seq(1:length(.global.coefficient.variables))) {
  				if (coef.var[i] == .global.coefficient.variables[j]) {
  					found <- TRUE
  					for (k in 1:how.many.gcv) {
  						if (coef.var[i]==temp.gcv[k]) {
  							position <- k
  						}
  					}
  				}
  		}

  		# If variable was found, no need to add it
  		if (found == FALSE) {
   
        # append new variable to list of regressors
        while ((position < how.many.gcv) & (!(temp.gcv[position+1] %in% coef.var))) {
          position <- position + 1
        }
        
  			temp.gcv <- append(temp.gcv, coef.var[i], after=position)
  			how.many.gcv <- how.many.gcv + 1
  			position <- position + 1
  		}
		
  	}
  
  	.global.coefficient.variables <<- temp.gcv[1:how.many.gcv]

  	# build up coefficients from scratch
  	temp.coefficients <- temp.std.errors <- temp.t.stats <- temp.p.values <- matrix(data = NA, nrow = length(.global.coefficient.variables), ncol = ncol(.global.coefficients)+1)
	  rownames(temp.coefficients) <- rownames(temp.std.errors) <- rownames(temp.t.stats) <- rownames(temp.p.values) <- .global.coefficient.variables

  	# fill in from previous iteration of .global coefficients
  	for (row in .global.coefficient.variables) {
		  for (col in seq(1:ncol(.global.coefficients))) {
  			if (sum(as.vector(names(.global.coefficients[,col])==row))!=0) { 
  				if (!is.null(.global.coefficients)) { temp.coefficients[row, col] <- .global.coefficients[row, col] }
  				if (!is.null(.global.std.errors)) { temp.std.errors[row, col] <- .global.std.errors[row, col] }
  				if (!is.null(.global.t.stats)) { temp.t.stats[row, col] <- .global.t.stats[row, col] }
  				if (!is.null(.global.p.values)) { temp.p.values[row, col] <- .global.p.values[row, col] }
  			}
  		}
      
      feed.coef <- NA
      # coefficients and standard errors
  		if (!is.null(.get.coefficients(object.name, user.coef)[row])) { temp.coefficients[row, ncol(temp.coefficients)] <- feed.coef <- .get.coefficients(object.name, user.coef)[row] }
  		if (!is.null(.get.standard.errors(object.name, user.se)[row])) { temp.std.errors[row, ncol(temp.std.errors)] <- .get.standard.errors(object.name, user.se)[row] }
      
      # t-stats and p-values
      if (!is.null(user.coef)) { feed.coef <- user.coef }   # feed user-defined coefficients, if available
      if (!is.null(.get.t.stats(object.name, user.t, auto.t, feed.coef, user.se)[row])) { temp.t.stats[row, ncol(temp.std.errors)] <- .get.t.stats(object.name, user.t, auto.t, feed.coef, user.se)[row] }
  		if (!is.null(.get.p.values(object.name, user.p, auto.p, feed.coef, user.se)[row])) { temp.p.values[row, ncol(temp.std.errors)] <- .get.p.values(object.name, user.p, auto.p, feed.coef, user.se)[row] }
  	}

  	if (!is.null(temp.coefficients)) { .global.coefficients <<- temp.coefficients }
  	if (!is.null(temp.std.errors)) { .global.std.errors <<- temp.std.errors }
  	if (!is.null(temp.t.stats)) { .global.t.stats <<- temp.t.stats }
  	if (!is.null(temp.p.values)) { .global.p.values <<- temp.p.values }

  }

  .adj.r.squared <-
  function(object.name) {

  	model.name <- .get.model.name(object.name)

  	if (!(model.name %in% c("arima", "lmer", "glmer", "nlmer"))) {
  		if (!is.null(suppressMessages(.summary.object$adj.r.squared))) {
  			return(as.vector(suppressMessages(.summary.object$adj.r.squared)))
  		}
  		else if (model.name %in% c("normal.gam", "logit.gam", "probit.gam", "poisson.gam", "gam()")) {
  			return(as.vector(.summary.object$r.sq))
  		}
      else if (model.name %in% c("plm")) {
        return(as.vector(.summary.object$r.squared["adjrsq"]))
      }
  	}
  	return(NA)
  }

  .adjust.settings.style <-
  function(what.style) {
    style <- tolower(what.style)
  
    if (style == "all") {
      .format.table.parts <<- c("=!","dependent variable label","dependent variables","models","numbers","-","coefficients","-","omit","-","additional","N","R-squared","adjusted R-squared","max R-squared","log likelihood","sigma2","theta(se)*(t)", "SER(df)","F statistic(df)*(p)","chi2(df)*(p)","Wald(df)*(p)","LR(df)*(p)","logrank(df)*(p)","AIC","BIC","UBRE","residual deviance","null deviance","=!","notes")  
      .format.coefficient.table.parts <<- c("variable name","coefficient*","standard error","t-stat","p-value")  
    }
  
    else if (style == "all2") {
      .format.table.parts <<- c("=!","dependent variable label","dependent variables","models","numbers","-","coefficients","-","omit","-","additional","N","R-squared","adjusted R-squared","max R-squared","log likelihood","sigma2","theta(se)*(t)", "SER(df)","F statistic(df)*(p)","chi2(df)*(p)","Wald(df)*(p)","LR(df)*(p)","logrank(df)*(p)","AIC","BIC","UBRE","residual deviance","null deviance","=!","notes")  
      .format.coefficient.table.parts <<- c("variable name","coefficient*","standard error")  
    }
  
    # aer = American Economic Review
    else if (style == "aer") {
      .format.table.parts <<- c("=!","dependent variables","models","numbers","-","coefficients","omit","additional","N","R-squared","adjusted R-squared","max R-squared","log likelihood","theta(se)*", "AIC","BIC","UBRE", "SER(df)","F statistic(df)*","chi2(df)*","Wald(df)*","LR(df)*","logrank(df)*","-!","notes")
      .format.models.skip.if.one <<- TRUE
    
      .format.until.nonzero.digit <<- FALSE
      .format.max.extra.digits <<- 0    
    
      .format.model.left <<- ""
      .format.model.right <<- ""
    
      .format.note <<- "\\textit{Notes:}"
      .format.note.alignment <<- "l"
      .format.note.content <<- c("$^{***}$Significant at the [***] percent level.","$^{**}$Significant at the [**] percent level.","$^{*}$Significant at the [*] percent level.")
    }
  
    # ajps = American Journal of Political Science
    else if (style == "ajps") {
      .format.table.parts <<- c("-!","dependent variables","models","numbers","-","coefficients","omit","additional","N","R-squared","adjusted R-squared","max R-squared","log likelihood","theta(se)*", "SER(df)","F statistic(df)*","chi2(df)*","Wald(df)*","LR(df)*","logrank(df)*","AIC","BIC","UBRE","-!","notes")    
      .format.models.skip.if.one <<- TRUE
      .format.digit.separator <<- ""
      .format.dependent.variables.left <<- "\\textbf{"
      .format.dependent.variables.right <<- "}"
      .format.models.left <<- "\\textbf{"
      .format.models.right <<- "}"
      .format.numbers.left <<- "\\textbf{Model "
      .format.numbers.right <<- "}"
      .format.coefficient.table.parts <<- c("variable name","coefficient*","standard error") 
      .format.N <<- "N"
      .format.AIC <<- "AIC"
      .format.BIC <<- "BIC"
      .format.chi.stat <<- "Chi-square"
      .format.R2 <<- "R-squared"
      .format.adj.R2 <<- "Adj. R-squared"
      .format.max.R2 <<- "Max. R-squared"
      .format.note <<- ""
      .format.note.content <<- c("$^{***}$p $<$ [.***]; $^{**}$p $<$ [.**]; $^{*}$p $<$ [.*]")
      .format.note.alignment <<- "l"
      .format.s.stat.parts <<- c("-!","stat names","-!","statistics1","-!")
    }  
  
    # ajs = American Journal of Sociology
    else if (style == "ajs") {
      .format.table.parts <<- c(" ","=!","dependent variables","models","numbers","-","coefficients","omit","additional","N","R-squared","adjusted R-squared","max R-squared","log likelihood","theta(se)*", "AIC","BIC","UBRE", "SER(df)","F statistic(df)*","chi2(df)*","Wald(df)*","LR(df)*","logrank(df)*","-!","notes")
      .format.models.skip.if.one <<- TRUE
      .format.dependent.variables.capitalize <<- TRUE
    
      .format.numbers.left <<- ""
      .format.numbers.right <<- ""
    
      .format.until.nonzero.digit <<- FALSE
      .format.max.extra.digits <<- 0    
    
      .format.model.left <<- ""
      .format.model.right <<- ""
    
      .format.note <<- "\\textit{Notes:}"
      .format.note.alignment <<- "l"
      .format.note.content <<- c("$^{*}$P $<$ [.*]","$^{**}$P $<$ [.**]","$^{***}$P $<$ [.***]")
      .format.one.star <<- 0.05
      .format.two.stars <<- 0.01
      .format.three.stars <<- 0.001
    
      .format.initial.zero <<- FALSE
    }
  
    # apsr = American Political Science Review
    else if (style == "apsr") {
      .format.table.parts <<- c("-!","dependent variables","models","numbers","-","coefficients","omit","additional","N","R-squared","adjusted R-squared","max R-squared","log likelihood","theta(se)*", "SER(df)","F statistic(df)*","chi2(df)*","Wald(df)*","LR(df)*","logrank(df)*","AIC","BIC","UBRE","-!","notes")    
      .format.models.skip.if.one <<- TRUE
      .format.models.left <<- ""
      .format.models.right <<- ""
      .format.coefficient.table.parts <<- c("variable name","coefficient*","standard error")
      .format.N <<- "N"
      .format.AIC <<- "AIC"
      .format.BIC <<- "BIC"
      .format.chi.stat <<- "chi$^{2}$"
      .format.note <<- ""
      .format.note.content <<- c("$^{*}$p $<$ [.*]; $^{**}$p $<$ [.**]; $^{***}$p $<$ [.***]")
      .format.note.alignment <<- "l"
      .format.s.stat.parts <<- c("-!","stat names","-!","statistics1","-!")
    }
    
    # asq = Administrative Science Quarterly
    else if (style == "asq") {
      .format.table.parts <<- c("-!","dependent variables","models","numbers","-","coefficients","omit","additional","N","R-squared","adjusted R-squared","max R-squared","log likelihood","theta(se)*", "SER(df)","F statistic(df)*","chi2(df)*","Wald(df)*","LR(df)*","logrank(df)*","AIC","BIC","UBRE","-!","notes")    
      .format.models.skip.if.one <<- TRUE
      .format.digit.separator <<- ""
      .format.dependent.variables.left <<- "\\textbf{"
      .format.dependent.variables.right <<- "}"
      .format.models.left <<- "\\textbf{"
      .format.models.right <<- "}"
      .format.numbers.left <<- "\\textbf{Model "
      .format.numbers.right <<- "}"
      .format.coefficient.table.parts <<- c("variable name","coefficient*","standard error") 
      .format.AIC <<- "AIC"
      .format.BIC <<- "BIC"
      .format.chi.stat <<- "Chi-square"
      .format.R2 <<- "R-squared"
      .format.adj.R2 <<- "Adj. R-squared"
      .format.max.R2 <<- "Max. R-squared"
      .format.note <<- ""
      .format.note.content <<- c("$^{\\bullet}$p $<$ [.*]; $^{\\bullet\\bullet}$p $<$ [.**]; $^{\\bullet\\bullet\\bullet}$p $<$ [.***]")
      .format.note.alignment <<- "l"
      .format.s.stat.parts <<- c("-!","stat names","-!","statistics1","-!")
      .format.significance.star <<- "\\bullet"
    }  
  
    # asr = American Sociological Review
    else if (style == "asr") {
      .format.table.parts <<- c("-!","dependent variables","models","numbers","-","coefficients","omit","additional","N","R-squared","adjusted R-squared","max R-squared","log likelihood","theta(se)", "SER(df)","F statistic(df)*","chi2(df)*","Wald(df)*","LR(df)*","logrank(df)*","AIC","BIC","UBRE","-!","notes")    
      .format.models.skip.if.one <<- TRUE
      .format.models.left <<- ""
      .format.models.right <<- ""
      .format.coefficient.table.parts <<- c("variable name","coefficient*")
      .format.N <<- "\\textit{N}"
      .format.AIC <<- "AIC"
      .format.BIC <<- "BIC"
      .format.chi.stat <<- "chi$^{2}$"
      .format.note <<- ""
      .format.note.content <<- c("$^{*}$p $<$ [.*]; $^{**}$p $<$ [.**]; $^{***}$p $<$ [.***]")
      .format.one.star <<- 0.05
      .format.two.stars <<- 0.01
      .format.three.stars <<- 0.001
      .format.note.alignment <<- "l"
      .format.s.stat.parts <<- c("-!","stat names","-!","statistics1","-!")
    }
  
    # "demography" = Demography
    else if (style == "demography") {
      .format.table.parts <<- c("-!","dependent variables","models","numbers","-","coefficients","omit","additional","N","R-squared","adjusted R-squared","max R-squared","log likelihood","theta(se)*", "SER(df)","F statistic(df)*","chi2(df)*","Wald(df)*","LR(df)*","logrank(df)*","AIC","BIC","UBRE","-!","notes")    
      .format.models.skip.if.one <<- TRUE
      .format.models.left <<- ""
      .format.models.right <<- ""
      .format.numbers.left <<- "Model "
      .format.numbers.right <<- ""
      .format.coefficient.table.parts <<- c("variable name","coefficient*","standard error")
      .format.N <<- "\\textit{N}"
      .format.AIC <<- "AIC"
      .format.BIC <<- "BIC"
      .format.chi.stat <<- "Chi-Square"
      .format.note <<- ""
      .format.note.content <<- c("$^{*}$p $<$ [.*]; $^{**}$p $<$ [.**]; $^{***}$p $<$ [.***]")
      .format.one.star <<- 0.05
      .format.two.stars <<- 0.01
      .format.three.stars <<- 0.001
      .format.note.alignment <<- "l"
      .format.s.stat.parts <<- c("-!","stat names","-!","statistics1","-!")
    }
  
    # io = International Organization
    else if (style == "io") {
      .format.table.parts <<- c("-!","dependent variables","models","numbers","-","coefficients","omit","additional","N","R-squared","adjusted R-squared","max R-squared","log likelihood","theta(se)*", "SER(df)","F statistic(df)*","chi2(df)*","Wald(df)*","LR(df)*","logrank(df)*","AIC","BIC","UBRE","-!","notes")    
      .format.models.skip.if.one <<- TRUE
      .format.coefficient.table.parts <<- c("variable name","coefficient*","standard error")
      .format.coefficient.variables.capitalize <<- TRUE
      .format.s.coefficient.variables.capitalize <<- TRUE
      .format.intercept.name <<- "\\textit{Constant}"
      .format.N <<- "\\textit{Observations}"
      .format.AIC <<- "\\textit{Akaike information criterion}"
      .format.BIC <<- "\\textit{Bayesian information criterion}"
      .format.chi.stat <<- "\\textit{Chi-square}"
      .format.logrank.stat <<- "\\textit{Score (logrank) test}"
      .format.lr.stat <<- "\\textit{LR test}"
      .format.max.R2 <<- "\\textit{Maximum R-squared}"
      .format.R2 <<- "\\textit{R-squared}"
      .format.adj.R2 <<- "\\textit{Adjusted R-squared}"
      .format.UBRE <<- "\\textit{UBRE}"
      .format.F.stat <<- "\\textit{F statistic}"
      .format.LL <<- "\\textit{Log likelihood}"
      .format.SER <<- "\\textit{Residual standard error}"
      .format.null.deviance <<- "\\textit{Null deviance}"
      .format.residual.deviance <<- "\\textit{Residual deviance}"
      .format.scale <<- "\\textit{Scale}"
      .format.wald.stat <<- "\\textit{Wald test}"
      .format.note <<- "\\textit{Notes:}"
      .format.note.content <<- c("$^{***}$p $<$ [.***]; $^{**}$p $<$ [.**]; $^{*}$p $<$ [.*]")
      .format.note.alignment <<- "l"
      .format.s.stat.parts <<- c("-!","stat names","-!","statistics1","-!")
    }
  
  
    # jpam = Journal of Policy Analysis and Management
    else if (style == "jpam") {
      .format.table.parts <<- c("-!","dependent variables","models","numbers","-","coefficients","omit","additional","N","R-squared","adjusted R-squared","max R-squared","log likelihood","theta(se)*", "SER(df)","F statistic(df)*","chi2(df)*","Wald(df)*","LR(df)*","logrank(df)*","AIC","BIC","UBRE","-!","notes")    
      .format.models.skip.if.one <<- TRUE
      .format.models.left <<- ""
      .format.models.right <<- ""
      .format.numbers.left <<- "Model "
      .format.numbers.right <<- ""
      .format.numbers.roman <<- TRUE
      .format.coefficient.table.parts <<- c("variable name","coefficient*","standard error")
      .format.intercept.bottom <<- FALSE
      .format.N <<- "N"
      .format.AIC <<- "AIC"
      .format.BIC <<- "BIC"
      .format.note <<- "\\textit{Note:}"
      .format.note.content <<- c("$^{***}$p $<$ [.***]; $^{**}$p $<$ [.**]; $^{*}$p $<$ [.*]")
      .format.note.alignment <<- "l"
      .format.s.stat.parts <<- c("-!","stat names","-!","statistics1","-!")
      .format.s.statistics.names <<- cbind(c("N","N"), c("nmiss","missing"), c("mean","Mean"), c("sd","SD"), c("median","Median"), c("min","Minimum"), c("max","Maximum"), c("mad","Median Abs. Dev."), c("p","Percentile(!)"))
      
    }
  
    # "qje" = Quarterly Journal of Economics
    else if (style=="qje") {
      .format.table.parts <<- c("=!","dependent variables","models","numbers","-","coefficients","omit","additional","N","R-squared","adjusted R-squared","max R-squared","log likelihood","theta(se)*", "AIC","BIC","UBRE", "SER(df)","F statistic(df)*","chi2(df)*","Wald(df)*","LR(df)*","logrank(df)*","=!","notes")    
      .format.s.stat.parts <<- c("-!","stat names","=!","statistics1","=!")
      .format.N <<- "\\textit{N}"
      .format.note <<- "\\textit{Notes:}"
      .format.note.content <<- c("$^{***}$Significant at the [***] percent level.", "$^{**}$Significant at the [**] percent level.", "$^{*}$Significant at the [*] percent level.") 
    }
  
    # find style based on journal ("default" or other)
    else if (style=="commadefault") {
      .format.table.parts <<- c("=!","dependent variable label","dependent variables","models","numbers","-","coefficients","-","omit","-","additional","N","R-squared","adjusted R-squared","max R-squared","log likelihood","sigma2","theta(se)*", "AIC","BIC","UBRE", "SER(df)","F statistic(df)*","chi2(df)*","Wald(df)*","LR(df)*","logrank(df)*","=!","notes")
      .format.digit.separator <<- " "
      .format.decimal.character <<- ","
    }
  
    else if (style=="default") {
      .format.table.parts <<- c("=!","dependent variable label","dependent variables","models","numbers","-","coefficients","-","omit","-","additional","N","R-squared","adjusted R-squared","max R-squared","log likelihood","sigma2","theta(se)*", "AIC","BIC","UBRE", "SER(df)","F statistic(df)*","chi2(df)*","Wald(df)*","LR(df)*","logrank(df)*","=!","notes")
    }
  }

  .AIC <-
  function(object.name) {
  
    model.name <- .get.model.name(object.name)
  
    if (model.name %in% c("lmer","glmer","nlmer", "ergm", "gls")) {
      return(as.vector(AIC(object.name)))
    }
    
    if (model.name %in% c("arima")) {
      return(as.vector(object.name$aic))
    }
    else if (!is.null(.summary.object$aic)) {
      return(as.vector(.summary.object$aic)) 
    }
    else if (!is.null(object.name$AIC)) {
      return(as.vector(object.name$AIC)) 
    }
    
    return(NA)
  }
  
  .BIC <-
    function(object.name) {
      
      model.name <- .get.model.name(object.name)
      
      if (model.name %in% c("lmer","glmer","nlmer", "ergm", "gls")) {
        return(as.vector(BIC(object.name)))
      }
      
      if (model.name %in% c("arima")) {
        return(as.vector(object.name$bic))
      }
      else if (!is.null(.summary.object$bic)) {
        return(as.vector(.summary.object$bic)) 
      }
      else if (!is.null(object.name$BIC)) {
        return(as.vector(object.name$BIC)) 
      }
      
      return(NA)
    }
  

  .chi.stat <-
  function(object.name) {
    chi.output <- as.vector(rep(NA,times=3))
  
    model.name <- .get.model.name(object.name)
  
    if (!(model.name %in% c("arima","lmer", "glmer", "nlmer", "normal.gam","logit.gam","probit.gam","poisson.gam","gam()"))) {
      if (!is.null(.summary.object$chi)) {
        chi.value <- suppressMessages(.summary.object$chi)
        df.value <- suppressMessages(.summary.object$df) - suppressMessages(.summary.object$idf)
        chi.p.value <- pchisq(chi.value, df.value, ncp=0, lower.tail = FALSE, log.p = FALSE)
        chi.output <- as.vector(c(chi.value, df.value, chi.p.value))
      }
    }
  
    names(chi.output) <- c("statistic","df1","p-value")
    return(cbind(chi.output))
  }

  .coefficient.table.part <-
  function(part, which.variable, variable.name=NULL) {
    
    # decide whether to enter math mode, or not
    if (.format.dec.mark.align==TRUE) {
      c <- ""
    }
    else {
      c <- "$"
    }

  	# coefficient variable name
  	if (part=="variable name") {
    
  		# use intercept name for intercept, otherwise variable name
      if (is.na(.format.covariate.labels[.which.variable.label])) {
        if (.global.coefficient.variables[which.variable] %in% .global.intercept.strings) { cat(" ", .format.coefficient.variables.left, .format.intercept.name, .format.coefficient.variables.right, sep="") }
        else if (.format.coefficient.variables.capitalize == TRUE) { cat(" ", .format.coefficient.variables.left, toupper(variable.name), .format.coefficient.variables.right, sep="") }
  		  else { cat(" ", .format.coefficient.variables.left, variable.name, .format.coefficient.variables.right, sep="") }
      }
      else { cat(" ", .format.coefficient.variables.left, .format.covariate.labels[.which.variable.label], .format.coefficient.variables.right, sep="") }
  	}
	
  	# coefficients and stars
  	else if ((part=="coefficient") | (part=="coefficient*")) {
  		for (i in seq(1:length(.global.models))) {
  			if (!is.na(.global.coefficients[.global.coefficient.variables[which.variable],i])) {
				
  				# report the coefficient
  				cat(" & ",c, .iround(.global.coefficients[.global.coefficient.variables[which.variable],i],.format.round.digits),sep="")

  				# add stars to denote statistical significance
  				if (part=="coefficient*") { 
  					p.value <- .global.p.values[.global.coefficient.variables[which.variable],i]
  					.enter.significance.stars(p.value, math.omit=TRUE) 
  				}
          cat(c)
			
  			}
  			else {
  				cat(" & ",sep="")
  			}
  		}
  		cat(" \\\\ \n ")
  	}

  	# standard errors
  	else if ((part=="standard error") | (part=="standard error*")) {
  		for (i in seq(1:length(.global.models))) {
  			if (!is.na(.global.std.errors[.global.coefficient.variables[which.variable],i])) {

  				# report standard errors
  				cat(" & ",c, .format.std.errors.left, .iround(.global.std.errors[.global.coefficient.variables[which.variable],i],.format.round.digits),.format.std.errors.right,sep="")

  				# add stars to denote statistical significance
  				if (part=="standard error*") { 
  					p.value <- .global.p.values[.global.coefficient.variables[which.variable],i]
  					.enter.significance.stars(p.value, math.omit=TRUE) 
  				}
          cat(c)

  			}
  			else {
  				cat(" & ",sep="")
  			}
  		}
  		cat(" \\\\ \n ")
  	}


  	# p-values
  	else if ((part=="p-value") | (part=="p-value*")) {
  		for (i in seq(1:length(.global.models))) {
  			if (!is.na(.global.p.values[.global.coefficient.variables[which.variable],i])) {

  				# report p-values
  				cat(" & ",c, .format.p.values.left, .iround(.global.p.values[.global.coefficient.variables[which.variable],i],.format.round.digits,round.up.positive=TRUE),.format.p.values.right,sep="")

  				# add stars to denote statistical significance
  				if (part=="p-value*") { 
  					p.value <- .global.p.values[.global.coefficient.variables[which.variable],i]
  					.enter.significance.stars(p.value, math.omit=TRUE) 
  				}
          cat(c)

  			}
  			else {
  				cat(" & ",sep="")
  			}
  		}
  		cat(" \\\\ \n ")
  	}

  	# t-statistics
  	else if ((part=="t-stat") | (part=="t-stat*")) {
  		for (i in seq(1:length(.global.models))) {
  			if (!is.na(.global.t.stats[.global.coefficient.variables[which.variable],i])) {
  				# report t-statistics
  				cat(" & ",c, .format.t.stats.left, .iround(.global.t.stats[.global.coefficient.variables[which.variable],i],.format.round.digits),.format.t.stats.right,sep="")

  				# add stars to denote statistical significance
  				if (part=="t-stat*") { 
  					p.value <- .global.p.values[.global.coefficient.variables[which.variable],i]
  					.enter.significance.stars(p.value, math.omit=TRUE) 
  				}
          cat(c)

  			}
  			else {
  				cat(" & ",sep="")
  			}
  		}
  		cat(" \\\\ \n ")
  	}


  	# empty line
  	else if (part==" ") {
  		.table.empty.line()
  	}

  	# horizontal line
  	else if (part=="-") {
  		cat("\\hline ")
  		.table.insert.space()
  		cat(" \n")
  	}

  	# double horizontal line
  	else if (part=="=") {
  		cat("\\hline \n") 
  		cat("\\hline ")
  		.table.insert.space()
  		cat(" \n")
  	}

  }

  .coefficient.variables <-
  function(object.name) {
	
  	model.name <- .get.model.name(object.name)

  	if (model.name %in% c("ls", "normal", "logit", "probit", "relogit", "poisson", "negbin", "normal.gee", "logit.gee", "probit.gee", "poisson.gee", "normal.gam", 
  				    "logit.gam", "probit.gam", "poisson.gam", "normal.survey", "poisson.survey", "probit.survey", "logit.survey", "gamma", "gamma.gee", "gamma.survey",
  				    "exp", "weibull", "coxph", "clogit", "lognorm", "tobit", "tobit(AER)", "glm()", "svyglm()", "gee()", "survreg()", "gam()", "plm", "ivreg", "pmg", "lmrob", "gls")) {
  		return(as.vector(names(object.name$coefficients)))
  	}
    else if (model.name %in% c("clm")) {
      if (.format.ordered.intercepts == FALSE) { return(as.vector(names(object.name$beta))) }
      else { return(c(as.vector(names(object.name$beta)), as.vector(names(object.name$alpha)))) }
    }
    else if (model.name %in% c("lmer", "glmer", "nlmer")) {
      return(as.vector(names(object.name@fixef)))
    }
    else if (model.name %in% c("ergm")) {
      return(as.vector(names(object.name$coef)))
    }
    else if (model.name %in% c("betareg")) {
      return(as.vector(names(object.name$coefficients$mean)))
    }
  	else if (model.name %in% c("zeroinfl", "hurdle")) {
      if (.global.zero.component==FALSE) {
        return(as.vector(names(object.name$coefficients$count)))
      }
      else {
        return(as.vector(names(object.name$coefficients$zero)))
      }
  	}
  	else if (model.name %in% c("cloglog.net", "gamma.net", "logit.net", "probit.net")) {
  		return(as.vector(rownames(.summary.object$coefficients))) 
  	}
    else if (model.name %in% c("rlm")) {
      return(as.vector(rownames(suppressMessages(.summary.object$coefficients))))
    }
  	else if (model.name %in% c("ologit", "oprobit", "polr()")) {
  		coef.temp <- as.vector(rownames(suppressMessages(.summary.object$coefficients)))
  		if (.format.ordered.intercepts == FALSE) { return(coef.temp[seq(from=1, to=length(coef.temp)-(length(suppressMessages(.summary.object$lev))-1))]) }
  		else { return(coef.temp) }
  	}
  	else if (model.name %in% c("arima")) {
  		return(as.vector(names(object.name$coef)))
  	}
    else if (model.name %in% c("multinom")) {
      return(as.vector(object.name$coefnames))
    }
  	
  	return(NULL)
  }

  .dependent.variable <-
  function(object.name) {
    
    model.name <- .get.model.name(object.name)
    
    if (model.name %in% c("lmer", "glmer", "nlmer", "gls")) {
      return(as.vector(as.character(formula(object.name))[2]))
    }
    
    if (!is.null(object.name$call$formula)) {
      if (is.symbol(object.name$call$formula)) {
        formula.temp <- as.formula(object.name)  
      }
      else {
        formula.temp <- object.name$call$formula
      }
      
      if (length(as.vector(as.character(formula.temp)))>1) {
        return(as.vector(as.character(formula.temp)[2]))
      }
    }
    if (!is.null(object.name$formula)) {
      if (is.symbol(object.name$formula)) {
        formula.temp <- as.formula(object.name)  
      }
      else {
        formula.temp <- object.name$formula
      }
      
      if (length(as.vector(as.character(formula.temp)))>1) {   # this is for zelig$result ones
        return(as.vector(as.character(formula.temp)[2])) 
      }
    }
    if (!is.null(object.name$formula2)) {
      if (is.symbol(object.name$formula2)) {
        formula.temp <- as.formula(object.name)  
      }
      else {
        formula.temp <- object.name$formula2
      }
      
      if (length(as.vector(as.character(formula.temp)))>1) {   # z.ls
        return(as.vector(as.character(formula.temp)[2])) 
      }      
    }
      
  }
  
  .dependent.variable.written <-
  function(object.name) {
	
  	model.name <- .get.model.name(object.name)

  	if (model.name %in% c("tobit","ologit","oprobit", "relogit", "coxph","exp","lognorm","weibull","survreg()","arima")) {
  		written.var <- .inside.bracket(.dependent.variable(object.name))[1] 
  	}
  	else if (model.name %in% c("clogit")) {
  	  written.var <- .inside.bracket(.dependent.variable(object.name))[2] 
  	}
  	else { written.var <- .dependent.variable(object.name) }
    
    # some formatting changes
  	# remove everything before and including he last dollar sign from variable name
  	temp <- strsplit(written.var,"$",fixed=TRUE)
  	written.var <- temp[[1]][length(temp[[1]])]
  	
  	# if underscore or ^ in variable name, then insert an escape \ before it
  	written.var <- .remove.special.chars(written.var)
  	
    return(written.var)
  }

  .enter.significance.stars <-
  function(p.value, math.omit) {
    if ((!is.na(p.value)) & (!is.null(p.value))) {
      
      if (math.omit==TRUE) {
        c <- ""
      }
      else {
        c <- "$"
      }
    
      # three stars
      if (!is.na(.format.three.stars)) {
        if (p.value < .format.three.stars) { cat(c,"^{",.format.significance.star,.format.significance.star,.format.significance.star,"}",c,sep="") }
      }
    
      if (!is.na(.format.two.stars)) {
        if (!is.na(.format.three.stars)) {
          if ((p.value >= .format.three.stars) & (p.value < .format.two.stars)) { cat(c,"^{",.format.significance.star,.format.significance.star,"}",c,sep="") }
        }
        else {
          if (p.value < .format.two.stars) { cat(c,"^{",.format.significance.star,.format.significance.star,"}",c,sep="") }
        }
      }
    
      if (!is.na(.format.one.star)) {
        if ((!is.na(.format.two.stars)) & (!is.na(.format.three.stars))) {
          if ((p.value >= .format.two.stars) & (p.value >= .format.three.stars) & (p.value < .format.one.star)) { cat(c,"^{",.format.significance.star,"}",c,sep="") }
        }
        else if ((!is.na(.format.two.stars)) & (is.na(.format.three.stars))) {
          if ((p.value >= .format.two.stars) & (p.value < .format.one.star)) { cat(c,"^{",.format.significance.star,"}",c,sep="") }
        }
        else if ((is.na(.format.two.stars)) & (!is.na(.format.three.stars))) {
          if ((p.value >= .format.three.stars) & (p.value < .format.one.star)) { cat(c,"^{",.format.significance.star,"}",c,sep="") }
        }
        else {
          if (p.value < .format.one.star) { cat(c,"^{",.format.significance.star,"}",c,sep="") }
        }
    
      }
    
    }

  }

  .F.stat <-
  function(object.name) {
  	F.stat.output <- as.vector(rep(NA,times=4))

  	model.name <- .get.model.name(object.name)

  	if (!(model.name %in% c("arima", "lmer", "glmer", "nlmer"))) {
      if (model.name %in% c("plm")) {
        F.stat.value <- .summary.object$fstatistic$statistic
        df.numerator <- .summary.object$fstatistic$parameter["df1"]
        df.denominator <- .summary.object$fstatistic$parameter["df2"]
        F.stat.p.value <- .summary.object$fstatistic$p.value
        
        F.stat.output <- as.vector(c(F.stat.value, df.numerator, df.denominator, F.stat.p.value))
      }
  		else if (!is.null(suppressMessages(.summary.object$fstatistic["value"]))) {
  			F.stat.value <- .summary.object$fstatistic["value"]
  			df.numerator <- .summary.object$fstatistic["numdf"]
  			df.denominator <- .summary.object$fstatistic["dendf"]
  			F.stat.p.value <- pf(F.stat.value, df.numerator, df.denominator, lower.tail=FALSE)

  			F.stat.output <- as.vector(c(F.stat.value, df.numerator, df.denominator, F.stat.p.value))
  		}
  	}

  	names(F.stat.output) <- c("statistic","df1","df2","p-value")
  	return(cbind(F.stat.output))
  }

  .gcv.UBRE <-
  function(object.name) {
  
    model.name <- .get.model.name(object.name)
  
    if (!(model.name %in% c("arima", "lmer", "glmer", "nlmer"))) {
      if (!is.null(object.name$gcv.ubre)) {
        return(as.vector(object.name$gcv.ubre))
      }
    }
    return(NA)
  }

  .get.model.name <-
  function(object.name) {
  	return.value <- .model.identify(object.name)
  	if (substr(return.value,1,5)=="glm()") { return.value <- "glm()" }
  	if (substr(return.value,1,8)=="svyglm()") { return.value <- "svyglm()" }
  	if (substr(return.value,1,5)=="gee()") { return.value <- "gee()" }
  	if (substr(return.value,1,5)=="gam()") { return.value <- "gam()" }
  	if (substr(return.value,1,6)=="polr()") { return.value <- "polr()" }
  	if (substr(return.value,1,9)=="survreg()") { return.value <- "survreg()" }
  	return(return.value)
  }

  .get.p.values.1 <-
  function(object.name, user.given=NULL, auto=TRUE, coef=NULL, se=NULL) {

    if (!is.null(user.given)) { return(user.given) }
        
    if (auto == TRUE) {
      if (!is.null(se)) {
        t <- (coef / se)
        return( 2*pnorm(abs(t), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE) )
      }
    }

    model.name <- .get.model.name(object.name)
    
  	if (model.name %in% c("ls", "normal", "logit", "probit", "relogit", "poisson", "negbin", "normal.survey", "poisson.survey", "probit.survey", "logit.survey", "gamma", "gamma.survey",
                            "cloglog.net", "gamma.net", "logit.net", "probit.net", "glm()", "svyglm()", "plm", "ivreg", "lmrob")) {
  		return(.summary.object$coefficients[,4])
  	}
  	if (model.name %in% c("lmer", "glmer", "nlmer")) {
  	  Vcov <- as.matrix(vcov(object.name, useScale = FALSE))
  	  coefs <- object.name@fixef
  	  se <- sqrt(diag(Vcov))
  	  tstat <- coefs / se
  	  pval <- 2 * pnorm(abs(tstat), lower.tail = FALSE)
      names(pval) <- names(coefs)
  	  return(pval)
  	}
  	if (model.name %in% c("ergm")) {
  	  return(.summary.object$coefs[,4])
  	}
    if (model.name %in% c("clm")) {
      if (.format.ordered.intercepts == FALSE) {
        return(.summary.object$coefficients[(length(object.name$alpha)+1):(length(object.name$coefficients)),4])
      }
      else {
        return(.summary.object$coefficients[,4])
      }
    }
    else if (model.name %in% c("pmg")) {
      coef.temp <- .summary.object$coefficients
      std.err.temp <- sqrt(diag(.summary.object$vcov))
      t.stat.temp <- coef.temp / std.err.temp
      df.temp <- length(.summary.object$residuals)
      return( 2 * pt(abs(t.stat.temp), df=df.temp, lower.tail = FALSE, log.p = FALSE) )
    }
    else if (model.name %in% c("zeroinfl", "hurdle")) {
      if (.global.zero.component==FALSE) {
        return(.summary.object$coefficients$count[,4])  
      }
      else {
        return(.summary.object$coefficients$zero[,4])
      }
      
    }
  	else if (model.name %in% c("normal.gee", "logit.gee", "poisson.gee", "probit.gee", "gamma.gee", "gee()")) {
  		return(2*pnorm(abs(.summary.object$coefficients[,"Robust z"]), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
  	}
  	else if (model.name %in% c("normal.gam", "logit.gam", "probit.gam", "poisson.gam", "gam()")) {
  		return(.summary.object$p.pv)
  	}
  	else if (model.name %in% c("coxph", "clogit")) {
  		return(.summary.object$coef[,"Pr(>|z|)"])
  	}
  	else if (model.name %in% c("exp","lognorm","weibull","tobit", "survreg()")) {
  		return(.summary.object$table[,"p"])
  	}
    else if (model.name %in% c("rlm")) {
      coef.temp <- suppressMessages(.summary.object$coefficients[,"t value"])
      coef.temp <- 2*pnorm(abs(coef.temp[seq(from=1, to=length(coef.temp))]), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
      return(coef.temp)
    }
  	else if (model.name %in% c("ologit", "oprobit", "polr()")) {
  		coef.temp <- suppressMessages(.summary.object$coefficients[,"t value"])
  		if (.format.ordered.intercepts == FALSE) { return(2*pnorm(abs(coef.temp[seq(from=1, to=length(coef.temp)-(length(suppressMessages(.summary.object$lev))-1))]), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)) }
  		else { 
  		  return( 2*pnorm(abs(coef.temp[seq(from=1, to=length(coef.temp))]), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE) ) 
      }
  		
  	}
  	else if (model.name %in% c("arima")) {
  		return(2*pnorm( abs(object.name$coef / (sqrt(diag(object.name$var.coef))) ), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
  	}
    else if (model.name %in% c("tobit(AER)")){
      return(.summary.object$coefficients[,"Pr(>|z|)"])
    }
    else if (model.name %in% c("multinom")) {
      if (is.null(nrow(.summary.object$coefficients))) {
        coef.temp <- .summary.object$coefficients
        se.temp <- .summary.object$standard.errors
      }
      else {
        coef.temp <- .summary.object$coefficients[1,]
        se.temp <- .summary.object$standard.errors[1,]
      }
      return( 2*pnorm( abs( (coef.temp) / (se.temp) ) , mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE) )
    }
  	else if (model.name %in% c("betareg")) {
  	  return(.summary.object$coefficients$mean[,"Pr(>|z|)"])
  	}
    else if (model.name %in% c("gls")) {
      coef.temp <- object.name$coefficients
      se.temp <- sqrt(diag(object.name$varBeta))
      t.temp <- coef.temp / se.temp
      p.temp <- 2*pnorm( abs( t.temp ) , mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
      return(p.temp)
    }
    return(NULL)
  }
  
  .get.p.values <-
  function(object.name, user.given=NULL, auto=TRUE, coef=NULL, se=NULL) {
      out <- .get.p.values.1(object.name, user.given, auto, coef, se)
      
      coef.vars <- .coefficient.variables(object.name)
      if (is.null(names(out))) {
        if (length(out) < length(coef.vars)) {
          out.temp <- rep(NA, times=length(coef.vars)-length(out))
          out <- c(out, out.temp)
        }
        else if (length(out) > length(coef.vars)) {
          out <- out[1:length(coef.vars)]
        }
        names(out) <- coef.vars
      }
      return(out)
  }
  

  .get.scale <-
  function(object.name) {
  
    model.name <- .get.model.name(object.name)
  
    if (!(model.name %in% c("arima", "lmer", "glmer", "nlmer"))) {
      if (!is.null(object.name$scale)) {
        if (model.name %in% c("normal.gee", "logit.gee", "poisson.gee", "probit.gee", "gamma.gee", "gee()", "exp","lognorm","weibull","tobit","survreg()","tobit(AER)")) {
          return(as.vector(object.name$scale))
        }
      }
    }
    return(NA)
  }

  .get.sigma2 <-
  function(object.name) {
  
    model.name <- .get.model.name(object.name)
    
    if (model.name %in% c("lmer", "glmer", "nlmer")) {
      return(NA)
    }
  
    if (!is.null(object.name$sigma2)) {
        return(as.vector(object.name$sigma2))
    }
    return(NA)
  }

  .get.standard.errors.1 <-
  function(object.name, user.given=NULL) {
    
    if (!is.null(user.given)) { return(user.given) }

  	model.name <- .get.model.name(object.name)

  	if (model.name %in% c("ls", "normal", "logit", "probit", "relogit", "poisson", "negbin", "normal.survey", "poisson.survey", "probit.survey", "logit.survey", "gamma", "gamma.survey",
                            "cloglog.net", "gamma.net", "logit.net", "probit.net", "glm()", "svyglm()", "plm", "ivreg", "lmrob")) {
  		return(.summary.object$coefficients[,"Std. Error"])
  	}
  	if (model.name %in% c("lmer", "glmer", "nlmer")) {
  	  Vcov <- as.matrix(vcov(object.name, useScale = FALSE))
  	  coefs <- object.name@fixef
  	  se <- sqrt(diag(Vcov))
      names(se) <- names(coefs)
  	  return(se)
  	}
  	if (model.name %in% c("ergm")) {
  	  return(.summary.object$coefs[,2])
  	}
  	if (model.name %in% c("clm")) {
  	  if (.format.ordered.intercepts == FALSE) {
  	    return(.summary.object$coefficients[(length(object.name$alpha)+1):(length(object.name$coefficients)),2])
  	  }
  	  else {
  	    return(.summary.object$coefficients[,2])
  	  }
  	}
  	else if (model.name %in% c("pmg")) {
  	  return (sqrt(diag(.summary.object$vcov)))
  	}
  	if (model.name %in% c("zeroinfl", "hurdle")) {
      if (.global.zero.component == FALSE) {
        return(.summary.object$coefficients$count[,"Std. Error"])  
      }
      else {
        return(.summary.object$coefficients$zero[,"Std. Error"])
      }
  	}
  	else if (model.name %in% c("normal.gee", "logit.gee", "poisson.gee",  "probit.gee", "gamma.gee", "gee()")) {
  		return(.summary.object$coefficients[,"Robust S.E."])
  	}
  	else if (model.name %in% c("normal.gam", "logit.gam", "probit.gam", "poisson.gam", "gam()")) {
  	  temp.se <- .summary.object$se
      names(temp.se) <- names(.summary.object$p.coeff)
      return(temp.se)
  	}
  	else if (model.name %in% c("coxph")) {
  		return(.summary.object$coef[,"robust se"])
  	}
    else if (model.name %in% c("clogit")) {
      return(.summary.object$coef[,"se(coef)"])
      
    }
  	else if (model.name %in% c("exp","lognorm","weibull","tobit","survreg()")) {
  		return(.summary.object$table[,"Std. Error"])
  	}
    else if (model.name %in% c("rlm")) {
      return(suppressMessages(.summary.object$coefficients[,"Std. Error"]))
    }
  	else if (model.name %in% c("ologit", "oprobit", "polr()")) {
  		se.temp <- suppressMessages(.summary.object$coefficients[,"Std. Error"])
  		if (.format.ordered.intercepts == FALSE) { return(se.temp[seq(from=1, to=length(se.temp)-(length(suppressMessages(.summary.object$lev))-1))]) }
  		else { return(se.temp) }
  	}
  	else if (model.name %in% c("arima")) {
  		return( sqrt(diag(object.name$var.coef)) )
  	}
  	else if (model.name %in% c("tobit(AER)")){
  	  return(.summary.object$coefficients[,"Std. Error"])
  	}
  	else if (model.name %in% c("multinom")) {
  	  if (is.null(nrow(.summary.object$coefficients))) {
  	    se.temp <- .summary.object$standard.errors
  	  }
  	  else {
  	    se.temp <- .summary.object$standard.errors[1,]
  	  }
  	  return(se.temp)
  	}
  	else if (model.name %in% c("betareg")) {
  	  return(.summary.object$coefficients$mean[,"Std. Error"])
  	}
    else if (model.name %in% c("gls")) {
      se.temp <- sqrt(diag(object.name$varBeta))
      return(se.temp)
    }
    return(NULL)
  }
  
  .get.standard.errors <-
  function(object.name, user.given=NULL) {
      out <- .get.standard.errors.1(object.name, user.given)
      
      coef.vars <- .coefficient.variables(object.name)
      if (is.null(names(out))) {
        if (length(out) < length(coef.vars)) {
          out.temp <- rep(NA, times=length(coef.vars)-length(out))
          out <- c(out, out.temp)
        }
        else if (length(out) > length(coef.vars)) {
          out <- out[1:length(coef.vars)]
        }
        names(out) <- coef.vars
      }
      return(out)
  }

  .get.t.stats.1 <-
  function(object.name, user.given=NULL, auto=TRUE, coef=NULL, se=NULL) {
    
    if (!is.null(user.given)) { return(user.given) }
    
    if (auto == TRUE) {
      if (!is.null(se)) {
        return(coef / se)
      }
    }

  	model.name <- .get.model.name(object.name)

  	if (model.name %in% c("ls", "normal", "logit", "probit", "relogit", "poisson", "negbin", "normal.survey", "poisson.survey", "probit.survey", "logit.survey", "gamma", "gamma.survey",
      				    "cloglog.net", "gamma.net", "logit.net", "probit.net", "glm()", "svyglm()","plm", "ivreg", "lmrob")) {
  		return(.summary.object$coefficients[,3])
  	}
  	if (model.name %in% c("lmer", "glmer", "nlmer")) {
  	  Vcov <- as.matrix(vcov(object.name, useScale = FALSE))
  	  coefs <- object.name@fixef
  	  se <- sqrt(diag(Vcov))
  	  tstat <- coefs / se
  	  names(tstat) <- names(coefs)
  	  
  	  return(tstat)
  	}
  	if (model.name %in% c("ergm")) {
  	  return((.summary.object$coefs[,1])/(.summary.object$coefs[,2]))
  	}
  	if (model.name %in% c("clm")) {
  	  if (.format.ordered.intercepts == FALSE) {
  	    return(.summary.object$coefficients[(length(object.name$alpha)+1):(length(object.name$coefficients)),3])
  	  }
  	  else {
  	    return(.summary.object$coefficients[,3])
  	  }
  	}
  	else if (model.name %in% c("pmg")) {
  	  coef.temp <- .summary.object$coef
  	  std.err.temp <- sqrt(diag(.summary.object$vcov))
  	  t.stat.temp <- coef.temp / std.err.temp
  	  return(t.stat.temp)
    }
    else if (model.name %in% c("zeroinfl", "hurdle")) {
      if (.global.zero.component == FALSE) {
        return(.summary.object$coefficients$count[,3])  
      }
      else {
        return(.summary.object$coefficients$zero[,3])
      }
      
    }
  	else if (model.name %in% c("normal.gee", "logit.gee", "poisson.gee",  "probit.gee", "gamma.gee", "gee()")) {
  		return(.summary.object$coefficients[,"Robust z"])
  	}
  	else if (model.name %in% c("normal.gam", "logit.gam", "probit.gam", "poisson.gam", "gam()")) {
  		return(.summary.object$p.t)
  	}
  	else if (model.name %in% c("coxph", "clogit")) {
  		return(.summary.object$coef[,"z"])
  	}
  	else if (model.name %in% c("exp","lognorm","weibull", "tobit","survreg()")) {
  		return(.summary.object$table[,"z"])
  	}
    else if (model.name %in% c("rlm")) {
      return(suppressMessages(.summary.object$coefficients[,"t value"]))
    }
  	else if (model.name %in% c("ologit", "oprobit", "polr()")) {
  		tstat.temp <- suppressMessages(.summary.object$coefficients[,"t value"])
  		if (.format.ordered.intercepts == FALSE) { return(tstat.temp[seq(from=1, to=length(tstat.temp)-(length(suppressMessages(.summary.object$lev))-1))]) }
  		else { return(tstat.temp) }
  	}
  	else if (model.name %in% c("arima")) {
  		return( object.name$coef / (sqrt(diag(object.name$var.coef))) )
  	}
  	else if (model.name %in% c("tobit(AER)")){
  	  return(.summary.object$coefficients[,"z value"])
  	}
  	else if (model.name %in% c("multinom")) {
  	  if (is.null(nrow(.summary.object$coefficients))) {
  	    coef.temp <- .summary.object$coefficients
  	    se.temp <- .summary.object$standard.errors
  	  }
  	  else {
  	    coef.temp <- .summary.object$coefficients[1,]
  	    se.temp <- .summary.object$standard.errors[1,]
  	  }
  	  return( (coef.temp) / (se.temp) )
  	}
  	else if (model.name %in% c("betareg")) {
  	  return(.summary.object$coefficients$mean[,"z value"])
  	}
    else if (model.name %in% c("gls")) {
      coef.temp <- object.name$coefficients
      se.temp <- sqrt(diag(object.name$varBeta))
      return(coef.temp / se.temp)
    }
  	
  	return(NULL)
  }
  
  .get.t.stats <-
  function(object.name, user.given=NULL, auto=TRUE, coef=NULL, se=NULL) {
    out <- .get.t.stats.1(object.name, user.given, auto, coef, se)

    coef.vars <- .coefficient.variables(object.name)
    if (is.null(names(out))) {
      if (length(out) < length(coef.vars)) {
        out.temp <- rep(NA, times=length(coef.vars)-length(out))
        out <- c(out, out.temp)
      }
      else if (length(out) > length(coef.vars)) {
        out <- out[1:length(coef.vars)]
      }
      names(out) <- coef.vars
    }
    return(out)
  }
  

  .get.theta <-
  function(object.name) {
    theta.output <- as.vector(rep(NA,times=4))
  
    model.name <- .get.model.name(object.name)
  
    if (!(model.name %in% c("arima", "lmer", "glmer", "nlmer"))) {
      if ((!is.null(object.name$theta)) & (!is.null(object.name$SE.theta))) {
        theta.value <- object.name$theta
        theta.se.value <- object.name$SE.theta
        theta.tstat.value <- theta.value / theta.se.value
        theta.p.value <- 2*pnorm(abs(theta.tstat.value), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
      
        theta.output <- as.vector(c(theta.value, theta.se.value, theta.tstat.value, theta.p.value))
      }
    }
  
    names(theta.output) <- c("statistic","se","tstat","p-value")
    return(cbind(theta.output))
  }

  .inside.bracket <-
  function(s) {
  	process.string <- ""
  	return.vector <- NULL
  
    if (is.null(s)) { return("") }
    if (is.na(s) | (s=="")) { return("") }

  	inside.inner.bracket <- 0
  	for (i in seq(from = (regexpr("(",s,fixed=TRUE)[1])+1, to = nchar(s))) {
  		letter <- substr(s,i,i)
  		if (letter == "(") { inside.inner.bracket <- inside.inner.bracket + 1 }
  		if (letter == ")") { inside.inner.bracket <- inside.inner.bracket - 1 }
		
  		if ((letter == ",") & (inside.inner.bracket == 0)) {
  			return.vector <- c(return.vector, process.string)
  			process.string <- ""
  		}
  		else if (inside.inner.bracket >= 0) { process.string <- paste(process.string, letter, sep="") }
  		else { break } 
  	}
  	if (process.string != "") { return.vector <- c(return.vector, process.string) }
  	return (.trim(return.vector))
  }

  .iround <-
  function(x, decimal.places=0, round.up.positive=FALSE, simply.output=FALSE) {
    if (is.na(x) | is.null(x)) { return("") }
  
    if (simply.output == TRUE) {
      if (!is.numeric(x)) { return(.remove.special.chars(x)) }
    }
    
    if (!is.na(decimal.places)) {
    
    	if ((.format.until.nonzero.digit == FALSE) | (decimal.places <= 0)) {
    		round.result <- round(x, digits=decimal.places)
    	}
    	else {
    		temp.places <- decimal.places
        if (!.is.all.integers(x)) {
    		  while ((round(x, digits=temp.places) == 0) & (temp.places < (decimal.places + .format.max.extra.digits))) {
      			temp.places <- temp.places + 1
    	  	}
    	  }
    		round.result <- round(x, digits=temp.places)
    		decimal.places <- temp.places
    	}
	
    	if ((round.up.positive==TRUE) & (round.result < x)) {       # useful for p-values that should be rounded up
    		if (x > (10^((-1)*(decimal.places+1)))) {
    			round.result <- round.result + 10^((-1)*decimal.places)
    		}
    		else { round.result <- 0 }
    	}
    }
    else {      # if the decimal place is NA
      round.result <- x
    }
    
    round.result.char <- as.character(format(round.result, scientific=FALSE))
  	split.round.result <- unlist(strsplit(round.result.char, "\\."))

  	## first deal with digit separator
  
    for (i in seq(from=1, to=length(.format.digit.separator.where))) {
      if (.format.digit.separator.where[i]<=0) {
        .format.digit.separator.where[i] <<- -1
      }
    }
  
  	first.part <- ""
  	separator.count <- 1
  	length.integer.part <- nchar(split.round.result[1])

  	digits.in.separated.unit <- 0
  	for (i in seq(from=length.integer.part, to=1)) {
  		if ((digits.in.separated.unit == .format.digit.separator.where[separator.count]) & (substr(split.round.result[1],i,i)!="-")){
  			first.part <- paste(.format.digit.separator,first.part,sep="")
  			if (separator.count < length(.format.digit.separator.where)) { separator.count <- separator.count + 1 }
  			digits.in.separated.unit <- 0	
  		}
  		first.part <- paste(substr(split.round.result[1],i,i),first.part,sep="")
  		digits.in.separated.unit <- digits.in.separated.unit + 1

  	}

  	# now deal with the decimal part
    if (!is.na(decimal.places)) {
  	  if (decimal.places <= 0) {
  	  	return(first.part) 
  	  }
    }

  	# remove initial zero, if that is requested
  	if (.format.initial.zero==FALSE) {
  		if ((round.result >= 0) & (round.result < 1)) {
  			first.part <- ""
  		}
  	}	

  	if (length(split.round.result)==2) {
  	  if (is.na(decimal.places)) { return(paste(first.part,.format.decimal.character,split.round.result[2],sep="")) }
  		if (nchar(split.round.result[2]) < decimal.places) {
  			decimal.part <- split.round.result[2]
  			for (i in seq(from = 1,to = (decimal.places - nchar(split.round.result[2])))) {
  				decimal.part <- paste(decimal.part,"0", sep="")
  			}
  			return(paste(first.part,.format.decimal.character,decimal.part,sep=""))
  		}
  		else { return(paste(first.part,.format.decimal.character,split.round.result[2],sep="")) }
  	}
  	else if (length(split.round.result)==1) { 
  	  if (is.na(decimal.places)) { return(paste(first.part,.format.decimal.character,decimal.part,sep="")) }
  		decimal.part <- ""
  		for (i in seq(from = 1,to = decimal.places)) {
  			decimal.part <- paste(decimal.part,"0", sep="")
  		}
  		return(paste(first.part,.format.decimal.character,decimal.part,sep=""))
  	}
  	else { return(NULL) }
  }
   
  is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  .is.all.integers <-
  function(x) {
      if (!is.numeric(x)) { return(FALSE) }
      if (length(x) == length(is.wholenumber(x)[(!is.na(x)) & (is.wholenumber(x)==TRUE)])) {
        return(TRUE)
      }
      else { return (FALSE) }
    }
  

  .log.likelihood <-
  function(object.name) {

  	model.name <- .get.model.name(object.name)

  	if (model.name %in% c("arima", "betareg", "zeroinfl", "hurdle")) {
  		return(as.vector(object.name$loglik))
  	}
  	if (model.name %in% c("lmer", "glmer", "nlmer")) { 
  	  return(as.vector(logLik(object.name)[1]))
  	}
  	if (model.name %in% c("clm", "gls")) {
  	  return(as.vector(object.name$logLik))
  	}
  	else if (model.name %in% c("coxph", "clogit", "exp", "weibull", "lognorm","tobit", "tobit(AER)", "survreg()")) {
  		return(as.vector(.summary.object$loglik[2]))
  	}
  	else if (!is.null(object.name$aic)) {
  	  return(as.vector(-(0.5)*(object.name$aic-2*length(.summary.object$coefficients[,"Estimate"]))))
  	}
  	return(NA)
  }

  .logrank.stat <-
  function(object.name) {
    logrank.output <- as.vector(rep(NA,times=3))
  
    model.name <- .get.model.name(object.name)
  
    if (!(model.name %in% c("arima", "lmer", "glmer", "nlmer"))) {
      if (!is.null(.summary.object$logtest)) {
        logrank.value <- suppressMessages(.summary.object$sctest[1])
        df.value <- suppressMessages(.summary.object$sctest[2])
        logrank.p.value <- suppressMessages(.summary.object$sctest[3])
        logrank.output <- as.vector(c(logrank.value, df.value, logrank.p.value))
      }
    
    }
  
    names(logrank.output) <- c("statistic","df1","p-value")
    return(cbind(logrank.output))
  }

  .lr.stat <-
  function(object.name) {
    log.output <- as.vector(rep(NA,times=3))
  
    model.name <- .get.model.name(object.name)
    
    if (!(model.name %in% c("arima","lmer","glmer","nlmer"))) {
      if (!is.null(.summary.object$logtest)) {
        log.value <- suppressMessages(.summary.object$logtest[1])
        df.value <- suppressMessages(.summary.object$logtest[2])
        log.p.value <- suppressMessages(.summary.object$logtest[3])
        log.output <- as.vector(c(log.value, df.value, log.p.value))
      }
      
    }
  
    names(log.output) <- c("statistic","df1","p-value")
    return(cbind(log.output))
  }

  .max.r.squared <-
  function(object.name) {
  
    model.name <- .get.model.name(object.name)
  
    if (!(model.name %in% c("arima", "lmer", "glmer", "nlmer"))) {
      if (model.name %in% c("coxph", "clogit")) {
        return(as.vector(.summary.object$rsq[2]))
      }
    }
    return(NA)
  }
  
  .model.identify <-
  function(object.name) {
    
    if (class(object.name)[1]=="ergm") {
      return("ergm")
    }
    
   if (.hasSlot(object.name,"Zt")) {
     if (object.name@call[1]=="lmer()") {
       return("lmer")
     }
     if (object.name@call[1]=="glmer()") {
       return("glmer")
     }
     if (object.name@call[1]=="nlmer()") {
       return("nlmer")
     }     
     return("unknown")
   }
       
   if (!is.null(object.name$call)) {
    
  	if (object.name$call[1]=="lm()") { return("ls") }
  	else if (object.name$call[1]=="glm()") {
  		if (object.name$family$family=="gaussian") {
  			if (object.name$family$link=="identity") {
  				return("normal")
  			}
  		}
  		else if (object.name$family$family=="binomial") {
  			if (object.name$family$link=="probit") {
  				return("probit")
  			}
  			if (object.name$family$link=="logit") {
  				return("logit")
	  		}

	  	}
	  	else if (object.name$family$family=="poisson") {
  			if (object.name$family$link=="log") {
  				return("poisson")
  			}
  		}
  		else if (object.name$family$family=="Gamma") {
  			if (object.name$family$link=="inverse") {
  				return("gamma")
  			}
  		}
  		return(paste("glm()#",object.name$family$family,"#",object.name$family$link, sep=""))
  	}

  	else if (object.name$call[1]=="svyglm()") {
  		if (object.name$family$family=="gaussian") {
  			if (object.name$family$link=="identity") {
  				return("normal.survey")
  			}
  		}
  		else if ((object.name$family$family=="binomial") | (object.name$family$family=="quasibinomial")) {
  			if (object.name$family$link=="probit") {
  				return("probit.survey")
  			}
  			if (object.name$family$link=="logit") {
  				return("logit.survey")
  			}

  		}
  		else if (object.name$family$family=="poisson") {
  			if (object.name$family$link=="log") {
  				return("poisson.survey")
  			}
  		}
  		else if (object.name$family$family=="Gamma") {
  			if (object.name$family$link=="inverse") {
  				return("gamma.survey")
  			}
  		}
  		return(paste("svyglm()#",object.name$family$family,"#",object.name$family$link, sep=""))
  	}

  	else if (object.name$call[1]=="gam()") {
  		if (object.name$family$family=="gaussian") {
  			if (object.name$family$link=="identity") {
  				return("normal.gam")
  			}
  		}
  		else if (object.name$family$family=="binomial")  {
  			if (object.name$family$link=="probit") {
  				return("probit.gam")
  			}
  			if (object.name$family$link=="logit") {
  				return("logit.gam")
  			}

  		}
  		else if (object.name$family$family=="poisson") {
  			if (object.name$family$link=="log") {
  				return("poisson.gam")
  			}
  		}
  		else if (object.name$family$family=="Gamma") {
  			if (object.name$family$link=="inverse") {
  				return("gamma.gam")
  			}
  		}
  		return(paste("gam()#",object.name$family$family,"#",object.name$family$link, sep=""))
  	}
	
  	else if (object.name$call[1]=="polr()") {
  		if (object.name$method=="logistic") {
  			return("ologit")
  		}
  		else if (object.name$method=="probit") {
  			return("oprobit")
  		}
  		return(paste("polr()#",object.name$method, sep=""))
  	}


  	else if (object.name$call[1]=="gee()") {
  		if (object.name$family$family=="gaussian") {
  			if (object.name$family$link=="identity") {
  				return("normal.gee")
  			}
  		}
  		else if (object.name$family$family=="binomial") {
  			if (object.name$family$link=="probit") {
  				return("probit.gee")
  			}
  			if (object.name$family$link=="logit") {
  				return("logit.gee")
  			}

  		}
  		else if (object.name$family$family=="poisson") {
  			if (object.name$family$link=="log") {
  				return("poisson.gee")
  			}
  		}
  		else if (object.name$family$family=="Gamma") {
  			if (object.name$family$link=="inverse") {
  				return("gamma.gee")
  			}
  		}
  		return(paste("gee()#",object.name$family$family,"#",object.name$family$link, sep=""))
  	}

  	else if (object.name$call[1]=="survreg()") {
  		if (object.name$dist=="exponential") {
  			return("exp")
  		}
  		else if (object.name$dist=="weibull") {
  			return("weibull")
  		}
  		else if (object.name$dist=="lognorm") {
  			return("lognormal")
  		}
      else if (object.name$dist=="gaussian") {
        return("tobit")
      }
  		return(paste("survreg()#",object.name$dist, sep=""))
  	}

  	else if (object.name$call[1]=="glm.nb()") {
  		return("negbin")
  	}
  	
    if (!is.null(object.name$userCall)) {
  	  if (object.name$userCall[1]=="clogit()") {
  	    return("clogit")
  	  }
  	}
  	
    if (object.name$call[1]=="coxph()") {
  		return("coxph")
  	}
  	if (object.name$call[1]=="pmg()") {
  	  return("pmg")
  	}
  	if (object.name$call[1]=="gls()") {
  	  return("gls")
  	}
  	if (object.name$call[1]=="clm()") {
  	  return("clm")
  	}
  	if (object.name$call[1]=="lmrob()") {
  	  return("lmrob")
  	}
  	if (object.name$call[1]=="rlm()") {
  	  return("rlm")
  	}
    else if (object.name$call[1]=="relogit()") {
      return("relogit")
    }
  	else if (object.name$call[1]=="netbinom()") {
  	  if (object.name$call$LF=="probit") { return("probit.net") }      
      if (object.name$call$LF=="logit") { return("logit.net") }
  	  if (object.name$call$LF=="cloglog") { return("cloglog.net") }
  	}
  	else if (object.name$call[1]=="netgamma()") {
  	  return("gamma.net")
  	}

  	else if (object.name$call[1]=="zelig()") {
        if (object.name$call$model %in% c("ls","normal","logit","probit","relogit","poisson","poisson.survey",
                                           "negbin","probit.survey","logit.survey","normal.gee","logit.gee","probit.gee",
                                           "poisson.gee","normal.gam","logit.gam","probit.gam","poisson.gam","exp",
                                           "coxph","weibull","lognorm","normal.survey","gamma","gamma.survey",
                                           "gamma.gee","cloglog.net","logit.net","probit.net","gamma.net","ologit",
                                           "oprobit","arima","tobit")) {
            return(object.name$call$model)
  		    }
          else { return("unsupported zelig") }
  	}
    
  	else if (object.name$call[1]=="tobit()") {
  	  return("tobit(AER)")
  	}
    
    else if (object.name$call[1]=="multinom()") {
      return("multinom")
    }
    
  	else if (object.name$call[1]=="betareg()") {
  	  return("betareg")
  	}
  	else if (object.name$call[1]=="zeroinfl()") {
  	  return("zeroinfl")
  	}
  	else if (object.name$call[1]=="hurdle()") {
  	  return("hurdle")
  	}  	
  	else if (object.name$call[1]=="plm()") {
  	  return("plm")
  	}  	
  	else if (object.name$call[1]=="ivreg()") {
  	  return("ivreg")
  	} 
   }
  	
   return("unknown")
    
  }

  .new.table <-
  function(object.name, user.coef=NULL, user.se=NULL, user.t=NULL, user.p=NULL, auto.t=TRUE, auto.p=TRUE) {
    
    .summary.object <<- summary(object.name)

  	.global.models <<- suppressMessages(as.vector(.model.identify(object.name)))
  	.global.dependent.variables <<- suppressMessages(.dependent.variable(object.name))
  	.global.dependent.variables.written <<- suppressMessages(.dependent.variable.written(object.name))
  	.global.coefficient.variables <<- suppressMessages(.coefficient.variables(object.name))
  	.global.coefficients <<- suppressMessages(cbind(.get.coefficients(object.name, user.coef)))
  	.global.std.errors <<- suppressMessages(cbind(.get.standard.errors(object.name, user.se)))
    
    feed.coef <- .global.coefficients
    if (!is.null(user.coef)) { feed.coef <- user.coef }   # feed user-defined coefficients, if available
    
  	.global.t.stats <<- suppressMessages(cbind(.get.t.stats(object.name, user.t, auto.t, feed.coef, user.se)))
  	.global.p.values <<- suppressMessages(cbind(.get.p.values(object.name, user.p, auto.p, feed.coef, user.se)))

  	.global.N <<- suppressMessages(.number.observations(object.name))
  	.global.LL <<- suppressMessages(.log.likelihood(object.name))
  	.global.R2 <<- suppressMessages(.r.squared(object.name))
  	.global.max.R2 <<- suppressMessages(.max.r.squared(object.name))
  	.global.adj.R2 <<- suppressMessages(.adj.r.squared(object.name))
  	.global.AIC <<- suppressMessages(.AIC(object.name))
    .global.BIC <<- suppressMessages(.BIC(object.name))
    .global.scale <<- suppressMessages(.get.scale(object.name))
    .global.UBRE <<- suppressMessages(.gcv.UBRE(object.name))
  	.global.sigma2 <<- suppressMessages(.get.sigma2(object.name))
  	.global.theta <<- suppressMessages(.get.theta(object.name))

  	.global.SER <<- suppressMessages(.SER(object.name))
  	.global.F.stat <<- suppressMessages(.F.stat(object.name))
    .global.chi.stat <<- suppressMessages(.chi.stat(object.name))
  	.global.wald.stat <<- suppressMessages(.wald.stat(object.name))
  	.global.lr.stat <<- suppressMessages(.lr.stat(object.name))
  	.global.logrank.stat <<- suppressMessages(.logrank.stat(object.name))
  	.global.null.deviance <<- suppressMessages(.null.deviance(object.name))
  	.global.residual.deviance <<- suppressMessages(.residual.deviance(object.name))

  }

  .null.deviance <-
  function(object.name) {
  	null.deviance.output <- as.vector(rep(NA,times=3))

  	model.name <- .get.model.name(object.name)

  	if (!(model.name %in% c("arima","lmer","glmer","nlmer", "ergm"))) {
  		if (!is.null(suppressMessages(.summary.object$null.deviance))) {
  			null.deviance.value <- suppressMessages(.summary.object$null.deviance)
  			df.value <- object.name$df.null

  			null.deviance.output <- as.vector(c(null.deviance.value, df.value, NA))
  		}
  		else if (!is.null(object.name$null.deviance)) {
  		  null.deviance.value <- object.name$null.deviance
  		  df.value <- object.name$df.null
		  
  		  null.deviance.output <- as.vector(c(null.deviance.value, df.value, NA))
  		}
  	}

  	names(null.deviance.output) <- c("statistic","df1","p-value")
  	return(cbind(null.deviance.output))
  }

  .number.observations <-
  function(object.name) {
  
    model.name <- .get.model.name(object.name)
  
    if (model.name %in% c("ls", "normal", "logit", "probit", "relogit",
                          "poisson", "negbin", "normal.survey", "poisson.survey",
                          "probit.survey", "logit.survey", "gamma", "gamma.survey",
                          "z.arima", "glm()", "svyglm()")) {
      return(as.vector(object.name$df.residual+length(object.name$coefficients)))
    }
    else if (model.name %in% c("lmer","glmer","nlmer")) {
      return(length(resid(object.name)))  
    }
    else if (model.name %in% c("plm", "pmg", "rlm", "lmrob")) {
      return(as.vector(length(object.name$residual)))
    }
    else if (model.name %in% c("hurdle", "zeroinfl")) {
      return(as.vector(object.name$n))
    }
    else if (model.name %in% c("ivreg","clm")) {
      return(as.vector(object.name$nobs))
    }
    if (model.name %in% c("normal.gee", "logit.gee", "poisson.gee",
                          "probit.gee", "gamma.gee", "gee()", "betareg")) {
      return(as.vector(.summary.object$nobs))
    }
    else if (model.name %in% c("normal.gam", "logit.gam", "probit.gam",
                               "poisson.gam", "coxph", "clogit", "exp", "lognorm", "weibull", "survreg()",
                               "gam()")) {
      return(as.vector(.summary.object$n))
    }
    else if (model.name %in% c("ologit", "oprobit", "polr()")) {
      return(as.vector(.summary.object$nobs))
    }
    else if (model.name %in% c("tobit(AER)")) {
      return(as.vector(.summary.object$n["Total"]))
    }
    return(NA)
  }

  .print.additional.lines <-
  function(part.number=NULL) {

  	# if no additional lines, then quit the function
  	if ((is.null(.format.additional.lines.titles)==TRUE) & (is.null(.format.additional.lines.values)==TRUE)) {
  		return(NULL)
  	}

  	## make matrix the right size - add columns if too few, remove if too many
  	if (length(.global.models) > ncol(.format.additional.lines.values)) {
  		for (i in seq(1:(ncol(.format.additional.lines.values)-length(.global.models)))) {
  			.format.additional.lines.values <<- cbind(.format.additional.lines.values,as.vector(rep(NA, times=nrow(.format.additional.lines.values))))
  		}		
  	}
  	else if (length(.global.models) < ncol(.format.additional.lines.values)) {
  		.format.additional.lines.values <<- .format.additional.lines.values[,seq(from=1, to=length(.global.models))]
  	}

  	# go column by column
  	for (i in seq(1:length(.format.additional.lines.titles))) {
  		cat(.format.additional.lines.titles[i])
  		for (j in seq(1:length(.global.models))) {

  			if (is.na(.format.additional.lines.values[i,j])==FALSE) { cat(" & ",.format.additional.lines.values[i,j], sep="") }
  			else { cat(" & ", sep="") }
  		}
  		cat(" \\\\ \n")
  	}
  	.table.part.published[part.number] <<- TRUE

	
  }

  .print.table.statistic <-
  function(.global.var.name, .format.var.name, decimal.digits=.format.round.digits, part.string="", part.number=NULL, type.se=FALSE) {
	
  	# default values
  	report.df <- FALSE
    report.p.value <- FALSE
    significance.stars <- FALSE
    report.se <- FALSE
    report.tstat <- FALSE
    intelligent.df <- .format.intelligent.df

  	# reporting of df, p-value, significance stars, standard errors, t-stats
  	if (length(grep("(df)", part.string,fixed=TRUE))!=0) { report.df <- TRUE } 
  	if (length(grep("(p)", part.string,fixed=TRUE))!=0) { report.p.value <- TRUE } 
  	if (length(grep("*", part.string,fixed=TRUE))!=0) { significance.stars <- TRUE } 
  
  	if (length(grep("(se)", part.string,fixed=TRUE))!=0) { report.se <- TRUE }
  	if (length(grep("(t)", part.string,fixed=TRUE))!=0) { report.tstat <- TRUE }

  	# first for vectors (statistics without, say, degrees of freedom)
  	if (is.vector(.global.var.name) == TRUE) {
  		if (sum(!is.na(.global.var.name))!=0) {
  			cat (.format.var.name)
  			for (i in seq(1:length(.global.models))) {
  	 			if (!is.na(.global.var.name[i])) { 
             if (.format.dec.mark.align == TRUE) {
                cat(" & \\multicolumn{1}{c}{$",.iround(.global.var.name[i], decimal.digits),"$}", sep="")
             }
             else {
               cat(" & $",.iround(.global.var.name[i], decimal.digits),"$", sep="")
             }
          }
  	 			else { cat(" & ", sep="") }
  			}
  			cat(" \\\\ \n")
  			.table.part.published[part.number] <<- TRUE
  		}
  	}
  	else if ((is.matrix(.global.var.name) == TRUE) & (type.se == FALSE)) {     # for statistics that have degrees of freedom
  		if (sum(!is.na(as.vector(.global.var.name["statistic",])))!=0) {

	  		# intelligent df reporting (figure out whether only report it on left side, or also)
	  		report.df.left.column <- FALSE
			
	  		# whittle down unique values
	  		df.all.together <- NULL
	  		for (i in seq(1:length(.global.models))) {
	  			df.string <- ""
	  			for (j in seq(1:(nrow(.global.var.name)- 2))) {
	  				df.string <- paste(df.string,";",as.character(.global.var.name[paste("df",as.character(j),sep=""),i]),sep="")
	  			}
	  			df.all.together <- append(df.all.together, df.string)
	  		}
	  		# remove.na.r
	  		df.all.together.no.NA <- NULL
	  		for (i in seq(1:length(df.all.together))) {
	  			if (substr(df.all.together[i],1,3)!=";NA") { df.all.together.no.NA <- c(df.all.together.no.NA, df.all.together[i]) }
	  		}
	  		df.all.together.no.NA.unique <- sort(unique(df.all.together.no.NA))

	  		# put df on the left if only one unique df in the table, and not just one column w/ given df
	  		if (intelligent.df == TRUE) {
	  			if ((length(df.all.together.no.NA.unique)==1) & (length(df.all.together.no.NA)>=2)) { report.df.left.column <- TRUE }				
	  		}

  			# write down the line	
  			cat (.format.var.name)

  			# report df on left side w/ intelligent reporting
  			if (report.df.left.column == TRUE) {
  				if (report.df == TRUE) {

  					cat(" ",.format.df.left,sep="")
  					df.list <- unlist(strsplit(df.all.together.no.NA.unique[1],";"))

  					for (i in seq(from=2, to=length(df.list))) {
  						if (i>=3) { cat(.format.df.separator) }
  						cat(df.list[i],sep="")
  					}
  					cat(.format.df.right,sep="")
  				}
  			}
		
  			# now, go column by column
  			for (i in seq(1:length(.global.models))) {
  	 			if (!is.na(.global.var.name["statistic",i])) {
             
  	 			  if (.format.dec.mark.align==TRUE) {
  					  cat(" & \\multicolumn{1}{c}{$",.iround(.global.var.name["statistic",i], decimal.digits), sep="") 
  	 			  }
            else {
              cat(" & $",.iround(.global.var.name["statistic",i], decimal.digits), sep="")
            }

  					# significance stars
  					if ((significance.stars == TRUE) & (!is.na(.global.var.name["p-value",i]))) { .enter.significance.stars(.global.var.name["p-value",i], math.omit=TRUE) }

										
  					# degrees of freedom - only report by statistics if not in the left column already
  					if (report.df.left.column == FALSE) {
  						if ((report.df == TRUE) & (!is.na(.global.var.name["df1",i]))) {
  							cat(" ",.format.df.left,sep="")
  							for (j in seq(1:(nrow(.global.var.name)- 2))) {
  								if (!is.na(.global.var.name[paste("df",as.character(j),sep=""),i])) {
  									if (j>=2) { cat(.format.df.separator) }
  									cat(.global.var.name[paste("df",as.character(j),sep=""),i],sep="")
  								}
  							}
  							cat(.format.df.right,sep="")
  						}
  					}

  					# p-values
  					if ((report.p.value == TRUE) & (!is.na(.global.var.name["p-value",i]))) {
  						cat(" ",.format.p.value.left,sep="")
  						if (!is.na(.global.var.name[paste("df",as.character(j),sep=""),i])) { 
  							cat(.iround(.global.var.name["p-value",i],.format.round.digits, round.up.positive=TRUE),sep="") 
  						}
  						cat(.format.p.value.right,sep="")
  					}
            
            if (.format.dec.mark.align==TRUE) {
              cat("$}")  
            }
            else {
              cat("$")
            }
            
  				}
  				else { cat(" & ", sep="") }
  			}
  			cat(" \\\\ \n")			
  			.table.part.published[part.number] <<- TRUE
  		}
  	}
  	else if ((is.matrix(.global.var.name) == TRUE) & (type.se == TRUE)) {       # for statistics that have a standard error
  	  if (sum(!is.na(as.vector(.global.var.name["statistic",])))!=0) {
	    
  	    # write down the line	
  	    cat (.format.var.name)
	    
  	    # now, go column by column
  	    for (i in seq(1:length(.global.models))) {
  	      if (!is.na(.global.var.name["statistic",i])) { 

            if (.format.dec.mark.align == TRUE) {
              cat(" & \\multicolumn{1}{c}{$",.iround(.global.var.name["statistic",i], decimal.digits), sep="")  
            }
            else {
              cat(" & $",.iround(.global.var.name["statistic",i], decimal.digits), sep="")
            }
  	        
	        
  	        # significance stars
  	        if ((significance.stars == TRUE) & (!is.na(.global.var.name["p-value",i]))) { .enter.significance.stars(.global.var.name["p-value",i], math.omit=TRUE) }
	        
  	        # standard errors
  	        if ((report.se == TRUE) & (!is.na(.global.var.name["se",i]))) { cat(" ",.format.se.left,.iround(.global.var.name["se",i], decimal.digits),.format.se.right,sep="") }
          
  	        # t-statistics
  	        if ((report.tstat == TRUE) & (!is.na(.global.var.name["tstat",i]))) { cat(" ",.format.tstat.left, .iround(.global.var.name["tstat",i], decimal.digits),.format.tstat.right,sep="") }
          
  	        # p-values
  	        if ((report.tstat == TRUE) & (!is.na(.global.var.name["p-value",i]))) { cat(" ",.format.p.value.left,.iround(.global.var.name["p-value",i], decimal.digits),.format.p.value.right,sep="") }
            
            if (.format.dec.mark.align == TRUE) {
              cat("$}")
            }
            else {
              cat("$")
            }
  	      }
  	      else { cat(" & ", sep="") }
  	    }
  	    cat(" \\\\ \n")			
  	    .table.part.published[part.number] <<- TRUE
  	  }
  	}
  }

  .publish.table <-
  function() {

  	.table.info.comment()

  	# table header
	
  	.table.header()
  	.table.insert.space()

  	.table.part.published <<- as.vector(rep(NA, times=length(.format.table.parts)))    # to keep track what has been published (to deal intelligently with horizontal lines)
  	.publish.horizontal.line <<- TRUE   # should non-compulsory horizontal lines be published? (yes, if something else published since the previous line)

  	if (length(.format.table.parts)>=1) {
  		for (i in seq(1:length(.format.table.parts))) {
  			.publish.table.part(part=.format.table.parts[i], which.part.number=i)

  			if (.table.part.published[i]==TRUE) { .publish.horizontal.line <<- TRUE }
  			if ((.format.table.parts[i]=="-") | (.format.table.parts[i]=="-!") | (.format.table.parts[i]=="=") | (.format.table.parts[i]=="=!")) { .publish.horizontal.line <<- FALSE }
  		}
  	}

  	cat("\\normalsize \n")
  	cat("\\end{tabular} \n")
  	cat("\\end{table} \n")

	
  }

  .publish.table.part <-
  function(part, which.part.number) {

  	.table.part.published[which.part.number] <<- FALSE

  	# dependent variable label line
  	if (part=="dependent variable label") {
  		if (.format.dependent.variable.text.on == TRUE) { cat(" & \\multicolumn{",length(.global.models),"}{c}{",.format.dependent.variable.text, "} \\\\ \n", sep="")}
  		if (.format.dependent.variable.text.underline == TRUE) { cat("\\cline{2-",length(.global.models)+1,"} \n", sep="")}
  		.table.part.published[which.part.number] <<- TRUE
  	}

  	# dependent variables
  	else if (part=="dependent variables") {
  		.table.insert.space()
  		cat(.format.dependent.variables.text)
  		how.many.columns <- 0
      label.counter <- 0
    
  		for (i in seq(1:length(.global.models))) {
        if (is.null(.format.dep.var.labels)) { .format.dep.var.labels <<- NA }
  			how.many.columns <- how.many.columns + 1

  			# write down if next column has different dependent variable, or if end of columns
  			different.dependent.variable <- FALSE
  			if (i == length(.global.models)) {different.dependent.variable <- TRUE}
  			else if ((as.character(.global.dependent.variables[i])) != (as.character(.global.dependent.variables[i+1])))  {different.dependent.variable <- TRUE}

  			if (different.dependent.variable == TRUE) {
          label.counter <- label.counter + 1 
  		 		if (how.many.columns == 1) {
            if (.format.dec.mark.align==TRUE) {
  		 		    if (is.na(.format.dep.var.labels[label.counter])) {
  		 		      if (.format.dependent.variables.capitalize == TRUE) { cat(" & \\multicolumn{1}{c}{",.format.dependent.variables.left,toupper(as.character(.global.dependent.variables.written[i])),.format.dependent.variables.right,"}", sep="") }
  		 		      else { cat(" & \\multicolumn{1}{c}{",.format.dependent.variables.left,as.character(.global.dependent.variables.written[i]),.format.dependent.variables.right,"}", sep="") }
  		 		    }
  		 		    else { cat(" & \\multicolumn{1}{c}{",.format.dependent.variables.left,.format.dep.var.labels[label.counter],.format.dependent.variables.right,"}", sep="") }
            }
            else {
              if (is.na(.format.dep.var.labels[label.counter])) {
  					    if (.format.dependent.variables.capitalize == TRUE) { cat(" & ",.format.dependent.variables.left,toupper(as.character(.global.dependent.variables.written[i])),.format.dependent.variables.right, sep="") }
  					    else { cat(" & ",.format.dependent.variables.left,as.character(.global.dependent.variables.written[i]),.format.dependent.variables.right, sep="") }
              }
              else { cat(" & ",.format.dependent.variables.left,.format.dep.var.labels[label.counter],.format.dependent.variables.right, sep="") }
            }
  				}
  				else {
            if (is.na(.format.dep.var.labels[label.counter])) {
  					  if (.format.dependent.variables.capitalize == TRUE) {cat(" & \\multicolumn{",how.many.columns,"}{c}{",.format.dependent.variables.left,toupper(as.character(.global.dependent.variables.written[i])),.format.dependent.variables.right,"}", sep="")}
  					  else {cat(" & \\multicolumn{",how.many.columns,"}{c}{",.format.dependent.variables.left,as.character(.global.dependent.variables.written[i]),.format.dependent.variables.right,"}", sep="")}
            }
            else {cat(" & \\multicolumn{",how.many.columns,"}{c}{",.format.dependent.variables.left,.format.dep.var.labels[label.counter],.format.dependent.variables.right,"}", sep="")}
  				}

	  			how.many.columns <- 0
	  		}
	  	}
	  	cat(" \\\\ \n")

  		.table.part.published[which.part.number] <<- TRUE
  	}

  	# models
  	else if (part=="models")  {
     	   if ((.format.model.names.include==TRUE) & ((.format.models.skip.if.one == FALSE) | ((.format.models.skip.if.one == TRUE) & (length(unique(.global.models))>=2)))) {
		
  		.table.insert.space()
  		cat(.format.models.text)
 
  		# rename models based on .formatting preferences
  		renamed.global.models <- as.matrix(rbind(.global.models, rep("", times=length(.global.models))))
  		for (i in seq(1:length(.global.models))) {
  			for (j in seq(1:ncol(.format.model.names))) {
  				model.strsplit <- unlist(strsplit(.global.models[i], split="#"))
  				if (.global.models[i]==.format.model.names[1,j]) { 
  					renamed.global.models[1,i] <- .format.model.names[2,j] 
  					renamed.global.models[2,i] <- .format.model.names[3,j]
  				}
  				else if ((model.strsplit[1]=="glm()") | (model.strsplit[1]=="svyglm()") | (model.strsplit[1]=="gee()") | (model.strsplit[1]=="gam()")) {
  					if ( .format.model.function == TRUE ) { renamed.global.models[1,i] <- paste(substr(model.strsplit[1],1,nchar(model.strsplit[1])-2),": ", .format.model.family, model.strsplit[2], sep="") }
  					else { renamed.global.models[1,i] <- paste(.format.model.family, model.strsplit[2], sep="")}

  					renamed.global.models[2,i] <- paste(.format.model.link, model.strsplit[3], sep="")
  				}
  				else if ((model.strsplit[1]=="survreg()") | (model.strsplit[1]=="polr()")) {
  					if ( .format.model.function == TRUE ) { renamed.global.models[1,i] <- paste(substr(model.strsplit[1],1,nchar(model.strsplit[1])-2),": ", .format.model.dist, model.strsplit[2], sep="") }
  					else { renamed.global.models[1,i] <- paste(.format.model.dist, model.strsplit[2], sep="")}
  					renamed.global.models[2,i] <- ""
  				}
  			}
  		}

  		if (sum(renamed.global.models[2,]==rep("", times=length(.global.models)))==length(.global.models)) { how.many.model.rows <- 1}
  		else { how.many.model.rows <- 2 }

  		for (row in seq(from=1, to=how.many.model.rows)) {
  			how.many.columns <- 0
  			for (i in seq(1:length(.global.models))) {
  				how.many.columns <- how.many.columns + 1

  				# write down if next column has different dependent variable, or if end of columns
  				different.model <- FALSE
  				if (i == length(.global.models)) {different.model <- TRUE}
  				else if ((as.character(.global.models[i])) != (as.character(.global.models[i+1]))) {different.model <- TRUE}
  				else if ((as.character(.global.dependent.variables[i])) != (as.character(.global.dependent.variables[i+1]))) {different.model <- TRUE}   # subsume models under dependent variables

  				if (different.model == TRUE) {
  			 		if (how.many.columns == 1) {
              if (.format.dec.mark.align == TRUE) {  
  			 		    cat(" & \\multicolumn{1}{c}{",.format.models.left,as.character(renamed.global.models[row,i]),.format.models.right,"}", sep="")
              }
              else {
                cat(" & ",.format.models.left,as.character(renamed.global.models[row,i]),.format.models.right, sep="")
              }
  			 		}
  					else {cat(" & \\multicolumn{",how.many.columns,"}{c}{",.format.models.left,as.character(renamed.global.models[row,i]),.format.models.right,"}", sep="")}

  					how.many.columns <- 0
  				}
  			}
  			cat(" \\\\ \n")	
  		}
	

  		# underline models
  		if (.format.underline.models == TRUE) {
  			how.many.columns <- 0
  			for (i in seq(1:length(.global.models))) {
  				how.many.columns <- how.many.columns + 1

  				# underline if next column has different dependent variable, or if end of columns
  				different.model <- FALSE
  				if (i == length(.global.models)) {different.model <- TRUE}
  				else if ((as.character(.global.models[i])) != (as.character(.global.models[i+1])))  {different.model <- TRUE}
  				else if ((as.character(.global.dependent.variables[i])) != (as.character(.global.dependent.variables[i+1]))) {different.model <- TRUE}   # subsume models under dependent variables

  				if (different.model== TRUE) {
  					cat("\\cline{",(i-how.many.columns+1)+1,"-",i+1,"} ",sep="")
	
  					how.many.columns <- 0
  				}
  			}
  		cat("\n")
  		}
  		.table.part.published[which.part.number] <<- TRUE
            }
  	}

  	# numbers
  	else if (part=="numbers") {
      if ((.format.model.numbers == TRUE) & (length(.global.models)>1)) {
  		  .table.insert.space()
  		  cat(.format.numbers.text)
  		  for (i in seq(1:length(.global.models))) {
          if (.format.dec.mark.align==TRUE) {
  	 		    if (.format.numbers.roman == TRUE) { cat(" & \\multicolumn{1}{c}{",.format.numbers.left,.roman.numeral(i),.format.numbers.right,"}", sep="") }
  	 		    else { cat(" & \\multicolumn{1}{c}{",.format.numbers.left,i,.format.numbers.right,"}", sep="") }
          }
          else {
            if (.format.numbers.roman == TRUE) { cat(" & ",.format.numbers.left,.roman.numeral(i),.format.numbers.right, sep="") }
            else { cat(" & ",.format.numbers.left,i,.format.numbers.right, sep="") }
          }

		    }
		    cat("\\\\ \n")
  		  .table.part.published[which.part.number] <<- TRUE
      }
  	}

  	## coefficients
  	else if (part=="coefficients") {
  		# first, find the position of the intercept
  		intercept.position <- NULL
  		for (i in seq(1:length(.global.coefficient.variables))) {
  			if (.global.coefficient.variables[i] %in% .global.intercept.strings) { intercept.position <- i }
  		}
    		
  		.which.variable.label <<- 0
  		if (is.null(.format.covariate.labels)) { .format.covariate.labels <<- NA }
		
      # then, enter the coefficients
  		if ((!is.null(intercept.position)) & (.format.intercept.bottom==FALSE)) { .table.enter.coefficients(intercept.position) }
  		for (i in seq(1:length(.global.coefficient.variables))) {
  			if (!(.global.coefficient.variables[i] %in% .global.intercept.strings)) { .table.enter.coefficients(i) }
  		}

  		if ((!is.null(intercept.position)) & (.format.intercept.bottom==TRUE)) { .table.enter.coefficients(intercept.position) }

  		.table.part.published[which.part.number] <<- TRUE
  	}

  	# number of observations
  	else if (part=="N") { .print.table.statistic(.global.var.name=.global.N, .format.var.name=.format.N, decimal.digits=0, part.number=which.part.number) }

  	# fixed effects table
  	else if (part=="omit") {
      if ((!is.null(.format.omit.regexp)) & (!is.null(.format.omit.labels))) {
  		  .format.omit.table <<- matrix(.format.omit.no, nrow=length(.format.omit.regexp), ncol=length(.global.models)) 
  		  for (i in seq(1:length(.global.models))) {
  				for (j in seq(1:length(.format.omit.regexp))) {
  					for (k in seq(1:length(.global.coefficients))) {
  						if (length(grep(.format.omit.regexp[j], .global.coefficient.variables[k],fixed=FALSE))!=0) {
  							if (!is.na(.global.coefficients[k,i])) {
  								.format.omit.table[j,i] <<- .format.omit.yes
  							}
  						}
  					}
  				}
  			}
        for (i in seq(1:length(.format.omit.regexp))) {
  			  cat (.format.omit.labels[i])
  			  for (j in seq(1:length(.global.models))) {
            if (.format.dec.mark.align == TRUE) {
  		 		    cat(" & \\multicolumn{1}{c}{",.format.omit.table[i,j],"}", sep="")
            }
            else {
              cat(" & ",.format.omit.table[i,j], sep="")
            }
    			}
  			  cat(" \\\\ \n")
  			}
  			.table.part.published[which.part.number] <<- TRUE
      }
  	}

  	# R-squared
  	else if (part=="R-squared") {	.print.table.statistic(.global.var.name=.global.R2, .format.var.name=.format.R2, part.number=which.part.number) }
  
  	# max R-squared
  	else if (part=="max R-squared") {	.print.table.statistic(.global.var.name=.global.max.R2, .format.var.name=.format.max.R2, part.number=which.part.number) }

  	# adjusted R-squared
  	else if (part=="adjusted R-squared") { .print.table.statistic(.global.var.name=.global.adj.R2, .format.var.name=.format.adj.R2, part.number=which.part.number) }

  	# log likelihood
  	else if (part=="log likelihood") { .print.table.statistic(.global.var.name=.global.LL, .format.var.name=.format.LL, part.number=which.part.number) }

  	# Akaike Information Criterion (AIC)
  	else if (part=="AIC") { .print.table.statistic(.global.var.name=.global.AIC, .format.var.name=.format.AIC, part.number=which.part.number) }
  
  	# Bayesian Information Criterion (BIC)
  	else if (part=="BIC") { .print.table.statistic(.global.var.name=.global.BIC, .format.var.name=.format.BIC, part.number=which.part.number) }
  	
  	# Scale Parameter
  	else if (part=="scale") { .print.table.statistic(.global.var.name=.global.scale, .format.var.name=.format.scale, part.number=which.part.number) }
  
  	# UBRE
  	else if (part=="UBRE") { .print.table.statistic(.global.var.name=.global.UBRE, .format.var.name=.format.UBRE, part.number=which.part.number) }
  
  	# sigma2
  	else if (part=="sigma2") { .print.table.statistic(.global.var.name=.global.sigma2, .format.var.name=.format.sigma2, part.number=which.part.number) }

  	## with degrees of freedom

	  # residual standard error (sigma); standard error of the regression
  	else if (substr(part,1,nchar("SER"))=="SER") { .print.table.statistic(.global.var.name=.global.SER, .format.var.name=.format.SER, part.string=part, part.number=which.part.number) }
	
  	# F-statitistic
  	else if (substr(part,1,nchar("F statistic"))=="F statistic") { .print.table.statistic(.global.var.name=.global.F.stat, .format.var.name=.format.F.stat, part.string=part, part.number=which.part.number) }
  
  	# theta
  	else if (substr(part,1,nchar("theta"))=="theta") { .print.table.statistic(.global.var.name=.global.theta, .format.var.name=.format.theta, part.string=part, part.number=which.part.number, type.se=TRUE) }
  
    # Chi-squared
  	else if (substr(part,1,nchar("chi2"))=="chi2") { .print.table.statistic(.global.var.name=.global.chi.stat, .format.var.name=.format.chi.stat, part.string=part, part.number=which.part.number) }
  
  	# Wald Test
  	else if (substr(part,1,nchar("Wald"))=="Wald") { .print.table.statistic(.global.var.name=.global.wald.stat, .format.var.name=.format.wald.stat, part.string=part, part.number=which.part.number) }
  
  	# LR Test
  	else if (substr(part,1,nchar("LR"))=="LR") { .print.table.statistic(.global.var.name=.global.lr.stat, .format.var.name=.format.lr.stat, part.string=part, part.number=which.part.number) }
  
  	# Score (Logrank) Test
  	else if (substr(part,1,nchar("logrank"))=="logrank") { .print.table.statistic(.global.var.name=.global.logrank.stat, .format.var.name=.format.logrank.stat, part.string=part, part.number=which.part.number) }

  	# null deviance
  	else if (substr(part,1,nchar("null deviance"))=="null deviance") { .print.table.statistic(.global.var.name=.global.null.deviance, .format.var.name=.format.null.deviance, part.string=part, part.number=which.part.number) }

  	# residual deviance
  	else if (substr(part,1,nchar("residual deviance"))=="residual deviance") { .print.table.statistic(.global.var.name=.global.residual.deviance, .format.var.name=.format.residual.deviance, part.string=part, part.number=which.part.number) }

  	##

  	# single horizontal line, no matter what
  	else if (part=="-!") {
  		cat("\\hline ")
  		.table.insert.space()
  		cat(" \n") 
  		.table.part.published[which.part.number] <<- TRUE
  	}

  	# single horizontal line, optional
  	else if (part=="-") {
  		if (.publish.horizontal.line==TRUE) {
  			cat("\\hline ")
  			.table.insert.space()
  			cat(" \n") 
  			.table.part.published[which.part.number] <<- TRUE
  		}
  	}

  	# double horizontal line, no matter what
  	else if (part=="=!") {
  		cat("\\hline \n") 
  		cat("\\hline ")
  		.table.insert.space()
  		cat(" \n")
  		.table.part.published[which.part.number] <<- TRUE
  	}

  	# double horizontal line
  	else if (part=="=") {
  		if (.publish.horizontal.line==TRUE) {
  			cat("\\hline \n") 
  			cat("\\hline ")
  			.table.insert.space()
  			cat(" \n") 
  			.table.part.published[which.part.number] <<- TRUE
  		}
  	}

  	# notes
  	else if (part=="notes") {
  		if (.format.note != "") { cat(.format.note) }
  		for (i in seq(1:length(.format.note.content))) {
        
        # print individual notes
  			if (.format.note == "") { cat("\\multicolumn{",length(.global.models)+1,"}{",.format.note.alignment,"}{",.format.note.content[i],"} \\\\ \n", sep="") }
  			else { cat(" & \\multicolumn{",length(.global.models),"}{",.format.note.alignment,"}{",.format.note.content[i],"} \\\\ \n", sep="") }
  		}
  		.table.part.published[which.part.number] <<- TRUE
  	}	

  	# empty line
  	else if (part==" ") {
  		.table.empty.line();
  		.table.part.published[which.part.number] <<- TRUE
  	}

  	# additional lines
  	else if (part=="additional") { .print.additional.lines(part.number=which.part.number) }
  }

  .r.squared <-
  function(object.name) {

  	model.name <- .get.model.name(object.name)

  	if (!(model.name %in% c("arima", "nlmer", "glmer", "lmer"))) {
      if (model.name %in% c("plm")) {
        return(as.vector(.summary.object$r.squared["rsq"]))
      }
      else if (model.name %in% c("betareg")) {
        return(as.vector(.summary.object$pseudo.r.squared))
      }
  		else if (!is.null(.summary.object$r.squared)) {
  			return(as.vector(.summary.object$r.squared)) 
  		}
  		else if (model.name %in% c("coxph", "clogit")) {
  			return(as.vector(.summary.object$rsq[1]))
  		}
      else if (model.name %in% c("pmg")) {
        return(as.vector(.summary.object$rsqr))
      }
  	}
  	return(NA)
  }
  
  .remove.special.chars <-
  function(s) {
    
    if (!is.character(s)) { s.out <- as.character(s) }
    else { s.out <- s }
    
    # this has to go first
    s.out <- gsub("\\","\\textbackslash ",s.out,fixed=TRUE)
    
    # basic special characters
    s.out <- gsub("_","\\_",s.out,fixed=TRUE)
    s.out <- gsub("^","\\^",s.out,fixed=TRUE)
    s.out <- gsub("#","\\#",s.out,fixed=TRUE)
    s.out <- gsub("~","\\~",s.out,fixed=TRUE)
    s.out <- gsub("{","\\{",s.out,fixed=TRUE)
    s.out <- gsub("}","\\}",s.out,fixed=TRUE)    
    s.out <- gsub("%","\\%",s.out,fixed=TRUE)
    
    # pre-defined text-mode commands (add more?)
    s.out <- gsub("*","\\textasteriskcentered ",s.out,fixed=TRUE)
    s.out <- gsub("|","\\textbar ",s.out,fixed=TRUE)
    s.out <- gsub(">","\\textgreater ",s.out,fixed=TRUE)
    s.out <- gsub("<","\\textless ",s.out,fixed=TRUE)
        
    return(s.out)
  }
  

  .residual.deviance <-
  function(object.name) {
  	residual.deviance.output <- as.vector(rep(NA,times=3))

  	model.name <- .get.model.name(object.name)

  	if (!(model.name %in% c("arima", "multinom","lmer","glmer","nlmer"))) {
  		if (!is.null(.summary.object$deviance)) {
  			residual.deviance.value <- suppressMessages(.summary.object$deviance)
  			df.value <- object.name$df.residual
  			residual.deviance.output <- as.vector(c(residual.deviance.value, df.value, NA))
  		}
      else if (!is.null(object.name$deviance)) {
  		  residual.deviance.value <- object.name$deviance
  		  df.value <- object.name$df.residual
  		  residual.deviance.output <- as.vector(c(residual.deviance.value, df.value, NA))
  		}
  	}

  	names(residual.deviance.output) <- c("statistic","df1","p-value")
  	return(cbind(residual.deviance.output))
  }

  .roman.numeral <-
  function(regular.number) {

  	# unique representation only for integers between 1 and 3899
  	if ((regular.number < 1) | (regular.number > 3899)) {
  		return(NULL)
  	}
  	else {
  		roman.output <- ""
  		number.remaining <- regular.number

  		while (number.remaining > 999) {
  			roman.output <- paste(roman.output, "M", sep="")
  			number.remaining <- number.remaining - 1000
  		}

  		if (number.remaining > 899) {
  			roman.output <- paste(roman.output, "CM", sep="")
  			number.remaining <- number.remaining - 900
  		}

  		if (number.remaining > 499) {
  			roman.output <- paste(roman.output, "D", sep="")
  			number.remaining <- number.remaining - 500
  		}

  		if (number.remaining > 399) {
  			roman.output <- paste(roman.output, "CD", sep="")
  			number.remaining <- number.remaining - 400
  		}

  		if (number.remaining > 399) {
  			roman.output <- paste(roman.output, "D", sep="")
  			number.remaining <- number.remaining - 400
  		}

  		while (number.remaining > 99) {
  			roman.output <- paste(roman.output, "C", sep="")
  			number.remaining <- number.remaining - 100
  		}

  		if (number.remaining > 89) {
  			roman.output <- paste(roman.output, "XC", sep="")
  			number.remaining <- number.remaining - 90
  		}

  		if (number.remaining > 49) {
  			roman.output <- paste(roman.output, "L", sep="")
  			number.remaining <- number.remaining - 50
  		}

  		if (number.remaining > 39) {
  			roman.output <- paste(roman.output, "XL", sep="")
  			number.remaining <- number.remaining - 40
  		}

  		while (number.remaining > 9) {
  			roman.output <- paste(roman.output, "X", sep="")
  			number.remaining <- number.remaining - 10
  		}

  		if (number.remaining > 8) {
  			roman.output <- paste(roman.output, "IX", sep="")
  			number.remaining <- number.remaining - 9
  		}

  		if (number.remaining > 4) {
  			roman.output <- paste(roman.output, "V", sep="")
  			number.remaining <- number.remaining - 5
  		}

  		if (number.remaining > 3) {
  			roman.output <- paste(roman.output, "IV", sep="")
  			number.remaining <- number.remaining - 4
  		}

  		if (number.remaining > 3) {
  			roman.output <- paste(roman.output, "IV", sep="")
  			number.remaining <- number.remaining - 4
  		}

  		while (number.remaining > 0) {
  			roman.output <- paste(roman.output, "I", sep="")
  			number.remaining <- number.remaining - 1
  		}

  		return(roman.output)
  	}
  }

  .SER <-
  function(object.name) {
  	SER.output <- as.vector(rep(NA,times=3))

  	model.name <- .get.model.name(object.name)

  	if (!(model.name %in% c("arima","lmer","glmer","nlmer","gls"))) {
  		if (!is.null(suppressMessages(.summary.object$sigma))) {
  			sigma.value <-suppressMessages(.summary.object$sigma)
        if (model.name %in% c("rlm")) {
          df.residual.value <- .summary.object$df[2]
        } 
        else {
  			  df.residual.value <- object.name$df.residual
        }
  			SER.output <- as.vector(c(sigma.value, df.residual.value, NA))
  		}
  	}

  	names(SER.output) <- c("statistic","df1","p-value")
  	return(cbind(SER.output))
  }

  .sg.add.line <-
  function(new.line.title="", new.line.values=NA) {
  	if ((is.null(new.line.title)==FALSE) & (is.null(new.line.values)==FALSE)) {
  		if (is.null(.format.additional.lines.values) == FALSE) {
  			if (length(new.line.values) < ncol(.format.additional.lines.values)) {
  					new.line.values <- c(new.line.values,as.vector(rep(NA, times=(ncol(.format.additional.lines.values)-length(new.line.values)))))
  			}
  			else if (length(new.line.values) > ncol(.format.additional.lines.values)) {
  				for (i in seq(1:(length(new.line.values) - ncol(.format.additional.lines.values)))) {
  					.format.additional.lines.values <<- cbind(.format.additional.lines.values,as.vector(rep(NA, times=nrow(.format.additional.lines.values))))
  				}
  			}
  		}
  		.format.additional.lines.titles <<- c(.format.additional.lines.titles, new.line.title[1])
  		.format.additional.lines.values <<- rbind(.format.additional.lines.values, new.line.values)
  	}
  }

  .stargazer.reg.table <-
  function(...) {

  	list.of.models <- as.list(list(...))
  	how.many.models <- length(list.of.models)
    
    # find how many models user wants to customize
    max.user <- max(length(coef),length(se),length(t),length(p))
  
  	if (how.many.models >= 1) {
  		suppressMessages(.new.table(list.of.models[[1]], user.coef=coef[[1]], user.se=se[[1]], user.t=t[[1]], user.p=p[[1]], auto.t=t.auto, auto.p=p.auto))
  		if (how.many.models >= 2) {
  			for (i in seq(from = 2,to = how.many.models)) { 
          if (i <= max.user) {
            suppressMessages(.add.model(list.of.models[[i]], user.coef=coef[[i]], user.se=se[[i]], user.t=t[[i]], user.p=p[[i]], auto.t=t.auto, auto.p=p.auto)) 
          }
          else {
            suppressMessages(.add.model(list.of.models[[i]], user.coef=NULL, user.se=NULL, user.t=NULL, user.p=NULL, auto.t=t.auto, auto.p=p.auto))
          }
  			}
  		}
  		suppressMessages(.publish.table())
  	}
  }
  
  .data.frame.table.header <-
    function(object) {
      cat("\\begin{table}[htb] \\centering \n",sep="")
      cat("  \\caption{", .format.title, "} \n",sep="")   
      cat("  \\label{", .format.label, "} \n",sep="")
      cat("\\footnotesize \n",sep="")
      cat("\n")
      
      .formatting.alignment <- "@{\\extracolsep{5pt}} "
      for (i in seq(1:(length(names(object))))) {
        if (.format.dec.mark.align == FALSE) {
          .formatting.alignment <- paste(.formatting.alignment, "c ", sep="")
        }
        else {
          .formatting.alignment <- paste(.formatting.alignment, "D{", .format.decimal.character,"}{", .format.decimal.character,"}{-", .format.s.round.digits,"} ", sep="")
        }
      }
      #
      
      
      cat("\\begin{tabular}{",.formatting.alignment,"} \n",sep="")
    }
  
  .stargazer.data.frame.table <-
  function(object) {
    if (length(names(object)) < 1) {
      cat("% Error: Data frame columns do not have any names.\n")
    }
    else if ((nrow(object) < 1) | (ncol(object) < 1)) {
      cat("% Error: Data frame must have at least one row and one column.\n")
    }
    else {
      .table.info.comment()
      
      #create table header
      .data.frame.table.header(object)
      .table.insert.space()
      
      for (i in seq(1:length(.format.s.stat.parts))) {
        .data.frame.table.part(object,.format.s.stat.parts[i])
      }
      
      cat("\\normalsize \n")
      cat("\\end{tabular} \n")
      cat("\\end{table} \n")
    }
  }
  
  .data.frame.table.part <-
  function(object, part) {
    if (part=="stat names") {
      
      x.which <- 0
      
      if (is.null(.format.covariate.labels)) { .format.covariate.labels <<- NA }
      
      for (x in seq(1:length(names(object)))) {
        
        omitted <- FALSE
        
        
        if (!is.null(.format.omit.regexp)) {
          for (j in seq(1:length(.format.omit.regexp))) {
            if (length(grep(.format.omit.regexp[j], names(object)[x], fixed=FALSE))!=0) { omitted <- TRUE  }
          }
        }
        
        if (omitted == FALSE) {
          
          x.which <- x.which + 1
          
          if (x >= 2) { cat(" & ", sep="")}
          
          # if underscore or ^ in variable name, then insert an escape \ before it
          name.printed <- .remove.special.chars(names(object)[x])
                
          if (is.na(.format.covariate.labels[x.which])) {
            if (.format.coefficient.variables.capitalize == TRUE) { name.printed <- toupper(name.printed) }
          }
          else { name.printed <- .format.covariate.labels[x.which] }
      
        
          if (.format.dec.mark.align==TRUE) {
            cat("\\multicolumn{1}{c}{",.format.s.coefficient.variables.left, name.printed,.format.s.coefficient.variables.right,"}", sep="")  
          }
          else {
            cat(.format.s.coefficient.variables.left, name.printed,.format.s.coefficient.variables.right, sep="")  
          }
        }
        
      }
      cat(" \\\\ \n")
    }
    
    if (substr(part,1,10)=="statistics") {
      for (y in seq(1:nrow(object))) {
        for (x in seq(1:length(names(object)))) {
          
          omitted <- FALSE
          
          if (!is.null(.format.omit.regexp)) {
            for (j in seq(1:length(.format.omit.regexp))) {
              if (length(grep(.format.omit.regexp[j], names(object)[x], fixed=FALSE))!=0) { omitted <- TRUE  }
            }
          }
          
          if (omitted == FALSE) {     
            if (x >= 2) { cat(" & ", sep="") }
            
            how.much.to.round <- .format.round.digits
            if (is.numeric(object[y,x])) {
              
              if (.is.all.integers(object[y,x])) { .how.much.to.round <- 0 }
              
              rounded.object <- .iround(object[y,x],.how.much.to.round)
              
              if (.format.dec.mark.align==TRUE) {
                cat(rounded.object, sep="")  
              }
              else {
                cat("$", rounded.object, "$",sep="")  
              }
            }
            else {
              adjusted.object <- .remove.special.chars(object[y, x])
              if (.format.dec.mark.align==TRUE) {
                cat("\\multicolumn{1}{c}{", adjusted.object, "}", sep="")  
              }
              else {
                cat(adjusted.object, sep="")  
              }
              
            }
            
          }
          
          
            
        }
        # add empty lines
        how.many.empty.lines <- as.numeric(substr(part,11,nchar(part)))
        if (is.na(how.many.empty.lines)) { how.many.empty.lines <- 1 } 
        
        for (j in seq(1:how.many.empty.lines)) {
          cat(" \\\\ \n") 
        }
      }
    
    }
        
    
    # notes
    else if (part=="notes") {
      if (.format.s.note != "") cat(.format.s.note)
      for (i in seq(1:length(.format.s.note.content))) {
        if (.format.s.note == "") { cat("\\multicolumn{",length(.format.s.statistics.list)+1,"}{",.format.s.note.alignment,"}{",.format.s.note.content[i],"} \\\\ \n", sep="") }
        else { cat(" & \\multicolumn{",length(.format.s.statistics.list),"}{",.format.s.note.alignment,"}{",.format.s.note.content[i],"} \\\\ \n", sep="") }
      }
    }	
    
    
    # empty line
    else if (part==" ") {
      .table.empty.line()
    }
    
    # horizontal line
    else if (part=="-!") {
      cat("\\hline ")
      .table.insert.space()
      cat(" \n")
    }
    
    # double horizontal line
    else if (part=="=!") {
      cat("\\hline \n") 
      cat("\\hline ")
      .table.insert.space()
      cat(" \n")
    }
    
  }

  .stargazer.summ.stat.table <-
  function(object) {
  
    if (length(names(object)) < 1) {
      cat("% Error: Data frame columns do not have any names.\n")
    }
    else if ((nrow(object) < 1) | (ncol(object) < 1)) {
      cat("% Error: Data frame must have at least one row and one column.\n")
    }
    else {
      
      .table.info.comment()
      
  	  # create table header
  	  .summ.stat.table.header()
  	  .table.insert.space()

  	  for (i in seq(1:length(.format.s.stat.parts))) {
    		.summ.stat.table.part(object,.format.s.stat.parts[i])
    	}

    	cat("\\normalsize \n")
    	cat("\\end{tabular} \n")
    	cat("\\end{table} \n")
    }
  }

  .summ.stat.publish.statistic <-
  function(object, which.variable, which.statistic) {
    
  	if ((is.numeric(object[,which.variable]) == TRUE) | ((is.logical(object[,which.variable])) & (.format.summ.logical==TRUE)))  {
      
      if ((is.logical(object[,which.variable])) & (.format.summ.logical==TRUE)) {
        temp.var <- rep(NA, time=length(object[,which.variable]))
        temp.var[object[,which.variable]==TRUE] <- 1
        temp.var[object[,which.variable]==FALSE] <- 0
      }
      else {
        temp.var <- object[,which.variable]
      }
      
  		if (which.statistic == "N") {
  			return(.iround(sum(!is.na(temp.var)), 0))
  		}
  		else if (which.statistic == "nmiss") {
  			return(.iround(sum(is.na(temp.var)), 0))
  		}
  		else if (which.statistic == "mean") {
  			return(.iround(mean(temp.var, na.rm=TRUE), .format.s.round.digits))
  		}
  		else if (which.statistic == "median") {
  			median.value <- median(temp.var, na.rm=TRUE)
     
  			if (.is.all.integers(temp.var) == FALSE) { how.much.to.round <- .format.s.round.digits }
  			else { 
  				if (.is.all.integers(median.value) == TRUE) { how.much.to.round <- 0 }
  				else { how.much.to.round <- 1 }
	  		}

	  		return(.iround(median.value, how.much.to.round))
  		}
  		else if (which.statistic == "sd") {
  			return(.iround(sd(temp.var, na.rm=TRUE), .format.s.round.digits))
  		}
  		else if (which.statistic == "min") {
  			if (.is.all.integers(temp.var) == FALSE) { how.much.to.round <- .format.s.round.digits }
  			else { how.much.to.round <- 0 }

  			return(.iround(min(temp.var, na.rm=TRUE), how.much.to.round))
  		}
  		else if (which.statistic == "max") {
  			if (.is.all.integers(temp.var) == FALSE) { how.much.to.round <- .format.s.round.digits }
  			else { how.much.to.round <- 0 }

  			return(.iround(max(temp.var, na.rm=TRUE), how.much.to.round))
  		}
  		else if (which.statistic == "mad") {
  			return(.iround(mad(temp.var, na.rm=TRUE), .format.s.round.digits))
  		}
  		else if (substr(which.statistic,1,1) == "p") {
      
        percentile.value <- quantile(temp.var, as.numeric(substr(which.statistic,2,nchar(which.statistic))) / 100, na.rm=TRUE)
		  
        if (.is.all.integers(temp.var) == FALSE) { how.much.to.round <- .format.s.round.digits }
  		  else { 
  		    if (.is.all.integers(percentile.value) == TRUE) { how.much.to.round <- 0 }
  		    else { how.much.to.round <- 1 }
  		  }
		  
  		  return(.iround(percentile.value, how.much.to.round))
  		}
  	}
  	else { return(NA) }
  }

  .summ.stat.table.header <-
  function() {
    cat("\\begin{table}[htb] \\centering \n",sep="")
    cat("  \\caption{", .format.title, "} \n",sep="")
    cat("  \\label{", .format.label, "} \n",sep="")
    cat("\\footnotesize \n",sep="")
    cat("\n")

   #
   .formatting.alignment <- "@{\\extracolsep{5pt}}l "
   for (i in seq(1:(length(.format.s.statistics.list)))) {
    	if (.format.dec.mark.align == FALSE) {
        .formatting.alignment <- paste(.formatting.alignment, "c ", sep="")
    	}
      else {
        .formatting.alignment <- paste(.formatting.alignment, "D{", .format.decimal.character,"}{", .format.decimal.character,"}{-", .format.s.round.digits,"} ", sep="")
      }
    }
    #

   cat("\\begin{tabular}{",.formatting.alignment,"} \n",sep="")
  }

  .summ.stat.table.part <-
  function(object, part) {

  	if (part=="stat names") {
  		cat(.format.s.statistics.names.label, sep="")
      if (length(.format.s.statistics.list)>=1) {
  		  for (i in seq(1:length(.format.s.statistics.list))) {
  			  for (j in seq(1:ncol(.format.s.statistics.names))) {
  				  if ((substr(.format.s.statistics.list[i],1,1)=="p") & (substr(.format.s.statistics.list[i],1,1)==.format.s.statistics.names[1,j])) {
  					  cat(" & \\multicolumn{1}{c}{", .format.s.statistics.names.left, sub("!", substr(.format.s.statistics.list[i],2,nchar(.format.s.statistics.list[i])), .format.s.statistics.names[2,j], ignore.case =FALSE, fixed=TRUE), .format.s.statistics.names.right,"}", sep="")
  				  }
  				  else if (.format.s.statistics.list[i]==.format.s.statistics.names[1,j]) {
  				  	cat(" & \\multicolumn{1}{c}{", .format.s.statistics.names.left, .format.s.statistics.names[2,j], .format.s.statistics.names.right, "}", sep="")
  			  	}
  			  }
  		  }
      }
  		cat(" \\\\ \n")
  	}

  	if (substr(part,1,10)=="statistics") {
  	  if (is.null(.format.covariate.labels)) { .format.covariate.labels <<- NA }
	  
      i.label <- 0
  		for (i in seq(1:length(names(object)))) {
      
		        
  			if ((is.numeric(object[,i]) == TRUE) | (is.logical(object[,i]) & (.format.summ.logical==TRUE))) {
        
          # skip all of this if omitted based on regular expression
          omitted <- FALSE
        
         if (!is.null(.format.omit.regexp)) {
            for (j in seq(1:length(.format.omit.regexp))) {
              if (length(grep(.format.omit.regexp[j], names(object)[i], fixed=FALSE))!=0) { omitted <- TRUE	}
            }
          }
        
          if (omitted==FALSE) {
            i.label <- i.label + 1
            
            # if underscore in variable name, then insert an escape \ before it
            name.printed <- .remove.special.chars(names(object)[i])
            
            if (is.na(.format.covariate.labels[i.label])) { 
  				    if ( .format.s.coefficient.variables.capitalize == TRUE) { cat(.format.s.coefficient.variables.left, toupper(name.printed), .format.s.coefficient.variables.right, sep="") }
  				    else { cat(.format.s.coefficient.variables.left, name.printed, .format.s.coefficient.variables.right, sep="") }
            }
            else { cat(.format.s.coefficient.variables.left, .format.covariate.labels[i.label], .format.s.coefficient.variables.right, sep="") }

            if (length(.format.s.statistics.list)>=1) {
	  			    for (j in seq(1:length(.format.s.statistics.list))) {
                
                # if aligning decimal marks, need to use multicolumn for anything w/o decimal mark
                if (.format.dec.mark.align == FALSE) {   # not aligning
                  cat(" & ", .summ.stat.publish.statistic(object, i, .format.s.statistics.list[j]), sep="")  
                }
                else {     # aligning
                  if (.is.all.integers(.summ.stat.publish.statistic(object, i, .format.s.statistics.list[j]))) {
                    cat(" & \\multicolumn{1}{c}{", .summ.stat.publish.statistic(object, i, .format.s.statistics.list[j]),"}", sep="")
                  }
                  else {
                    cat(" & ", .summ.stat.publish.statistic(object, i, .format.s.statistics.list[j]), sep="")
                  }
                }
	  	 			    
	  			    }
            }
            
	
  				  # add empty lines
  				  how.many.empty.lines <- as.numeric(substr(part,11,nchar(part)))
  				  if (is.na(how.many.empty.lines)) { how.many.empty.lines <- 1 } 
				
  				  for (j in seq(1:how.many.empty.lines)) {
  					  cat(" \\\\ \n")
  				  }
  			  }
		
  			}
  		}
  	}

  	# notes
  	else if (part=="notes") {
  		if (.format.s.note != "") cat(.format.s.note)
  		for (i in seq(1:length(.format.s.note.content))) {
  			if (.format.s.note == "") { cat("\\multicolumn{",length(.format.s.statistics.list)+1,"}{",.format.s.note.alignment,"}{",.format.s.note.content[i],"} \\\\ \n", sep="") }
  			else { cat(" & \\multicolumn{",length(.format.s.statistics.list),"}{",.format.s.note.alignment,"}{",.format.s.note.content[i],"} \\\\ \n", sep="") }
  		}
  	}	

	
  	# empty line
  	else if (part==" ") {
  		.table.empty.line()
  	}

  	# horizontal line
  	else if (part=="-!") {
  		cat("\\hline ")
  		.table.insert.space()
  		cat(" \n")
  	}

  	# double horizontal line
  	else if (part=="=!") {
  		cat("\\hline \n") 
  		cat("\\hline ")
  		.table.insert.space()
  		cat(" \n")
  	}
  }

  .table.empty.line <-
  function() {
  	cat(" ")
  	for (i in seq(1:length(.global.models))) {
  		cat("& ")
  	}
  	cat("\\\\ \n")
  }

  .table.enter.coefficients <-
  function(which.variable) {
    
    if (which.variable > length(.global.coefficients)) {
      return();
    }

  	local.coefficient.var.name <- .global.coefficient.variables[which.variable]

  	#skip all of this if omitted based on regular expression
  	omitted <- FALSE
  	if (!is.null(.format.omit.regexp)) {
  		for (i in seq(1:length(.format.omit.regexp))) {
  			if (length(grep(.format.omit.regexp[i], local.coefficient.var.name, fixed=FALSE))!=0) { omitted <- TRUE	}
  		}
  	}

  	if (omitted == FALSE) {
    
      .which.variable.label <<- .which.variable.label + 1

  		# remove final -TRUE (added by Zelig) from dummy variables
  		if (substr(local.coefficient.var.name, nchar(local.coefficient.var.name)-3, nchar(local.coefficient.var.name)) == "TRUE") {

  			### only remove TRUE if added by Zelig, rather than pre-existing in the formula name
  			if (length(grep(local.coefficient.var.name, .global.formulas.rhs,fixed=TRUE))==0) {    
  				local.coefficient.var.name <- substr(local.coefficient.var.name, 1, nchar(local.coefficient.var.name)-4)
  			}
  		}
      
      # remove everything before and including he last dollar sign from variable name
      temp <- strsplit(local.coefficient.var.name,"$",fixed=TRUE)
      local.coefficient.var.name <- temp[[1]][length(temp[[1]])]
      
      # if underscore or ^ in variable name, then insert an escape \ before it
      local.coefficient.var.name <- .remove.special.chars(local.coefficient.var.name)
      
  		if (length(.format.coefficient.table.parts)>=1) {
  			for (i in seq(1:length(.format.coefficient.table.parts))) {
  				.coefficient.table.part(part=.format.coefficient.table.parts[i], which.variable, variable.name=local.coefficient.var.name)
  			}
  		}
  	}
  }

  .table.header <-
  function() {
    cat("\\begin{table}[htb] \\centering \n",sep="")
    cat("  \\caption{", .format.title, "} \n",sep="")   
    cat("  \\label{", .format.label, "} \n",sep="")
    cat("\\footnotesize \n",sep="")
    cat("\n")

   #
  .formatting.alignment <- "@{\\extracolsep{5pt}}l"
  for (i in seq(1:length(.global.models))) {
  	if (.format.dec.mark.align==FALSE) {
      .formatting.alignment <- paste(.formatting.alignment, "c", sep="")
  	}
    else {
  	  .formatting.alignment <- paste(.formatting.alignment, "D{", .format.decimal.character,"}{", .format.decimal.character,"}{-", .format.round.digits,"} ", sep="")
    }
  }
  #

  cat("\\begin{tabular}{",.formatting.alignment,"} \n",sep="")
  }

  .table.info.comment <-
  function() {
    cat("\n")
    cat("% Table created by ", .global.package.name, " v.", .global.package.version, " by ", .global.package.author.name, ", ", .global.package.author.affiliation, ". E-mail: ", .global.package.author.email, "\n", sep="")  
    cat("% Date and time:", format(Sys.time(), "%a, %b %d, %Y - %X"))
    cat("\n")
    if (.format.dec.mark.align==TRUE) {
      cat("% Requires \\usepackage{dcolumn} in LaTeX preamble.")
      cat("\n")
    }
  }

  .table.insert.space <-
  function() {
  	cat("\\\\[",.format.space.size,"]",sep="")
  }

  .trim <-
  function (x) gsub("^\\s+|\\s+$", "", x)

  .wald.stat <-
  function(object.name) {
    wald.output <- as.vector(rep(NA,times=3))
  
    model.name <- .get.model.name(object.name)
  
    if (!(model.name %in% c("arima", "ivreg","lmer","glmer","nlmer"))) {
      if (!is.null(.summary.object$waldtest)) {
        wald.value <- suppressMessages(.summary.object$waldtest[1])
        df.value <- suppressMessages(.summary.object$waldtest[2])
        wald.p.value <- suppressMessages(.summary.object$waldtest[3])
        wald.output <- as.vector(c(wald.value, df.value, wald.p.value))
      }
      else if (model.name %in% c("tobit(AER)")) {
        wald.value <- .summary.object$wald
        df.value <- .summary.object$df - .summary.object$idf
        wald.p.value <- pchisq(wald.value, df.value, lower.tail=FALSE)
        wald.output <- as.vector(c(wald.value, df.value, wald.p.value))
        
      }
    
    }
  
    names(wald.output) <- c("statistic","df1","p-value")
    return(cbind(wald.output))
  }

  .get.coefficients.1 <-
  function(object.name, user.given=NULL) {
    
    if (!is.null(user.given)) { return(user.given) }

  	model.name <- .get.model.name(object.name)
	
  	if (model.name %in% c("ls", "normal", "logit", "probit", "relogit", "poisson", "negbin", "normal.survey", "poisson.survey", "probit.survey", "logit.survey", "gamma", "gamma.survey",
     				    "cloglog.net", "gamma.net", "logit.net", "probit.net", "glm()", "svyglm()", "plm", "ivreg", "lmrob")) {
  		return(.summary.object$coefficients[,"Estimate"])
  	}
    if (model.name %in% c("lmer","glmer","nlmer")) {
      coefs <- object.name@fixef
      return(coefs)
    }
    if (model.name %in% c("ergm")) {
      return(.summary.object$coefs[,1])
    }
  	if (model.name %in% c("clm")) {
  	  if (.format.ordered.intercepts == FALSE) {
  	    return(.summary.object$coefficients[(length(object.name$alpha)+1):(length(object.name$coefficients)),1])
  	  }
  	  else {
  	    return(.summary.object$coefficients[,1])
  	  }
  	}
  	else if (model.name %in% c("pmg")) {
  	  return(.summary.object$coefficients)
  	}
    else if (model.name %in% c("zeroinfl", "hurdle")) {
      if (.global.zero.component==FALSE) {
        return(.summary.object$coefficients$count[,"Estimate"])  
      }
      else {
        return(.summary.object$coefficients$zero[,"Estimate"])
      }
    }
  	else if (model.name %in% c("normal.gee", "logit.gee", "probit.gee", "poisson.gee", "gamma.gee", "gee()")) {
  		return(.summary.object$coefficients[,"Estimate"])
  	}
  	else if (model.name %in% c("normal.gam", "logit.gam", "probit.gam", "poisson.gam", "gam()")) {
  		return(.summary.object$p.coeff)
  	}
  	else if (model.name %in% c("coxph", "clogit")) {
  		return(.summary.object$coef[,"coef"])
  	}
  	else if (model.name %in% c("exp","lognorm","weibull","tobit","survreg()")) {
  		return(.summary.object$table[,"Value"])
  	}
    else if (model.name %in% c("rlm")) {
      return(suppressMessages(.summary.object$coefficients[,"Value"]))
    }
  	else if (model.name %in% c("ologit", "oprobit", "polr()")) {
  		coef.temp <- suppressMessages(.summary.object$coefficients[,"Value"])
  		if (.format.ordered.intercepts == FALSE) { return(coef.temp[seq(from=1, to=length(coef.temp)-(length(suppressMessages(.summary.object$lev))-1))]) }
  		else { return(coef.temp) }
  	}
  	else if (model.name %in% c("arima")) {
  		return( object.name$coef )
  	}
  	else if (model.name %in% c("tobit(AER)")){
  	  return(.summary.object$coefficients[,"Estimate"])
  	}
    else if (model.name %in% c("multinom")){
      if (is.null(nrow(.summary.object$coefficients))) {
        coef.temp <- .summary.object$coefficients
      }
      else {
        coef.temp <- .summary.object$coefficients[1,]
      }
      return(coef.temp)
    }
  	else if (model.name %in% c("betareg")){
  	  return(.summary.object$coefficients$mean[,"Estimate"])
  	}
    else if (model.name %in% c("gls")) {
      coef.temp <- object.name$coefficients
      return(coef.temp)
    }
  	else { return(NULL) }

  }
  
  .get.coefficients <-
  function(object.name, user.given=NULL) {
    out <- .get.coefficients.1(object.name, user.given)
                               
    coef.vars <- .coefficient.variables(object.name)
    if (is.null(names(out))) {
      if (length(out) < length(coef.vars)) {
        out.temp <- rep(NA, times=length(coef.vars)-length(out))
        out <- c(out, out.temp)
      }
      else if (length(out) > length(coef.vars)) {
        out <- out[1:length(coef.vars)]
      }
      names(out) <- coef.vars
    }
    return(out)
  }
  
  .is.list.numeric <- 
  function(x) {
    # tolerate NA or NULL
    if (is.null(x)) { return(TRUE) }
    if (!is.list(x)) { return(FALSE) }
    for (i in 1:length(x)) {
      elem <- x[[i]]
      if (!is.null(elem)) {
        if (length(elem) != length(elem[is.numeric(elem) | (is.na(elem))])) { return(FALSE) }
      }
    }
    return(TRUE)
  }
  


###########################################

    ## invisible output
    invisible.output <- NULL
    
    ## error handling
    error.present <- NULL
    
    # get objects
    objects <- as.list(list(...))
    how.many.objects <- length(objects)
    
    ## check if argument input is ok
    if (how.many.objects < 1) { error.present <- c(error.present, "% Error: At least one object is required.\n") }
  
    # identify objects
    for (i in seq(1:how.many.objects)) {
      if (!is.data.frame(objects[[i]])) {
        
        # if zelig$result relevant, identify this automatically
        if (.hasSlot(objects[[i]],"Zt")) {  # use this to eliminate lmer, glmer, nlmer
          if (.model.identify(objects[[i]])=="unknown") { error.present <- c(error.present, "% Error: Unrecognized object type.\n") }
        }
        else {
          if (!is.null(objects[[i]]$zelig.call)) {
            if (!is.null(objects[[i]]$formula)) { formula <- objects[[i]]$formula }
            objects[[i]] <- objects[[i]]$result          
            if (!is.null(formula)) { objects[[i]]$formula2 <- formula }
          }
        
          ###
          if (is.atomic(objects[[i]])) { error.present <- c(error.present, "% Error: Unrecognized object type.\n") }
          else if (.model.identify(objects[[i]])=="unknown") { error.present <- c(error.present, "% Error: Unrecognized object type.\n") }
          else if (.model.identify(objects[[i]])=="unsupported zelig") { error.present <- c(error.present, "% Error: Unsupported 'zelig' model.\n") }
        }  
      }
    }
    
    
    if (!is.character(title)) { error.present <- c(error.present, "% Error: Argument 'title' must be of type 'character.'\n") }
    
    if (!is.character(style)) { error.present <- c(error.present, "% Error: Argument 'style' must be of type 'character.'\n") }
    if (length(style) != 1) { error.present <- c(error.present, "% Error: Argument 'style' must be of length 1.'\n") }
    if (is.character(style)) {
      if (!(tolower(style) %in% c("all","all2","default","commadefault","aer","ajps","ajs","asq","asr","apsr","demography","io","jpam","qje"))) {
        error.present <- c(error.present, "% Error: 'style' not recognized'\n")
      }
    }
    
    if (!is.logical(summary)) { error.present <- c(error.present, "% Error: Argument 'summary' must be of type 'logical' (TRUE/FALSE) \n") }
    if ((length(summary) != 1) & (!is.null(summary))) { error.present <- c(error.present, "% Error: Argument 'summary' must be of length 1.'\n") }
  
    if ((!is.character(covariate.labels)) & (!is.null(covariate.labels))) { error.present <- c(error.present, "% Error: Argument 'covariate.labels' must be NULL (default), or a vector of type 'character.'") }
    if ((!is.character(dep.var.labels)) & (!is.null(dep.var.labels))) { error.present <- c(error.present, "% Error: Argument 'dep.var.labels' must be NULL (default), or a vector of type 'character.'") }
    
    if ((!.is.list.numeric(coef)))  { error.present <- c(error.present, "% Error: Argument 'coef' must be NULL (default), or a list of numerical vectors.") }
    if ((!.is.list.numeric(se)))  { error.present <- c(error.present, "% Error: Argument 'se' must be NULL (default), or a list of numerical vectors.") }
    if ((!.is.list.numeric(t)))  { error.present <- c(error.present, "% Error: Argument 't' must be NULL (default), or a list of numerical vectors.") }
    if ((!.is.list.numeric(p)))  { error.present <- c(error.present, "% Error: Argument 'p' must be NULL (default), or a list of numerical vectors.") }
  
    if (!is.logical(t.auto)) { error.present <- c(error.present, "% Error: Argument 't.auto' must be of type 'logical' (TRUE/FALSE) \n") }
    if (length(t.auto) != 1) { error.present <- c(error.present, "% Error: Argument 't.auto' must be of length 1.'\n") }

    if (!is.logical(p.auto)) { error.present <- c(error.present, "% Error: Argument 't.auto' must be of type 'logical' (TRUE/FALSE) \n") }
    if (length(p.auto) != 1) { error.present <- c(error.present, "% Error: Argument 't.auto' must be of length 1.'\n") }

    if (!is.logical(align)) { error.present <- c(error.present, "% Error: Argument 'align' must be of type 'logical' (TRUE/FALSE) \n") }
    if (length(align) != 1) { error.present <- c(error.present, "% Error: Argument 'align' must be of length 1.'\n") }

    if ((!is.character(decimal.mark)) & (!is.null(decimal.mark))) { error.present <- c(error.present, "% Error: Argument 'decimal.mark' must be NULL (default), or of type 'character.'\n") }
    if ((length(decimal.mark) != 1) & (!is.null(decimal.mark))) { error.present <- c(error.present, "% Error: Argument 'decimal.mark' must be of length 1.'\n") }
  
    if (!is.logical(align)) { error.present <- c(error.present, "% Error: Argument 'align' must be of type 'logical' (TRUE/FALSE) \n") }
    if ((length(align) != 1) & (!is.null(align))) { error.present <- c(error.present, "% Error: Argument 'align' must be of length 1.'\n") }
  
    if ((!is.numeric(digit.separate)) & (!is.null(digit.separate))) { error.present <- c(error.present, "% Error: Argument 'digit.separate' must be NULL (default), a vector of type 'numeric,' or \"lakh\"/\"china\"/\"japan\".\n") }
    
    if ((!is.character(digit.separator)) & (!is.null(digit.separator))) { error.present <- c(error.present, "% Error: Argument 'digit.separator' must be NULL (default), or of type 'character.'\n") }
    if ((length(digit.separator) != 1) & (!is.null(digit.separator))) { error.present <- c(error.present, "% Error: Argument 'digit.separator' must be of length 1.'\n") }
    
    if ((!is.numeric(digits)) & (!is.null(digits))) { 
      if (!is.na(digits)) { error.present <- c(error.present, "% Error: Argument 'digits' must be NULL (default), or of type 'numeric.'\n") }
    }
    if ((length(digits) != 1) & (!is.null(digits))) { 
      if (!is.na(digits)) { error.present <- c(error.present, "% Error: Argument 'digits' must be of length 1.'\n") }
    }
    if (!is.null(digits)) {
      if (!is.na(digits)) {
        if ((digits<0) & (is.numeric(digits))) { error.present <- c(error.present, "% Error: Argument 'digits' must be >= 0.'\n") }
      }
    }
    
    if ((!is.numeric(digits.extra)) & (!is.null(digits.extra))) { error.present <- c(error.present, "% Error: Argument 'digits.extra' must be NULL (default), or of type 'numeric.'\n") }
    if ((length(digits.extra) != 1) & (!is.null(digits.extra))) { error.present <- c(error.present, "% Error: Argument 'digits.extra' must be of length 1.'\n") }
    if (!is.null(digits.extra)) {
      if ((digits.extra<0) & (is.numeric(digits.extra))) { error.present <- c(error.present, "% Error: Argument 'digits.extra' must be >= 0.'\n") }
    }
    
    if ((!is.logical(initial.zero)) & (!is.null(initial.zero))) { error.present <- c(error.present, "% Error: Argument 'initial.zero' must be NULL (default), or of type 'logical' (TRUE/FALSE) \n") }
    if ((length(initial.zero) != 1) & (!is.null(initial.zero))) { error.present <- c(error.present, "% Error: Argument 'initial.zero' must be of length 1.'\n") }
    
    if ((!is.logical(intercept.top)) & (!is.null(intercept.top))) { error.present <- c(error.present, "% Error: Argument 'intercept.top' must be NULL (default), or of type 'logical' (TRUE/FALSE) \n") }
    if ((length(intercept.top) != 1) & (!is.null(intercept.top))) { error.present <- c(error.present, "% Error: Argument 'intercept.top' must be of length 1.'\n") }
    
    if (!is.character(label)) { error.present <- c(error.present, "% Error: Argument 'label' must be of type 'character.'\n") }
  
    if ((!is.logical(model.names)) & (!is.null(model.names))) { error.present <- c(error.present, "% Error: Argument 'model.names' must be of type 'logical' (TRUE/FALSE) \n") }
    if ((length(model.names) != 1) & (!is.null(model.names))) { error.present <- c(error.present, "% Error: Argument 'model.names' must be of length 1.'\n") }
    
    if ((!is.logical(model.numbers)) & (!is.null(model.numbers))) { error.present <- c(error.present, "% Error: Argument 'model.numbers' must be of type 'logical' (TRUE/FALSE) \n") }
    if ((length(model.numbers) != 1) & (!is.null(model.numbers))) { error.present <- c(error.present, "% Error: Argument 'model.numbers' must be of length 1.'\n") }
    
    if ((!is.character(notes)) & (!is.null(notes))) { error.present <- c(error.present, "% Error: Argument 'notes' must be NULL (default), or a vector of type 'character.'") }
    
    if (!is.null(notes.align)) {
      if (!(tolower(notes.align) %in% c("l","c","r"))) { error.present <- c(error.present, "% Error: Argument 'notes.align' must be NULL (default), or \"l\"/\"c\"/\"r\".") }  
    }
    if ((length(notes.align) != 1) & (!is.null(notes.align))) { error.present <- c(error.present, "% Error: Argument 'notes.align' must be of length 1.'\n") }
    
    if ((!is.character(notes.label)) & (!is.null(notes.label))) { error.present <- c(error.present, "% Error: Argument 'notes.label' must be NULL (default), or of type 'character.'") }
    if ((length(notes.label) != 1) & (!is.null(notes.label))) { error.present <- c(error.present, "% Error: Argument 'notes.label' must be of length 1.'\n") }
    
    if ((!is.character(omit)) & (!is.null(omit))) { error.present <- c(error.present, "% Error: Argument 'omit' must be NULL (default; no omissions), or a vector of regular expressions of type 'character.'") }
    if ((!is.character(omit.labels)) & (!is.null(omit.labels))) { error.present <- c(error.present, "% Error: Argument 'omit' must be NULL (default; no omissions), or a vector of type 'character.'") }
    if (!is.null(omit.labels)) {
      if (length(omit) != length(omit.labels)) { error.present <- c(error.present, "% Error: Arguments 'omit.labels' must be NULL (default; no omissions), or equal in length to 'omit.labels'.'\n") }
    }
  
    if ((!is.character(omit.stat)) & (!is.null(omit.stat))) { error.present <- c(error.present, "% Error: Argument 'omit.stat' must be NULL (default), or a vector of type 'character.'") }
    omit.stat.acceptable <- c("all","n","rsq","adj.rsq","max.rsq","ll","aic","bic","scale","ubre","sigma2","ser","f","theta","chi2","wald","lr","logrank","null.dev","res.dev")     # list of statistic codes that are acceptable
    if (is.character(omit.stat)) {
      is.acceptable <- unique(tolower(omit.stat) %in% omit.stat.acceptable)
      if (length(is.acceptable)>1) { is.acceptable <- FALSE }
      if (!is.acceptable) { error.present <- c(error.present, "% Error: Unknown statistic in 'omit.stat' argument.") }
    } 
  
    if ((!is.character(omit.yes.no)) & (!is.null(omit.yes.no))) { error.present <- c(error.present, "% Error: Argument 'omit.yes.no' must be a vector of type 'character.'") }
    if ((length(omit.yes.no) != 2) & (!is.null(omit.yes.no))) { error.present <- c(error.present, "% Error: Argument 'omit.yes.no' must be of length 2.'\n") }
    
    if (!is.logical(ord.intercepts)) { error.present <- c(error.present, "% Error: Argument 'ord.intercepts' must be of type 'logical' (TRUE/FALSE) \n") }
    if (length(ord.intercepts) != 1) { error.present <- c(error.present, "% Error: Argument 'ord.intercepts' must be of length 1.'\n") }
    
    if ((!is.character(star.char)) & (!is.null(star.char))) { error.present <- c(error.present, "% Error: Argument 'star.char' must be NULL (default), or of type 'character.'\n") }
    if ((length(star.char) != 1) & (!is.null(star.char))) { error.present <- c(error.present, "% Error: Argument 'star.char' must be of length 1.'\n") }
    
    if (!is.null(star.cutoffs)) {
      if (sum(is.na(star.cutoffs)) != length(star.cutoffs)) {
        if (!is.numeric(star.cutoffs)) { error.present <- c(error.present, "% Error: Argument 'star.cutoffs' must be NULL (default), or a vector of type 'numeric.'\n") }
      }
      if ( !((length(star.cutoffs) >= 1) & (length(star.cutoffs) <= 3)) & (!is.null(star.cutoffs))) { error.present <- c(error.present, "% Error: Argument 'star.cutoffs' must be a vector with at least one and at most three elements.\n") }
      if (sum(star.cutoffs[!is.na(star.cutoffs)] == sort(star.cutoffs, decreasing = TRUE, na.last=NA)) != length(star.cutoffs[!is.na(star.cutoffs)])) { error.present <- c(error.present, "% Error: The elements of 'star.cutoffs' must be in weakly decreasing order.\n") }
    }
  
    if (!is.logical(zero.component)) { error.present <- c(error.present, "% Error: Argument 'zero.component' must be of type 'logical' (TRUE/FALSE) \n") }
    if ((length(zero.component) != 1) & (!is.null(zero.component))) { error.present <- c(error.present, "% Error: Argument 'zero.component' must be of length 1.'\n") }
    
    if (!is.logical(summary.logical)) { error.present <- c(error.present, "% Error: Argument 'summary.logical' must be of type 'logical' (TRUE/FALSE) \n") }
    if (length(summary.logical) != 1) { error.present <- c(error.present, "% Error: Argument 'summary.logical' must be of length 1.'\n") }

    if (!is.logical(nobs)) { error.present <- c(error.present, "% Error: Argument 'nobs' must be of type 'logical' (TRUE/FALSE) \n") }
    if (length(nobs) != 1) { error.present <- c(error.present, "% Error: Argument 'nobs' must be of length 1.'\n") }
    
    if (!is.logical(mean.sd)) { error.present <- c(error.present, "% Error: Argument 'mean.sd' must be of type 'logical' (TRUE/FALSE) \n") }
    if (length(mean.sd) != 1) { error.present <- c(error.present, "% Error: Argument 'mean.sd' must be of length 1.'\n") }
    
    if (!is.logical(min.max)) { error.present <- c(error.present, "% Error: Argument 'min.max' must be of type 'logical' (TRUE/FALSE) \n") }
    if (length(min.max) != 1) { error.present <- c(error.present, "% Error: Argument 'min.max' must be of length 1.'\n") }
    
    if (!is.logical(median)) { error.present <- c(error.present, "% Error: Argument 'median' must be of type 'logical' (TRUE/FALSE) \n") }
    if (length(median) != 1) { error.present <- c(error.present, "% Error: Argument 'median' must be of length 1.'\n") }
    
    if (!is.logical(iqr)) { error.present <- c(error.present, "% Error: Argument 'iqr' must be of type 'logical' (TRUE/FALSE) \n") }
    if (length(iqr) != 1) { error.present <- c(error.present, "% Error: Argument 'iqr' must be of length 1.'\n") }
    
    ## decide what style to use here: start with all settings, and then make adjustment based on desired journal

    # initialize pseudo-global variables at NULL
    .summary.object <- NULL
    .global.dependent.variables.written <- NULL
    .global.coefficients <- NULL
    .format.model.left <- NULL
    .format.model.right <- NULL
    .which.variable.label <- NULL
    .return.value <- NULL
    .publish.horizontal.line <- NULL
    .table.part.published <- NULL
    .format.omit.table <- NULL
  
    # should we include a summary statistics table when given a data frame
    .global.summary <- TRUE

    # info about the package and author
    .global.package.name <- "StarGazer"
    .global.package.version <- "3.0"
    .global.package.author.name <- "Marek Hlavac"
    .global.package.author.affiliation <- "Harvard University"
    .global.package.author.email <- "hlavac at fas.harvard.edu"
    
    # statistics (.global variables)
    .global.formulas.rhs <- NULL
    .global.models <- NULL
    .global.dependent.variables <- NULL
    .global.coefficient.variables <- NULL
    .global.std.errors <- NULL
    .global.t.stats <- NULL
    .global.p.values <- NULL
    .global.N <- NULL
    .global.LL <- NULL
    .global.R2 <- NULL
    .global.max.R2 <- NULL # maximum possible R2
    .global.adj.R2 <- NULL
    .global.AIC <- NULL
    .global.BIC <- NULL
    .global.scale <- NULL   # estimated scale parameter (gee)
    .global.UBRE <- NULL    # UBRE score (GAM)
    .global.sigma2 <- NULL  # sigma2 from arima
    .global.theta <- NULL   # theta from negative binomial
    .global.zero.component <- NULL # zero, as opposed to count, component in hurdle and zeroinfl
    
    # with degrees of freedom
    .global.SER <- NULL   # residual standard error; standard error of the regression
    .global.F.stat <- NULL # F-statistic for the regression
    .global.chi.stat <- NULL  # chi-squared statistic
    .global.wald.stat <- NULL # Wald test statistic (for coxph)
    .global.lr.stat <- NULL  # LR test statistic (for coxph)
    .global.logrank.stat <- NULL # Score (logrank) test (for coxph)
    .global.null.deviance <- NULL 
    .global.residual.deviance <- NULL
    
    # intercept strings
    .global.intercept.strings <- c("(Intercept)", "(intercept)")
    
    # .formatting: Default
    .format.space.size <- "-1.8ex"
    
    .format.dependent.variable.text <- "\\textit{Dependent variable:}"
    .format.dependent.variable.text.underline <- TRUE
    .format.dependent.variable.text.on <- TRUE
    
    .format.dep.var.labels <- NULL
    .format.covariate.labels <- NULL
    
    .format.dependent.variables.text <- ""
    .format.underline.dependent.variables <- TRUE
    .format.dependent.variables.left <- ""
    .format.dependent.variables.right <- ""
    .format.dependent.variables.capitalize <- FALSE
    
    .format.ordered.intercepts <- TRUE
    
    # model numbers
    .format.model.numbers <- TRUE
    
    # names for models
    .format.model.names.include <- TRUE
    .format.model.names <- NULL
    .format.model.names <- cbind(c("aov","ANOVA",""), c("arima","ARIMA",""), c("blogit","bivariate","logistic"))
    .format.model.names <- cbind(.format.model.names, c("bprobit","bivariate","probit"), c("betareg", "beta",""), c("chopit","compound hierarchical","ordered probit"))
    .format.model.names <- cbind(.format.model.names, c("clm","cumulative","link"), c("cloglog.net","network compl.","log log"), c("clogit","conditional","logistic"), c("coxph","Cox","prop. hazards"))
    .format.model.names <- cbind(.format.model.names, c("ei.dynamic","Quinn dynamic","ecological inference"), c("ei.hier","$2 \times 2$ hierarchical","ecological inference"))
    .format.model.names <- cbind(.format.model.names, c("ei.RxC","hierarchical multinominal-Dirichlet","ecological inference"), c("exp","exponential",""), c("ergm","exponential family","random graph"))
    .format.model.names <- cbind(.format.model.names, c("factor.bayes","Bayesian","factor analysis"), c("factor.mix","mixed data","factor analysis"))
    .format.model.names <- cbind(.format.model.names, c("factor.ord","ordinal data","factor analysis"), c("gamma","gamma",""))
    .format.model.names <- cbind(.format.model.names, c("gamma.gee","gamma generalized","estimating equation"), c("gamma.mixed","mixed effects","gamma"))
    .format.model.names <- cbind(.format.model.names, c("gamma.net","network","gamma"), c("gamma.survey","survey-weighted","gamma"), c("gls","generalized","least squares"))
    .format.model.names <- cbind(.format.model.names, c("irt1d","IRT","(1-dim.)"), c("irtkd","IRT","(k-dim.)"))
    .format.model.names <- cbind(.format.model.names, c("logit","logistic",""), c("logit.bayes","Bayesian","logistic"))
    .format.model.names <- cbind(.format.model.names, c("logit.gam","GAM","(logistic)"), c("logit.gee","logistic generalized","estimating equation"))
    .format.model.names <- cbind(.format.model.names, c("logit.mixed","mixed effects","logistic"), c("logit.net","network","logistic"))
    .format.model.names <- cbind(.format.model.names, c("logit.survey","survey-weighted","logistic"), c("lognorm","log-normal",""))
    .format.model.names <- cbind(.format.model.names, c("lmer","linear","mixed-effects"), c("glmer","generalized linear","mixed-effects"), c("nlmer","non-linear","mixed-effects"))
    .format.model.names <- cbind(.format.model.names, c("ls","OLS",""), c("ls.mixed","mixed effect","linear"), c("lmrob","MM-type","linear"))
    .format.model.names <- cbind(.format.model.names, c("ls.net","network","least squares"), c("mlogit","multinomial","logistic"))
    .format.model.names <- cbind(.format.model.names, c("mlogit.bayes","Bayesian","multinomial logistic"), c("negbin","negative","binomial"), c("normal","normal",""))
    .format.model.names <- cbind(.format.model.names, c("multinom","multinomial log-linear","(neural networks)"))
    .format.model.names <- cbind(.format.model.names, c("normal.bayes","Bayesian","normal"), c("normal.gam","GAM","(continuous)"))
    .format.model.names <- cbind(.format.model.names, c("normal.gee","normal generalized","estimating equation"), c("normal.net","network","normal"))
    .format.model.names <- cbind(.format.model.names, c("normal.survey","survey-weighted","normal"), c("ologit","ordered","logistic"))
    .format.model.names <- cbind(.format.model.names, c("oprobit","ordered","probit"), c("oprobit.bayes","Bayesian","ordered probit"))
    .format.model.names <- cbind(.format.model.names, c("pmg","mean","groups"), c("poisson","Poisson",""), c("poisson.bayes","Bayesian","Poisson"))
    .format.model.names <- cbind(.format.model.names, c("poisson.gam","GAM","(count)"), c("poisson.mixed","mixed effects","Poisson"))
    .format.model.names <- cbind(.format.model.names, c("poisson.survey","survey-weighted","Poisson"), c("poisson.gee","Poisson generalized","estimation equation"))
    .format.model.names <- cbind(.format.model.names, c("probit","probit",""), c("probit.bayes","Bayesian","probit"))
    .format.model.names <- cbind(.format.model.names, c("probit.gam","GAM","(probit)"), c("probit.gee","probit generalized","estimating equation"))
    .format.model.names <- cbind(.format.model.names, c("probit.mixed","mixed effects","probit"), c("probit.net","network","probit"))
    .format.model.names <- cbind(.format.model.names, c("probit.survey","survey-weighted","probit"), c("relogit","rare events","logistic"))
    .format.model.names <- cbind(.format.model.names, c("rlm","robust","linear"), c("sur","SUR",""), c("threesls","3SLS",""))
    .format.model.names <- cbind(.format.model.names, c("tobit","Tobit",""), c("tobit(AER)","Tobit",""), c("tobit.bayes","Bayesian","Tobit"))
    .format.model.names <- cbind(.format.model.names, c("twosls","2SLS",""), c("weibull","Weibull",""))
    .format.model.names <- cbind(.format.model.names, c("zeroinfl","zero-inflated","count data"), c("hurdle","hurdle",""))
    .format.model.names <- cbind(.format.model.names, c("plm","panel","linear"), c("ivreg","instrumental","variable"))
  
    # if you use, say, glm() that does not correspond to one of the pre-defined models, put this as family and link
    .format.model.function <- TRUE
    .format.model.family <- ""
    .format.model.dist <- ""
    .format.model.link <- "link = "
    
    ## names for journal/output styles
    # economics
    .journal.style.names <- cbind(c("aer","American Economic Review"), c("qje","Quarterly Journal of Economics"), c("econometrica","Econometrica"))
    .journal.style.names <- cbind(.journal.style.names, c("jpe","Journal of Political Economy"), c("jel","Journal of Economic Literature"))
    .journal.style.names <- cbind(.journal.style.names, c("jep","Journal of Economic Perspestives"))
    
    .format.coefficient.variables.capitalize <- FALSE
    .format.coefficient.variables.left <- ""
    .format.coefficient.variables.right <- ""
    .format.coefficient.table.parts <- c("variable name","coefficient*","standard error"," ")
    
    ## .formatting of numerical output
    # keep initial zeros?
    .format.initial.zero <- TRUE
    # if all zeros, keep going until you find a non-zero digit
    .format.until.nonzero.digit <- TRUE
    .format.max.extra.digits <- 2
    
    ## threshholds for the stars
    .format.significance.star <- "*"
    .format.three.stars <- 0.01
    .format.two.stars <- 0.05
    .format.one.star <- 0.1
    
    .format.std.errors.left <- "("
    .format.std.errors.right <- ")"
    
    .format.p.values.left <- "p = "
    .format.p.values.right <- ""
    
    .format.t.stats.left <- "t = "
    .format.t.stats.right <- ""
    
    .format.models.text <- ""
    .format.models.left <- "\\textit{"
    .format.models.right <- "}"
    .format.underline.models <- FALSE
    .format.models.skip.if.one <- TRUE # skip models section if only one model in table?
    
    .format.numbers.text <- ""
    .format.numbers.left <- "("
    .format.numbers.right <- ")"
    .format.numbers.roman <- FALSE
    
    .format.digit.separator.where <- c(3)    # how 'often' to separate digits (e.g., thousands separator = 3)
    .format.digit.separator <- ","
    .format.round.digits <- 3
    # for decimal comma use: .format.decimal.character <- "{,}"
    .format.decimal.character <- "."
    .format.dec.mark.align <- FALSE
    
    .format.table.parts <- c("=!","dependent variable label","dependent variables","models","numbers","-","coefficients","-","omit","-","additional","N","R-squared","adjusted R-squared","max R-squared","log likelihood","scale","sigma2","theta(se)*", "AIC","BIC","UBRE", "SER(df)","F statistic(df)*(p)","chi2(df)*(p)","Wald(df)*(p)","LR(df)*(p)","logrank(df)*(p)","null deviance(df)","residual deviance(df)","=!","notes")
    
    .format.omit.regexp <- NULL
    .format.omit.labels <- NULL
    .format.omit.yes <- "Yes"
    .format.omit.no <- "No"
    
    .format.N <- "Observations"
    .format.LL <- "Log likelihood"
    .format.R2 <- "R$^{2}$"
    .format.max.R2 <- "Max. Possible R$^{2}$"
    .format.adj.R2 <- "Adjusted R$^{2}$"
    .format.scale <- "Scale Parameter"
    .format.UBRE <- "UBRE"
    .format.AIC <- "Akaike Inf. Crit."
    .format.BIC <- "Bayesian Inf. Crit."
    .format.sigma2 <- "$\\sigma^{2}$"
    .format.theta <- "$\\theta$"
    
    .format.SER <- "Residual Std. Error"
    .format.F.stat <- "F statistic"
    .format.chi.stat <- "$\\chi^{2}$"
    .format.wald.stat <- "Wald Test"
    .format.lr.stat <- "LR Test"
    .format.logrank.stat <- "Score (Logrank) Test"
    .format.null.deviance <- "Null Deviance"
    .format.residual.deviance <- "Residual Deviance"
    
    .format.df.left <- "(df = "
    .format.df.right <- ")"
    .format.df.separator <- "; "
    .format.intelligent.df <- TRUE
    
    # this is for se, tstat, p.values at the bottom of the table, by statistics
    .format.se.left <- " ("
    .format.se.right <- ")"
    .format.tstat.left <- " (z = "
    .format.tstat.right <- ")"
    .format.p.value.left <- " (p = "
    .format.p.value.right <- ")"
    
    .format.intercept.name <- "Constant"
    .format.intercept.bottom <- TRUE
    .format.note <- "\\textit{Note:} "
    .format.note.alignment <- "r"
    .format.note.content <- c("$^{*}$p$<$[0.*]; $^{**}$p$<$[0.**]; $^{***}$p$<$[0.***]")
    
    # additional lines
    .format.additional.lines.titles <- NULL
    .format.additional.lines.values <- NULL
    
    
    #### summary statistic table
    .format.s.statistics.names <- cbind(c("N","N"), c("nmiss","missing"), c("mean","Mean"), c("sd","St. Dev."), c("median","Median"), c("min","Min"), c("max","Max"), c("mad","Median Abs. Dev."), c("p","Pctl(!)"))
    .format.s.stat.parts <- c("=!","stat names","-!","statistics1","-!")
    .format.s.statistics.list <- c("N","mean","sd","min","p25","median","p75","max")
    
    .format.s.statistics.names.left <- ""
    .format.s.statistics.names.right <- ""
    .format.s.statistics.names.label <- "Statistic"
    
    .format.s.coefficient.variables.capitalize <- FALSE
    .format.s.coefficient.variables.left <- ""
    .format.s.coefficient.variables.right <- ""
    
    .format.s.round.digits <- 3
    
    .format.s.note <- ""
    .format.s.note.alignment <- "l"
    .format.s.note.content <- c("")

    ####
    .adjust.settings.style(style)
    
    # continue only if no errors
    if (is.null(error.present)) {
      
      # summary statistic table or regular table of data frame contents
      if (!is.null(.global.summary)) { .global.summary <- summary }
      
      ## use formatting arguments
      
      # if not case-sensitive, transfer to lower case
      if (!is.null(digit.separate)) { digit.separate <- tolower(digit.separate) }
      
      # column, dependent variable and covariate labels
      .format.covariate.labels <- covariate.labels
      .format.dep.var.labels <- dep.var.labels
      
      # omit
      .format.omit.regexp <- omit
      .format.omit.labels <- omit.labels
      if (!is.null(omit.yes.no)) { 
        .format.omit.yes <- omit.yes.no[1]
        .format.omit.no <- omit.yes.no[2]
      }
      
      # remove omitted statistics from table parts
      if (!is.null(omit.stat)) {
        .lower.omit.stat <- tolower(omit.stat)    # make it all lower-case
        if ("all" %in% .lower.omit.stat) { .lower.omit.stat <- omit.stat.acceptable }
        if ("n" %in% .lower.omit.stat) { .format.table.parts <- .format.table.parts[.format.table.parts!="N"] }
        if ("rsq" %in% .lower.omit.stat) { .format.table.parts <- .format.table.parts[.format.table.parts!="R-squared"] }
        if ("adj.rsq" %in% .lower.omit.stat) { .format.table.parts <- .format.table.parts[.format.table.parts!="adjusted R-squared"] }
        if ("max.rsq" %in% .lower.omit.stat) { .format.table.parts <- .format.table.parts[.format.table.parts!="max R-squared"] }
        if ("ll" %in% .lower.omit.stat) { .format.table.parts <- .format.table.parts[.format.table.parts!="log likelihood"] }
        if ("scale" %in% .lower.omit.stat) { .format.table.parts <- .format.table.parts[.format.table.parts!="scale"] }
        if ("sigma2" %in% .lower.omit.stat) { .format.table.parts <- .format.table.parts[.format.table.parts!="sigma2"] }        
        if ("theta" %in% .lower.omit.stat) { .format.table.parts <- .format.table.parts[substr(.format.table.parts,1,5)!="theta"] }
        if ("aic" %in% .lower.omit.stat) { .format.table.parts <- .format.table.parts[.format.table.parts!="AIC"] }
        if ("bic" %in% .lower.omit.stat) { .format.table.parts <- .format.table.parts[.format.table.parts!="BIC"] }
        if ("ubre" %in% .lower.omit.stat) { .format.table.parts <- .format.table.parts[.format.table.parts!="UBRE"] }
        if ("ser" %in% .lower.omit.stat) { .format.table.parts <- .format.table.parts[substr(.format.table.parts,1,3)!="SER"] }
        if ("f" %in% .lower.omit.stat) { .format.table.parts <- .format.table.parts[substr(.format.table.parts,1,11)!="F statistic"] }
        if ("chi2" %in% .lower.omit.stat) { .format.table.parts <- .format.table.parts[substr(.format.table.parts,1,4)!="chi2"] }
        if ("wald" %in% .lower.omit.stat) { .format.table.parts <- .format.table.parts[substr(.format.table.parts,1,4)!="Wald"] }
        if ("lr" %in% .lower.omit.stat) { .format.table.parts <- .format.table.parts[substr(.format.table.parts,1,2)!="LR"] }
        if ("logrank" %in% .lower.omit.stat) { .format.table.parts <- .format.table.parts[substr(.format.table.parts,1,7)!="logrank"] }
        if ("null.dev" %in% .lower.omit.stat) { .format.table.parts <- .format.table.parts[substr(.format.table.parts,1,13)!="null deviance"] }
        if ("res.dev" %in% .lower.omit.stat) { .format.table.parts <- .format.table.parts[substr(.format.table.parts,1,17)!="residual deviance"] }
      }
      
      # digits, initial.zeros, decimal characters
      if (!is.null(decimal.mark)) { .format.decimal.character <- decimal.mark }
      if (!is.null(align)) { .format.dec.mark.align <- align }
      if (!is.null(digit.separator)) { .format.digit.separator <- digit.separator }
      if (!is.null(initial.zero)) { .format.initial.zero <- initial.zero }
      
      if (!is.null(digit.separate)) { 
        if (digit.separate=="lakh") { .format.digit.separator.where <- c(3,2) }  # lakhs 
        else if ((digit.separate=="china") | (digit.separate=="japan")) { .format.digit.separator.where <- 4 }
        else { .format.digit.separator.where <- digit.separate}
      }
      
      if (!is.null(digits)) { 
        .format.round.digits <- digits 
        .format.s.round.digits <- digits
      }
      
      if (!is.null(digits.extra)) { 
        .format.max.extra.digits <- digits.extra
        if (digits.extra>=1) { .format.until.nonzero.digit <- TRUE }
        else ( .format.until.nonzero.digit <- FALSE )
      }
      
      # intercept top
      if (!is.null(intercept.top)) { .format.intercept.bottom <- !(intercept.top) }    
      
      # intercept top
      if (!is.null(model.names)) { 
        .format.model.names.include <- model.names 
        if (model.names == TRUE) { .format.models.skip.if.one <- FALSE }
      }    
      if (!is.null(model.numbers)) { .format.model.numbers <- model.numbers }    
      
      # significance stars
      if (!is.null(star.char)) { .format.significance.star <- star.char }
      if (!is.null(star.cutoffs)) { 
        # assign cutoff values
        
        .format.one.star <- star.cutoffs[1]
        .format.two.stars <- star.cutoffs[2]
        .format.three.stars <- star.cutoffs[3]
      }
      
      # zero vs. count component
      .global.zero.component <- zero.component
      
      # notes
      
      # replace star cutoffs in the notes section
      if (!is.na(.format.one.star)) {
        .format.note.content <- gsub("[.*]", gsub("^[0]+", "",.format.one.star), .format.note.content, fixed=TRUE)  
        .format.note.content <- gsub("[0.*]", .format.one.star, .format.note.content, fixed=TRUE)
        .format.note.content <- gsub("[*]", .format.one.star*100, .format.note.content, fixed=TRUE)        
      }
      if (!is.na(.format.two.stars)) {
        .format.note.content <- gsub("[.**]", gsub("^[0]+", "",.format.two.stars), .format.note.content, fixed=TRUE)        
        .format.note.content <- gsub("[0.**]", .format.two.stars, .format.note.content, fixed=TRUE)
        .format.note.content <- gsub("[**]", .format.two.stars*100, .format.note.content, fixed=TRUE)        
      }
      if (!is.na(.format.three.stars)) {
        .format.note.content <- gsub("[.***]", gsub("^[0]+", "",.format.three.stars), .format.note.content, fixed=TRUE)
        .format.note.content <- gsub("[0.***]", .format.three.stars, .format.note.content, fixed=TRUE)
        .format.note.content <- gsub("[***]", .format.three.stars*100, .format.note.content, fixed=TRUE)
      }
      
      if (!is.null(notes)) { 
        .format.note.content <- notes
        .format.s.note.content <- notes
      }
      if (!is.null(notes.align)) { 
        .format.note.alignment <- notes.align 
        .format.s.note.alignment <- notes.align
      }
      if (!is.null(notes.label)) { 
        .format.note <- notes.label 
        .format.s.note <- notes.label
      }    
      
      # ordered probit/logit, etc. - report intercepts?
      .format.ordered.intercepts <- ord.intercepts
      
      # report logical variables in summary statistics tables?
      .format.summ.logical <- summary.logical
      
      # summary statistics - what statistics to report
      if (!nobs) { .format.s.statistics.list <- .format.s.statistics.list[.format.s.statistics.list!="N"] }
      if (!mean.sd) { .format.s.statistics.list <- .format.s.statistics.list[(.format.s.statistics.list!="mean")&(.format.s.statistics.list!="sd")]}
      if (!min.max) { .format.s.statistics.list <- .format.s.statistics.list[(.format.s.statistics.list!="min")&(.format.s.statistics.list!="max")]}
      if (!median) { .format.s.statistics.list <- .format.s.statistics.list[.format.s.statistics.list!="median"] }
      if (!iqr) { .format.s.statistics.list <- .format.s.statistics.list[(.format.s.statistics.list!="p25")&(.format.s.statistics.list!="p75")]}
      
      # intelligent division of regression tables vs. summary statistics tables
      regression.table.objects <- NULL
      number.of.table <- 0
      title.table <- NULL
      label.table <- NULL
      for (i in seq(1:how.many.objects)) {
        if (is.data.frame(objects[[i]])==TRUE) {
          if (!is.null(regression.table.objects)) { 
            number.of.table <- number.of.table + 1    # allows for multiple table titles and labels
            
            if (!is.na(title[number.of.table])) { .format.title <- title[number.of.table] }
            else { .format.title <- title[length(title)] }
            
            if (!is.na(label[number.of.table])) { .format.label <- label[number.of.table] }
            else { .format.label <- label[length(label)] }
            
            do.call(.stargazer.reg.table, as.list(objects[regression.table.objects]))
            invisible.output <- c(invisible.output, invisible(capture.output(do.call(.stargazer.reg.table, as.list(objects[regression.table.objects])),file=NULL)) )
          }
          
          number.of.table <- number.of.table + 1
          if (!is.na(title[number.of.table])) { .format.title <- title[number.of.table] }
          else { .format.title <- title[length(title)] }
          
          if (!is.na(label[number.of.table])) { .format.label <- label[number.of.table] }
          else { .format.label <- label[length(label)] }
          
          if (.global.summary==TRUE) {
            .stargazer.summ.stat.table(objects[[i]])
            invisible.output <- c(invisible.output, invisible(capture.output(.stargazer.summ.stat.table(objects[[i]]),file=NULL)) )
          }
          else {
            .stargazer.data.frame.table(objects[[i]])
            invisible.output <- c(invisible.output, invisible(capture.output(.stargazer.data.frame.table(objects[[i]]),file=NULL)) )
          }
          regression.table.objects <- NULL
        }
        else {
          regression.table.objects <- c(regression.table.objects, i)
        }
      }
      
      if (!is.null(regression.table.objects)) {	
        number.of.table <- number.of.table + 1
        if (!is.na(title[number.of.table])) { .format.title <- title[number.of.table] }
        else { .format.title <- title[length(title)] }
        
        if (!is.na(label[number.of.table])) { .format.label <- label[number.of.table] }
        else { .format.label <- label[length(label)] }
        
        do.call(.stargazer.reg.table, as.list(objects[regression.table.objects]))  
        invisible.output <- c(invisible.output, invisible(capture.output(do.call(.stargazer.reg.table, as.list(objects[regression.table.objects])),file=NULL)) )
      }
    }
    else { 
      cat(error.present, sep="")
      invisible.output <- error.present 
    }
    return(invisible(invisible.output))
}

