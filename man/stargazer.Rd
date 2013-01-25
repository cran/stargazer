\name{stargazer}
\alias{stargazer}
\alias{stargazer-package}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{stargazer
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The \code{stargazer} command produces LaTeX code for well-formatted tables that hold regression analysis results from several models side-by-side, as well as summary statistics. It supports model objects from  \code{betareg} (betareg), \code{coxph} (survival), \code{clm} (ordinal), \code{clogit} (survival), \code{ergm} (ergm), \code{gam} (mgcv), \code{gee} (gee), \code{glm} (stats), \code{glmer} (lme4), \code{gls} (nlme), \code{hurdle} (pscl),  \code{ivreg} (AER), \code{lm} (stats), \code{lmer} (lme4), \code{lmrob} (robustbase), \code{multinom} (nnet), \code{nlmer} (lme4), \code{plm} (plm), \code{pmg} (plm), \code{polr} (MASS), \code{rlm} (MASS), \code{svyglm} (survey), \code{survreg} (survival), \code{tobit} (AER), \code{zeroinfl} (pscl), as well as from the implementation of these in \code{zelig}. It also supports the following \code{zelig} models for social network analysis: \code{"cloglog.net"}, \code{"gamma.net"}, \code{"probit.net"} and \code{"logit.net"}.
}
\usage{
stargazer(  ..., 
            title = "", style = "default", summary = TRUE,
            covariate.labels = NULL, dep.var.labels = NULL, align = FALSE,
            coef = NULL, se = NULL, t = NULL, p = NULL,
            t.auto = TRUE, p.auto = TRUE,
            decimal.mark = NULL, digit.separate = NULL, digit.separator = NULL,
            digits = NULL, digits.extra = NULL, initial.zero = NULL, intercept.top = NULL, 
            label = "", model.names = NULL, model.numbers = NULL, 
            notes = NULL, notes.align = NULL, notes.label = NULL, 
            omit = NULL, omit.labels = NULL, omit.stat = NULL, omit.yes.no = c("Yes", "No"), 
            ord.intercepts = FALSE, star.char = NULL, star.cutoffs = NULL, zero.component = FALSE,
            summary.logical = TRUE,
            nobs = TRUE, mean.sd = TRUE, min.max = TRUE, median = FALSE, iqr = FALSE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{one or more model objects (for regression analysis tables) or data frames (for summary statistics, or direct output of content).
%%     ~~Describe \code{\dots} here~~
}
  \item{title}{character string vector with titles for the tables.
%%     ~~Describe \code{title} here~~
}
  \item{style}{character string that specifies what style, typically designed to resemble an existing academic journal, should be used in producing the tables. This argument is not case-sensitive. See \link{list of supported styles}. 
%%     ~~Describe \code{style} here~~
}

  \item{summary}{a logical value indicating whether the package should output a summary statistics table when given a data frame. If \code{FALSE}, the package will instead output the contents of the data frame.}

\item{covariate.labels}{a character string vector of labels for covariates in regression tables. A value of \code{NA} for any element means that \code{stargazer} will print the corresponding variable name. In the default case of \code{NULL}, variable names are printed.%%     ~~Describe \code{covariate.labels} here~~
}

\item{dep.var.labels}{a character string vector of labels for the dependent variables in regression tables. A value of \code{NA} for any element means that \code{stargazer} will print the corresponding variable name. In the default case of \code{NULL}, variable names are printed.%%     ~~Describe \code{dep.var.labels} here~~
}

  \item{align}{a logical value indicating whether numerical values in the same column should be aligned at the decimal mark. Requires \code{\\usepackage{dcolumn}} in LaTeX preamble.}
  
    \item{coef}{a list of numerical vectors that will replace the default coefficient values for each model. Element names will be used to match coefficients to individual covariates, and should therefore match covariate names. A \code{NULL} vector indicates that, for a given model, the default set of coefficients should be used. By contrast, an \code{NA} vector means that all of the model's coefficients should be left blank.}
    
    \item{se}{a list of numerical vectors that will replace the default coefficient values for each model. Behaves exactly like the argument \code{coef}.}
    
    \item{t}{a list of numerical vectors that will replace the default test statistics (e.g., t-scores, or z-scores) for each model. Like \code{coef} and \code{se}, test statistics are matched to covariates by their element names.}
    
    \item{p}{a list of numerical vectors that will replace the default p-values for each model. Matched by element names. These will form the basis of decisions about significance stars.}
    
  \item{t.auto}{a logical value that indicates whether \code{stargazer} should calculate the test statistics (i.e., the z-scores) automatically if standard errors are supplied by the user (from argument \code{se}). If \code{FALSE}, the package will use model's default values if \code{t} is \code{NULL}. }
 
  \item{p.auto}{a logical value that indicates whether \code{stargazer} should calculate the p-values, using the standard normal distribution, standard errors are supplied by the user (from argument \code{se}). If \code{FALSE}, the package will use model's default values if \code{p} is \code{NULL}. }
 
  \item{decimal.mark}{character string that will serve as the decimal mark. For instance, the string \code{","} will represent decimal commas, while \code{"."} means tables will use decimal points.}


  \item{digit.separate}{a numerical vector that indicates where digit separators should be placed. The first element of the vector indicates the number of digits (counted from the decimal mark to the left) that will be separated. The second element indicates the number of digits that will be separated from that 'first' separator, and so on. A value of \code{3} corresponds to a thousands separator, while a value of \code{0} indicates no separation. Alternatively, \code{digit.separate} can be one of the following character strings: \code{"lakh"} (equivalent to \code{c(4,3)}), \code{"china"} or \code{"japan"} (both equivalent to a value of \code{4}).
%%     ~~Describe \code{digit.separate} here~~
}
  \item{digit.separator}{character string that will serve as the digit (e.g., thousands) separator. Commonly used strings include \code{","} for a comma separator, \code{" "} for a single space separator, and \code{""} for no separation.
%%     ~~Describe \code{digit.separator} here~~
}
  \item{digits}{integer that indicates how many decimal places should be used. A value of \code{NA} indicates that no rounding should be done at all, and that all available decimal places should be reported.
%%     ~~Describe \code{digits} here~~
}
  \item{digits.extra}{integer indicating the maximum number of additional decimal places to be used if a number, rounded to \code{digits} decimal places, is equal to zero.
%%     ~~Describe \code{digits.extra} here~~
}
  \item{initial.zero}{a logical value indicating whether an initial zero should be printed before the decimal mark if a number is between 0 and 1. 
%%     ~~Describe \code{initial.zero} here~~
}
  \item{intercept.top}{a logical value indicating whether the intercept (or constant) coefficients should be on top, rather than on the bottom, of the table.
%%     ~~Describe \code{intercept.top} here~~
}
  \item{label}{a character string containing the \code{\\label{}} TeX markers for the tables.
%%     ~~Describe \code{intercept.top} here~~
}
 
 \item{model.names}{a logical value indicating whether model names should be included in the table.
%%     ~~Describe \code{model.names} here~~
}
 \item{model.numbers}{a logical value indicating whether models should be numbered. No number is used whenever a regression table includes only one model.
%%     ~~Describe \code{model.numbers} here~~
}
  \item{notes}{a character string vector containing notes to be included below the table. The character strings can include special substrings that will be replaced by the corresponding cutoffs for statistical significance 'stars': \code{[*]}, \code{[**]}, and \code{[***]} will be replaced by the cutoffs, in percentage terms, for one, two and three 'stars,' respectively (e.g., 10, 5, and 1). Similarly, \code{[0.*]}, \code{[0.**]} and \code{[0.***]} will be replaced by the numerical value of cutoffs for one, two and three 'stars' (e.g., 0.1, 0.05, and 0.01). \code{[.*]}, \code{[.**]} and \code{[.***]} will omit the leading zeros (e.g., .1, .05, .01).
%%     ~~Describe \code{notes} here~~
}
  \item{notes.align}{a character string that specifies how notes should be aligned under the table. One of three strings can be used: \code{"l"} for left alignment, \code{"r"} for right alignment, and \code{"c"} for centering. This argument is not case-sensitive.
%%     ~~Describe \code{notes.align} here~~
}
  \item{notes.label}{a character string containing a label for the notes section of the table.
%%     ~~Describe \code{notes.label} here~~
}
  \item{omit}{a vector of regular expressions that specifies which of the explanatory variables should be omitted from presentation in the table. This argument might be used, for instance, to exclude fixed effects dummies from being presented. The default value of \code{NULL} means that no variables will be excluded.
%%     ~~Describe \code{omit} here~~
}
  \item{omit.labels}{a character string vector of labels that correspond to each of the regular expressions in \code{omit}, and that will be used in a sub-table that indicates whether variables have been omitted from a given model. \code{omit} and \code{omit.labels} must be equal in length.%%     ~~Describe \code{omit.labels} here~~
}
  \item{omit.stat}{a character string vector that specifies which statistics should be omitted from regression table output. For instance \code{omit = c("ll","rsq")} will omit the log-likelihood and the R squared statistics. See the \link{list of statistic codes}. This argument is not case-sensitive.%%     ~~Describe \code{omit.stat} here~~
} 
  \item{omit.yes.no}{a character string vector of length 2 that contains the 'yes' and 'no' strings to indicate whether, in any specific model, variables were omitted from the table, as specified by \code{"omit"}.
%%     ~~Describe \code{omit.yes.no} here~~
}
  \item{ord.intercepts}{a logical value indicating whether intercepts for models with ordered dependent variables (such as ordered probit, or ordered logit) are included in the table.
%%     ~~Describe \code{ord.intercepts} here~~
}
  \item{star.char}{character string to be used as the 'star' to denote statistical significance.
%%     ~~Describe \code{star.char} here~~
}
  \item{star.cutoffs}{a numeric vector with a length of at least \code{1} and at most \code{3} that indicates the statistical signficance cutoffs for one, two and three 'stars,' respectively. For elements with \code{NA} values, the corresponding 'star' will not be used.
%%     ~~Describe \code{star.cutoffs} here~~
}
  \item{zero.component}{a logical value indicating whether to report coefficients for the \code{zero} component of \code{zeroinfl} and \code{hurdle} estimation results. If \code{FALSE}, the \code{count} component is displayed.
%%     ~~Describe \code{zero.component} here~~
}
 \item{summary.logical}{a logical value indicating whether logical variables should be reported in summary statistics table. If so, they will be treated as if  they had values of 0 (corresponding to \code{FALSE}) and 1 (\code{TRUE}). 
%%     ~~Describe \code{zero.component} here~~
}
  \item{nobs}{a logical value that toggles whether the number of observations (N) for each variable is shown in summary statistics tables.}
  \item{mean.sd}{a logical value that toggles whether variable means and standard deviations are shown in summary statistics tables.}
  \item{min.max}{a logical value that toggles whether variable minima and maxima are shown in summary statistics tables.}
  \item{median}{a logical value that toggles whether variable medians are shown in summary statistics tables.}
  \item{iqr}{a logical value that toggles whether the 25th and 75th percentiles for each variable are shown in summary statistics tables. ('iqr' stands for interquartile range.)}

}
\details{Arguments with a value of \code{NULL} will use the default settings of the requested \code{style}.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\code{stargazer} uses \code{cat()} to output LaTeX code for the table. To allow for further processing of this output, \code{stargazer} also returns the same output invisibly as a character string vector. You can include the produced tables in your paper by inserting \code{stargazer} output into your publication's TeX source.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
%%  ~~who you are~~
Marek Hlavac   < hlavac at fas.harvard.edu > \cr
Harvard University, Political Economy and Government \cr
}
\section{Acknowledgments and New Features}{
I would like to thank everyone who has tested this package, or provided useful comments and suggestions. Please see \link{stargazer package acknowledgments}. \cr
\cr
See \link{stargazer news} for a list of new models and features in each release of stargazer.
}

\section{Please cite as:}{
Hlavac, Marek (2013). stargazer: LaTeX code for well-formatted regression and summary statistics tables.\cr
 R package version 3.0. http://CRAN.R-project.org/package=stargazer 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## create summary statistics table based for 'attitude' data frame
stargazer(attitude)

## list the content of the data frame 'attitude'
stargazer(attitude, summary=FALSE)

##  2 OLS models
linear.1 <- lm(rating ~ complaints + privileges + learning + raises + critical, data=attitude)
linear.2 <- lm(rating ~ complaints + privileges + learning, data=attitude)

## create an indicator dependent variable, and run a probit model
 
attitude$high.rating <- (attitude$rating > 70)
probit.model <- glm(high.rating ~ learning + critical + advance, data=attitude, family = binomial(link = "probit"))
 
stargazer(linear.1, linear.2, probit.model, title="Regression Results")

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
\keyword{ package }
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ univar }
\keyword{ interface }
\keyword{ programming }
\keyword{ IO }% __ONLY ONE__ keyword per line
