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
The \code{stargazer} command produces LaTeX code for well-formatted tables that hold regression analysis results from several models side-by-side, as well as summary statistics. It supports model objects from \code{lm}, \code{glm}, \code{svyglm}, \code{plm}, \code{betareg}, \code{gee}, \code{gam}, \code{polr}, \code{survreg}, \code{coxph}, \code{tobit} (AER), \code{ivreg} (AER), \code{zeroinfl} (pscl), \code{hurdle} (pscl), \code{multinom} (nnet) as well as from the implementation of these in \code{zelig}. It also supports the following \code{zelig} models for social network analysis: \code{"cloglog.net"}, \code{"gamma.net"}, \code{"probit.net"} and \code{"logit.net"}.
}
\usage{
stargazer(  ..., 
            title = "", style = "default", summary = TRUE,
            covariate.labels = NULL, dep.var.labels = NULL,
            align = FALSE, decimal.mark = NULL,
            digit.separate = NULL, digit.separator = NULL,
            digits = NULL, digits.extra = NULL, initial.zero = NULL, intercept.top = NULL, 
            model.names = NULL, model.numbers = NULL, 
            notes = NULL, notes.align = NULL, notes.label = NULL, 
            omit = NULL, omit.labels = NULL, omit.stat = NULL, omit.yes.no = c("Yes", "No"), 
            ord.intercepts = FALSE, star.char = NULL, star.cutoffs = NULL, zero.component = FALSE,
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

  \item{decimal.mark}{character string that will serve as the decimal mark. For instance, the string \code{","} will represent decimal commas, while \code{"."} means tables will use decimal points.}


  \item{digit.separate}{a numerical vector that indicates where digit separators should be placed. The first element of the vector indicates the number of digits (counted from the decimal mark to the left) that will be separated. The second element indicates the number of digits that will be separated from that 'first' separator, and so on. A value of \code{3} corresponds to a thousands separator, while a value of \code{0} indicates no separation. Alternatively, \code{digit.separate} can be one of the following character strings: \code{"lakh"} (equivalent to \code{c(4,3)}), \code{"china"} or \code{"japan"} (both equivalent to a value of \code{4}).
%%     ~~Describe \code{digit.separate} here~~
}
  \item{digit.separator}{character string that will serve as the digit (e.g., thousands) separator. Commonly used strings include \code{","} for a comma separator, \code{" "} for a single space separator, and \code{""} for no separation.
%%     ~~Describe \code{digit.separator} here~~
}
  \item{digits}{integer that indicates how many decimal places should be used.
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
Marek Hlavac
<hlavac at fas.harvard.edu>
}
\section{Acknowledgments}{
I would like to thank everyone who has tested this package, or provided useful comments and suggestions. Please see \link{stargazer package acknowledgments}.
}

\section{Please cite as:}{
 Marek Hlavac (2013). stargazer: LaTeX code for well-formatted regression and summary statistics tables.\cr
 R package version 2.0.1. http://CRAN.R-project.org/package=stargazer 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Several packages by other authors provide functionality similar to that of \code{stargazer}. These include \code{xtable} by David B. Dahl and \code{apsrtable} by Michael Malecki.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## create summary statistics table based for 'attitude' data frame
stargazer(attitude)

## list the content of the data frame 'attitude'
stargazer(attitude, summary=FALSE)

## estimate and report results of 3 OLS models
m1 <- lm(rating ~ complaints + privileges + learning + raises + critical, data=attitude)
m2 <- lm(rating ~ complaints + privileges + learning, data=attitude)
m3 <- lm(rating ~ learning + critical + advance, data=attitude)

stargazer(m1, m2, m3)

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
