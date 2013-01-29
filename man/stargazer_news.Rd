\name{stargazer_news}
\alias{stargazer news}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{What's new?}
%%  ~~function to do ... ~~
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This page summarizes the models, features, and bug fixes that were introduced in each release of \code{stargazer}. Please note that I am always looking for comments and suggestions. Do not hesitate to contact me at hlavac [at] fas.harvard.edu.
}

\section{What's new?}{
v. 3.0.1: \cr
\cr
- Fixed bugs: \cr
-- corrected an issue that led to problems when printing data frames with multiple decimal places\cr
-- some reporting and formatting changes to summary statistics and data frame tables \cr
\cr
v. 3.0: \cr
\cr
- New models: \cr
-- clm (ordinal) \cr
-- clogit (survival) \cr
-- ergm (ergm) \cr
-- glmer (lme4) \cr
-- gls (nlme) \cr
-- lmer (lme4) \cr
-- lmrob (robustbase) \cr
-- nlmer (lme4) \cr
-- pmg (plm) \cr
-- rlm (MASS) \cr
\cr
- New features: \cr
-- users can customize coefficients (argument \code{coef}), standard errors (\code{se}), test statistics (\code{t}), and p-values (\code{p}) \cr
-- automatic calculation of z-scores and p-values when the user supplies custom standard errors (can be toggled using arguments \code{t.auto} and \code{p.auto}) \cr
-- ability to set \code{\\label{}} markers in TeX for each table (using argument \code{label}) \cr
-- summary statistics table can now report logical (i.e., dummy, indicator) variables as if they had values of 0 (corresponding to \code{FALSE}) and 1 (\code{TRUE}). See argument \code{summary.logical}. \cr
-- user can choose between rounding to a set number of decimal places (argument \code{digits}) and reporting all available decimal places (by setting \code{digits} equal to \code{NA}) \cr
-- can omit all test statistics \cr
\cr
- Fixed bugs: \cr
-- package runs much faster \cr
-- corrected p-values for polr() and zelig oprobit, ologit models \cr
-- coxph now, by default, reports robust standard errors, for consistency with its summary() output \cr
-- multinom models now report the first set of coefficients (when multiple sets are present), instead of producing an error\cr
-- singularities in model regression no longer disrupt the output \cr
-- model formulas can now be symbols \cr
-- no more warnings when creating summary statistics tables \cr
-- minor formatting changes \cr
\cr
v. 2.0.1: \cr
\cr
- Fixed bugs: \cr
-- stargazer now reports the correct number of observations for plm() models \cr
\cr
v. 2.0: \cr
\cr
- New models: \cr
-- betareg (betareg) \cr
-- hurdle (pscl) \cr
-- ivreg (AER) \cr
-- multinom (nnet) \cr
-- plm (plm) \cr
-- tobit (AER) \cr
-- zeroinfl (pscl) \cr
    \cr
- New features: \cr
-- direct output of data frames into LaTeX (\code{summary = FALSE}) \cr
-- ability to omit the reporting of selected statistics (using argument \code{omit.stat})
-- alignment of columns at the decimal mark (\code{align = TRUE}) \cr
-- automatic coordination of star cutoffs and values in regression table notes \cr
\cr
- Fixed bugs: \cr
-- argument \code{digits} now works properly \cr
-- \code{stargazer} updated to reflect, and work properly with, recent changes to \code{Zelig} \cr
-- variable names can now contain dollar signs and underscores \cr
-- some minor formatting fixes \cr
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ package }
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ univar }
\keyword{ interface }
\keyword{ programming }
\keyword{ IO }% __ONLY ONE__ keyword per line
