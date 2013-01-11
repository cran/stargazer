\name{stargazer_acknowledge}
\alias{stargazer package acknowledgments}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Acknowledgments}
%%  ~~function to do ... ~~
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
I would like to thank everyone who has tested this package, or provided useful suggestions. I am especially grateful to the following people, listed in alphabetical order (based on last names):
}

\section{Thank you:}{
Erin Baggott, Harvard University, USA\cr
Ruben de Bliek, Erasmus University Rotterdam, the Netherlands\cr
Volha Charnych, Harvard University, USA\cr
John Coglianese, Harvard University, USA\cr
Stefan Dimitriadis, Harvard University, USA\cr
Bernard Fraga, Harvard University, USA\cr
Cheng Gao, Harvard University, USA\cr
Rebecca Goldstein, Harvard University, USA\cr
Connor Huff, Harvard University, USA\cr
Melissa Kline, Massachussetts Institute of Technology, USA\cr
Akos Lada, Harvard University, USA\cr
Carl Jacob Liebersohn, Massachussetts Institute of Technology, USA\cr
Daniel Yew Mao Lim, Harvard University, USA\cr
Christopher Lucas, Harvard University, USA\cr
Tamim Mohammad, University of Massachussetts-Boston, USA\cr
Stephen Pettigrew, Harvard University, USA\cr
Giuseppe Ragusa, Luiss University, Italy\cr
Jason Sclar, Harvard University, USA\cr
Jennifer Sheehy-Skeffington, Harvard University, USA\cr
Emily Stephen, Boston University, USA\cr
Beth Truesdale, Harvard University, USA\cr
Anna Weisfeiler, University of Wisconsin-Madison, USA\cr
Steffen Zittlau, University of Mannheim, Germany\cr
}

\section{What's new?}{
v. 2.0.1: \cr
\cr
- Fixed bugs: \cr
-- stargazer now reports the correct number of observations for plm() models \cr
\cr
v. 2.0: \cr
\cr
- New models: \cr
-- plm (plm) \cr
-- tobit (AER) \cr
-- ivreg (AER) \cr
-- multinom (nnet) \cr
-- betareg (betareg) \cr
-- zeroinfl (pscl) \cr
-- hurdle (pscl) \cr
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
