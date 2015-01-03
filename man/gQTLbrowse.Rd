\name{gQTLbrowse}
\alias{gQTLbrowse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
interactive exploration of eQTL results
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gQTLbrowse(store, baseSE, stateGR, phenGR, FDRsupp, orgDbObj = Homo.sapiens)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{store}{
%%     ~~Describe \code{store} here~~
}
  \item{baseSE}{
%%     ~~Describe \code{baseSE} here~~
}
  \item{stateGR}{
%%     ~~Describe \code{stateGR} here~~
}
  \item{phenGR}{
%%     ~~Describe \code{phenGR} here~~
}
  \item{FDRsupp}{
%%     ~~Describe \code{FDRsupp} here~~
}
  \item{orgDbObj}{
%%     ~~Describe \code{orgDbObj} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(ggr19)
library(GenomeInfoDb)
seqlevelsStyle(ggr19) = "NCBI"
library(geuvPack)
data(geuFPKM)
library(gQTLBase)
library(geuvStore)
pr = partialRegistry()
prst = ciseStore(pr, TRUE, FALSE)
data(ggr19)
data(filtFDR)
library(gQTLstats)
data(hmm878)
seqlevelsStyle(hmm878) = "NCBI"
library(Homo.sapiens)
print( gQTLbrowse( prst, geuFPKM, stateGR = hmm878, 
     phenGR = ggr19, FDRsupp = filtFDR ) )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ graphics }
