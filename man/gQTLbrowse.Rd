\name{gQTLbrowse}
\alias{gQTLbrowse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
interactive exploration of eQTL results
}
\description{
interactive exploration of eQTL results; code as of jan 6 2015 is
written with selectizeInput to throw an error on gene symbol deletion
}
\usage{
gQTLbrowse(store, baseSE, stateGR, phenGR, FDRsupp, orgDbObj = Homo.sapiens,
  selector=selectizeInput)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{store}{
instance of \code{\link[gQTLBase]{ciseStore-class}}
}
  \item{baseSE}{
a \code{\link[GenomicRanges]{SummarizedExperiment-class}} instance
that will be used to get gene-level metadata; the rowData component
should include fields \code{gene\_name}, \code{gene\_type}
}
  \item{stateGR}{
a \code{\link[GenomicRanges]{GRanges-class}} instance with metadata field
\code{name} that is used to define `state' of range; used with chromHmm
labeling by default
}
  \item{phenGR}{
a GRanges instance with disease association labeling in a metadata
field \code{Disease.Trait}; currently the examples depend on the
GEUVADIS data and the GEUVADIS SNP identifiers are assumed given
in metadata field \code{geuvvid}
}
  \item{FDRsupp}{
an instance of \code{\link[gQTLstats]{FDRsupp-class}} that provides
a function mapping from association scores to FDR estimates
}
  \item{orgDbObj}{
an AnnotationDbi OrgDb instance such as \code{Homo.sapiens} to be
used for transcript location
}
}
\details{
currently coded to simplify exposure of an error when
selectizeInput is used
}
\value{
side effects: will open browser window and allow interactive
gene selection and hover-over tool tip for locus- or trait-specific
metadata reporting
}
%\references{
%%% ~put references to the literature/web site here ~
%}
%\author{
%%%  ~~who you are~~
%}
%\note{
%%%  ~~further notes~~
%}
%
%%% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
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
#
# will throw error when gene symbol in input box is deleted
# but will simply select properly if it is not; change selector to
# selectInput for simpler interface
#
print( gQTLbrowse( prst, geuFPKM, stateGR = hmm878, 
     phenGR = ggr19, FDRsupp = filtFDR, selector=selectizeInput ) )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ graphics }
