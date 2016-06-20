\name{aggregateList}
\alias{aggregateListByAllele}
\alias{aggregateListByPosition}
\title{Aggregate variant lists}
\description{Generate lists of variants for input to association tests}
\usage{
aggregateListByAllele(gds, variants, indexOnly=FALSE)
aggregateListByPosition(gds, groups, indexOnly=FALSE)
}
\arguments{
  \item{gds}{A \code{\link[SeqArray]{SeqVarGDSClass}} object}
  \item{variants}{A data.frame of variants with columns "group_id", "chromosome", "position", "ref", "alt".}
  \item{groups}{A data.frame of groups with column "group_id", "chromosome", "start", "end".}
  \item{indexOnly}{Logical for whether to return only the "variant.id" and "allele.index" columns in the output (see Value).}
}
\details{
These functions produce output suitable for providing to \code{\link[GENESIS]{assocTestSeq}} in the \pkg{\link[GENESIS]{GENESIS}} package.
}
\value{
A list of data frames, one for each group. A single variant may have multiple rows if multiple alternate alleles are selected. Each data frame contains the following columns:
  \item{variant.id}{Unique identifier for the variant}
  \item{chromosome}{Chromosome}
  \item{position}{Position in base pairs}
  \item{ref}{Reference allele}
  \item{nAlleles}{Total number of alleles for this variant}
  \item{allele}{Alternate allele}
  \item{allele.index}{Integer index of this allele (1=first alternate, 2=second alternate, etc.)}
}
\author{Stephanie Gogarten}
%\seealso{}
\examples{
library(SeqVarTools)
gds <- seqOpen(seqExampleFileName())
seqSetFilter(gds, variant.sel=seqGetData(gds, "chromosome") == 22)
variants <- data.frame(chromosome=seqGetData(gds, "chromosome"),
                       position=seqGetData(gds, "position"),
                       ref=refChar(gds),
                       alt=altChar(gds, n=1),
                       stringsAsFactors=FALSE)
variants$group_id <- sample(LETTERS[1:2], nrow(variants), replace=TRUE)
aggregateListByAllele(gds, variants)

groups <- data.frame(group_id=LETTERS[1:2],
                     chromosome=22,
                     start=c(16000000, 2900000), 
                     end=c(30000000, 49000000),
		     stringsAsFactors=FALSE)
aggregateListByPosition(gds, groups)

seqClose(gds)
}
\keyword{utilities}
