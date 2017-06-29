
gmod2 = function (sym, genome = "hg19", orgDb, #=Homo.sapiens,
   collector=exonsBy, verbose=FALSE) 
{
    if (missing(orgDb)) {
       require("Homo.sapiens")
       orgDb = Homo.sapiens
       }   
    rend = suppressPackageStartupMessages
    if (verbose) rend = force
    rend({
    require(txn <- gsub("%%G%%", genome, "TxDb.Hsapiens.UCSC.%%G%%.knownGene"),
      character.only=TRUE)
    require(Homo.sapiens)
    })  
    txdb = get(txn)
    suppressWarnings({
    num = AnnotationDbi::select(orgDb, keys=sym, keytype="SYMBOL",
          columns="ENTREZID")$ENTREZID
    })
    collector(txdb, by = "gene")[[num]]
}

