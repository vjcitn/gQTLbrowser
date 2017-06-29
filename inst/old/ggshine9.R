
gQTLbrowse = function( store, baseSE, 
   stateGR, phenGR, FDRsupp, orgDbObj=Homo.sapiens, selector=selectizeInput ) {
#
# interface to shiny/ggvis eqtl exploration
# assumes central identifier is the gene symbol
#
   #
   # get all available gene symbols
   #
   
   allsyms = keys(orgDbObj, keytype="SYMBOL")
   
#
#  filter the symbols relevant to the input SummarizedExperiment baseSE,
#  which will generally not use gene symbols as rownames ... in fact
#  FIXME -- we are assuming existence of gene_name and gene_type in
#  rowData, as for geuvPack geuFPKM
#
#  this sort of symbol mapping exercise is common and should be
#  substantially abstracted -- perhaps a mapping between GSEABase entities
#
   
   availProbes = store@probemap$probeid
   availSyms = rowData(baseSE[ availProbes, ])$gene_name
   sorts = sort(availSyms)
   symok = which(availSyms %in% allsyms)
   availTypes = rowData(baseSE[ availProbes, ])$gene_type
   p2g = availSyms = availSyms[symok]
   p2t = availTypes[symok]
   availProbes = availProbes[symok]
   names(p2g) = availProbes
   names(p2t) = availProbes
   
   #
   # we want to get the symbols from shiny, but use
   # ggvis for tooltips
   #
   
   ui = fluidPage(
   #
   # FIXME should be selectize
   #
      fluidRow(selector('sym', 'Gene symbol', choices=c("", sorts), 
             selected=sorts[2],
             multiple=FALSE)), 
      fluidRow(verbatimTextOutput('ens_out')) ,
      fluidRow( ggvisOutput('p') )
      )
   
   server = function(input, output) {
   
   ## first, extract the ENSEMBL ID for selected symbol
   
      output$ens_out = renderText(
                         paste0("GEUVADIS ENSEMBL ID: ",
                            availProbes[ which(availSyms == input$sym)[1] ] ))
   
   ## second, acquire the appropriate eQTL testing results
   ##   and bind information on chromatin state, transcript location
   
      filteredData = reactive( {
        validate( 
          need( input$sym != "", "provide symbol" )
        )
   #
   # get the GRanges with eQTL results
   #
        n1 = extractByProbes( store, 
                  availProbes[ which(availSyms == input$sym)[1] ] )
   #
   # obtain chromatin state labels
   #
        n1$st878 = rep("none", length(n1))
        fo = findOverlaps(n1, stateGR)
        n1$st878[ queryHits(fo) ] = as.character(stateGR$name)[ subjectHits(fo) ]
   #     uniqst = unique(stateGR$name) # useless effort at persistent colormap
   #     nuniqst = length(uniqst)
   #     cmap = colorRampPalette(c("red", "blue"))(nuniqst) # vector of codes
   #     n1$col878 = cmap[as.numeric(factor(n1$st878))]
   #
   # execute the FDR filter
   #
        n1 = FDRsupp@filterUsed(n1)
   #
   # compute FDR
   #
        n1$ml10FDR = pmin(6, -log10(getFDRfunc(FDRsupp)(n1$chisq)))
   #     n1 <<- n1
   #
   # build data frame for visualization
   #
        mydf <- data.frame(chr=as.character(seqnames(n1)), pos=start(n1),
     	  MAF = n1$MAF,
               probeid=n1$probeid, snp=n1$snp, ml10FDR = n1$ml10FDR,
               stringsAsFactors=FALSE, Mb=start(n1)/1e6, st878=n1$st878)
   #
   # use global maps to recover symbol and GEUVADIS 'gene type'
   #
        mydf$gene = as.character( p2g[ mydf$probeid ] )
        mydf$type = as.character( p2t[ mydf$probeid ] )
   #
   # obtain a 'gene model' for the selected symbol, so that
   # locations of transcripts can be given
   #
        mod = gmod2( input$sym )
   #
   # add the location information, "faking" fields for eQTL results
   #
        extra = tail(mydf, length(mod))
        extra$st878 = paste0("TXLOC(", input$sym,")")
        extra$Mb = start(mod)/1e6
        extra$ml10FDR = -.25
        extra$snp = extra$MAF = NA
        mydf = rbind(mydf, extra)
   #
   # get disease loci
   #
        disdat = phenGR[ which(phenGR$geuvvid %in% mydf$snp) ]
        if (length(disdat) > 0) {
           extra2 = tail(mydf, length(disdat))
           extra2$ml10FDR = -.4
           extra2$MAF = NA
           extra2$Mb = start(disdat)/1e6
           extra2$st878 = paste0("  trait: ", disdat$Disease.Trait)
           mydf = rbind(mydf, extra2)
           }
   #     mydf = mydf[ order(mydf$st878), ]
   #
   # construct the key for tooltip
   #
        mydf$rowid = 1:nrow(mydf)
        mydf <<- mydf
        vals = mydf %>% dplyr::filter( gene == input$sym )  # do earlier !?!
        return(vals)
        } )
   
      P1 = reactive( {
        validate( 
          need( input$sym != "", "provide symbol" )
        )
         all_values <- function(x) {
             if(is.null(x)) return(NULL)
             row <- mydf[mydf$rowid == x$rowid, ]
             paste0(names(row), ": ", format(row), collapse = "<br />")
           }
         filteredData %>% ggvis(~Mb, ~ml10FDR, key := ~rowid,
                  fill = ~st878) %>% 
               layer_points() %>%
               add_tooltip(all_values, "hover") %>% layer_points() %>%
               add_legend("fill", title=paste0(deparse(substitute(stateGR)), " state"), values=unique(mydf$st878)) %>% 
   
               add_axis("y", title=paste0("-log10 FDR assoc w/ ", input$sym, " expr" ))
         } )
      P1 %>% bind_shiny("p")
   }
   
   shinyApp(ui = ui, server=server)
}
