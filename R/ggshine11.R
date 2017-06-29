gQTLbrowse = function (store, 
    probeGR, stateGR, phenGR, AssocTag.in = "p.value", 
    AssocGen, orgDbObj = Homo.sapiens, AssocTag.out = "ml10p", 
    SNP.tag.in = "SNP", phen.snp.tag = "SNPs", extractTag = "gene", 
    chromStateTag = "name",
    datafilter = function(x) {
        seqlevelsStyle(x) = "UCSC"
        x
    }) 
#
# to call with gtexLite infrastructure
# library(gtexWbl); ww = wbLite()
# library(gwascat); data(ebicat37)
# gQTLbrowse( ww, gencodeV12, hmm878, ebicat37, AssocTag.in="pvalue",
#  phen.snp.tag="SNPS")
# 
#
{
    allsyms = keys(orgDbObj, keytype = "SYMBOL")
    availProbes = store@probemap$probeid
    availProbes = intersect(availProbes, names(probeGR))
    availSyms = probeGR[availProbes]$gene_name
    symok = which(availSyms %in% allsyms)
    availTypes = probeGR[availProbes]$gene_type
    p2g = availSyms = availSyms[symok]
    p2t = availTypes[symok]
    availProbes = availProbes[symok]
    names(p2g) = availProbes
    names(p2t) = availProbes
    ui = fluidPage(
          fluidRow(helpText("genomic QTL visualizer: mouseover for metadata")),
          fluidRow(selectInput("sym", "Gene symbol", 
        choices = sort(availSyms), multiple = FALSE)), fluidRow(verbatimTextOutput("ens_out")), 
          fluidRow(ggvisOutput("p"))
          )
    server = function(input, output) {
        output$ens_out = renderText(paste0("GENCODE V12 ENSEMBL ID: ", 
            availProbes[which(availSyms == input$sym)[1]]))
        filteredData = reactive({
            if (is.null(input$sym) || input$sym == "") 
                return(NULL)
            n1 = datafilter(extractByProbes(store, availProbes[which(availSyms == 
                input$sym)[1]], extractTag = extractTag))
            n1$stateAnno = rep("none", length(n1))
            fo = findOverlaps(n1, stateGR)
            n1$stateAnno[queryHits(fo)] = as.character(mcols(stateGR)[[chromStateTag]])[subjectHits(fo)]
            mcols(n1)[[AssocTag.out]] = pmin(7, -log10(mcols(n1)[[AssocTag.in]]))
            mydf <- data.frame(chr = as.character(seqnames(n1)), 
                pos = start(n1), probeid = mcols(n1)[[extractTag]], 
                snp = mcols(n1)[[SNP.tag.in]], assoc = mcols(n1)[[AssocTag.out]], 
                Mb = start(n1)/1e+06, stateAnno = n1$stateAnno, stringsAsFactors = FALSE)
            mydf$gene = as.character(p2g[mydf$probeid])
            mydf$type = as.character(p2t[mydf$probeid])
            mod = gmod2(input$sym)
            extra = tail(mydf, length(mod))
            extra$stateAnno = paste0("TXLOC(", input$sym, ")")
            extra$Mb = start(mod)/1e+06
            extra$assoc = -0.25
            extra$snp = NA
            mydf = rbind(mydf, extra)
            disdat = phenGR[which(mcols(phenGR)[[phen.snp.tag]] %in% 
                mydf$snp)]
            if (length(disdat) > 0) {
                extra2 = tail(mydf, length(disdat))
                extra2$assoc = -0.4
                extra2$Mb = start(disdat)/1e+06
                extra2$stateAnno = paste0("  trait: ", disdat$"DISEASE/TRAIT")
                mydf = rbind(mydf, extra2)
            }
            mydf$rowid = 1:nrow(mydf)
            mydf <<- mydf
            vals = mydf %>% dplyr::filter(gene == input$sym)
            return(vals)
        })
        P1 = reactive({
            if (is.null(filteredData)) 
                return(NULL)
            all_values <- function(x) {
                if (is.null(x)) 
                  return(NULL)
                row <- mydf[mydf$rowid == x$rowid, ]
                paste0(names(row), ": ", format(row), collapse = "<br />")
            }
            stateN = metadata(stateGR)$displayTag
            filteredData %>% ggvis(~Mb, ~assoc, `:=`(key, ~rowid), 
                fill = ~stateAnno) %>% layer_points() %>% add_tooltip(all_values, 
                "hover") %>% layer_points() %>% add_legend("fill", 
                title = paste0(stateN, " state"), values = unique(mydf$stateAnno)) %>% 
                add_axis("y", title = paste0("assoc (-log10p) w/ ", 
                  input$sym, " expr"))
        })
        P1 %>% bind_shiny("p")
    }
    shinyApp(ui = ui, server = server)
}
