Rscript runChicago.R --input-files raw/PCHIC_Data/GM_rep1.chinput  raw/PCHIC_Data/GM_rep2.chinput raw/PCHIC_Data/GM_rep3.chinput --output-prefix output --design-dir raw/Design --en-feat-list raw/Features/featuresGM.txt --export-format output_washU_text

1.2: Enrichments are suggesting the significant number of overlaps for each feature mark with the bait. It makes sense because the interactions are with histone modification markers for heterochromatin or euchromatin. 

promoter-promoter:
      chrom chromStart  chromEnd name score  value exp color sourceChrom sourceStart sourceEnd sourceName sourceStrand targetChrom targetStart targetEnd targetName targetStrand
1645  chr20   44438565  44565593    .  1000  34.77   .     0          20    44562442  44565593     414690            +          20    44438565  44442365     414672            +
1655  chr20   44438565  44607204    .   986  34.29   .     0          20    44596299  44607204     414697            +          20    44438565  44442365     414672            +
2819  chr21   26837918  26939577    .   978  34.02   .     0          21    26837918  26842640     423425            +          21    26926437  26939577     423450            +
1646  chr20   44452862  44565593    .   974  33.89   .     0          20    44562442  44565593     414690            +          20    44452862  44471524     414675            +
475   chr20   17660712  17951709    .   973  33.85   .     0          20    17946510  17951709     408664            +          20    17660712  17672229     408599            +
577   chr20   24972345  25043735    .   973  33.84   .     0          20    24972345  24985047     410679            +          20    25036380  25043735     410693            +

promoter-enhancer:
      chrom chromStart  chromEnd name score  value exp color sourceChrom sourceStart sourceEnd sourceName sourceStrand targetChrom targetStart targetEnd targetName targetStrand
2842  chr21   26797667  26939577    .   952  33.13   .     0          21    26926437  26939577     423450            +       chr21    26797667  26799364          .            -
2254  chr20   55957140  56074932    .   928  32.29   .     0          20    55957140  55973022     417632            +       chr20    56067414  56074932          .            -
2838  chr21   26790966  26939577    .   838  29.17   .     0          21    26926437  26939577     423450            +       chr21    26790966  26793953          .            -
231   chr20    5585992   5628028    .   830  28.88   .     0          20     5585992   5601172     404738            +       chr20     5625693   5628028          .            -
2839  chr21   26793954  26939577    .   754  26.23   .     0          21    26926437  26939577     423450            +       chr21    26793954  26795680          .            -
278   chr20    5515866   5933156    .   750  26.08   .     0          20     5929472   5933156     404832            +       chr20     5515866   5523933          .            


gene1: LINC00158 encodes a six-lncRNA signature based on a competing endogenous RNA network for predicting the risk of tumour recurrence in bladder cancer patient, which makes sense since this is a cancer cell line.

gene2: MIR155HG The long RNA transcribed from this gene is expressed at high levels in lymphoma and may function as an oncogene. This makes sense since this is a lymnphoma cell line

