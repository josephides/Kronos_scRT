#parse input
suppressPackageStartupMessages(library(optparse, quietly = TRUE))

# Sys.setenv('R_MAX_VSIZE'=32000000000)

options(stringsAsFactors = FALSE,
        dplyr.summarise.inform=FALSE,
        warn = 1,
        scipen = 999)

option_list = list(
  make_option(
    c("-F", "--File"),
    type = "character",
    default = NULL,
    help = "Replication timing files separated by a comma. Format: chr <TAB> start <TAB> end <TAB> group",
    metavar = "character"
  ),
  make_option(
    c("-s", "--sort"),
    type = "character",
    default = NULL,
    help = "Group names orders",
    metavar = "character"
  ),
  make_option(
    c("-o", "--out"),
    type = "character",
    default = "output",
    help = "Output directory [default= %default]",
    metavar = "character"
  ),
  make_option(
    c("-f", "--output_file_base_name"),
    type = "character",
    default = "out",
    help = "Base name for the output file [default= %default]",
    metavar = "character"
  )
)

#recover inputs
opt = parse_args(object = OptionParser(option_list = option_list))

#load module
suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
suppressPackageStartupMessages(library(GGally, quietly = TRUE))
suppressPackageStartupMessages(library(ggcorrplot, quietly = TRUE))
suppressPackageStartupMessages(library(foreach, quietly = TRUE))

theme_set(theme_linedraw())


if (!'File' %in% names(opt)) {
  stop("RT files were not provided. See script usage (--help)")
  
} else{
  opt$File = unlist(str_split(opt$File, ',') )
}
if ('sort' %in% names(opt)) {
  opt$sort = str_split(opt$sort, ',')[[1]]
  
}

#create directory
if (!dir.exists(opt$out)) {
  dir.create(opt$out,recursive = T)
}

scRT = foreach(
  i = 1:length(opt$File),
  .packages = 'tidyverse',
  .combine = 'rbind'
) %do% {
  tmp = read_tsv(opt$File[i], col_types = cols(chr = 'c'))
  if ('sort' %in% names(opt)) {
    tmp %>%
      mutate(group = factor(group, levels = opt$sort))
    
  } else{
    tmp
  }
}

scRT = scRT[,c("chr","start","end",'RT',"group")]

scRT$group <- gsub("_", " ", scRT$group)
scRT$group <- gsub(" $", "", scRT$group)

old_names = unique(scRT$group)

scRT = scRT %>% spread(group, RT) %>%
  drop_na() %>%
  dplyr::select(-chr, -start, -end)

new_names =  c( "hTERT-RPE1 Connolly2022 (retina)",
                     "HCT-116 DKO1 (colon)",
                     "HCT-116 WT S1 (colon)",
                     "HCT-116 WT S2 (colon)",
                     "HeLa Gnan2022 (cervix)",
                     "JEFF (lymphocyte)",
                     "MCF-7 Gnan2022 S1 (breast)",
                     "MCF-7 Gnan2022 S2 (breast)",
                     "184-hTERT SA039 S1 (breast)",
                     "184-hTERT SA039 S2 (breast)",
                     "184-hTERT SA039 S3 (breast)",
                     "184-hTERT SA1101 S1 (breast)",
                     "184-hTERT SA906 (breast)",
                     "ERpos-PDX SA532X2 S2 (breast)",
                     "ERpos-PDX SA532X4 (breast)",
                     "FNA-TNBC SA1135 S2 (breast)",
                     "FNA-TNBC SA1135 S4 (breast)",
                     "GM18507 (lymphocyte)",
                     "HGSOC-OV2295 S1 (ovary)",
                     "HGSOC-OV2295 S2 (ovary)",
                     "HeLa Laks2019 (cervix)",
                     "T-47D S1 (breast)",
                     "TNBC-PDX SA501X11 S2 (breast)",
                     "TNBC-PDX SA501X2 S1 (breast)",
                     "TNBC-PDX SA501X5 (breast)",
                     "TNBC-PDX SA501X6 (breast)",
                     "TNBC-PDX SA535X5 S1 (breast)",
                     "TNBC-PDX SA535X8 S1 (breast)",
                     "TNBC-PDX SA604X6 (breast)",
                     "TNBC-PDX SA609X4 S2 (breast)",
                     "TNBC-PDX SA609X5 S1 (breast)",
                     "TNBC-PDX SA609X6 B01898 (breast)",
                     "TNBC-PDX SA609X6 B01899 S2 (breast)",
                     "TNBC-PDX SA609X7 (breast)",
                     "GM12878 (lymphocyte)",
                     "GM12891 (lymphocyte)",
                     "GM12892 (lymphocyte)",
                     "H7 (hESC)",
                     "MCF-7 Massey2022 (breast)",
                     "RKO (colon)",
                     "hTERT-RPE1 Takahashi2019 (retina)" )

colnames(scRT) <- new_names[match(colnames(scRT), old_names)]

print(colnames(scRT))

plot = ggcorrplot(
  scRT %>%
    cor(method = 'spearman'),
  lab = T,
  lab_col = 'white',
  outline.color = 'white',
  ggtheme = ggplot2::theme_bw , 
  tl.cex = 11, lab_size=3, pch.cex=10, hc.order = T, hc.method='ward.D2'
) +  scale_fill_gradient2(breaks=c(0, 1), limit=c(0, 1),
                          low='#00008B', mid='#1A85FF', high='#D41159',
                          midpoint = 0.5, name= 'Spearman\ncorrelation') +
  theme(text=element_text(size=12, family="mono"), legend.position="top", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggsave(plot = plot,
       filename = paste0(opt$out,'/Rplot_spearman_correlation.pdf'),
       width = 15,
       height = 16)

print('done')
