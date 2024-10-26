suppressPackageStartupMessages(library(optparse, quietly = TRUE))

options(stringsAsFactors = FALSE,
        dplyr.summarise.inform=FALSE,
        warn = 1,
        scipen = 999)

option_list = list(
    make_option(
        c("-F", "--file"),
        type = "character",
        default = NULL,
        help = "Per cell stat file , if multiple files are provided they have to be separated by a comma",
        metavar = "character"
    ),
    make_option(
        c("-T", "--tracks"),
        type = "character",
        default = NULL,
        help = "Tracks file,  if multiple files are provided they have to be separated by a comma",
        metavar = "character"
    ),
    make_option(
        c("-o", "--out"),
        type = "character",
        default = "output",
        help = "Output directory [default= %default]",
        metavar = "character"
    )
)

#recover inputs
opt = parse_args(object = OptionParser(option_list = option_list))

suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))

#create directory

if(!dir.exists(opt$out)){
    dir.create(opt$out,recursive = T)
}

#load files
opt$file = str_split(opt$file, ',')[[1]]

opt$tracks = str_split(opt$tracks, ',')[[1]]


if ('file' %in% names(opt)) {
    for (file in opt$file) {
        #check if the file exists
        if (!file.exists(file)) {
            warning(paste0(file, " does not exist"))
            
            # check if it is the right file
        } else if (all(
            c(
                'cell_id',
                'normalized_dimapd',
                'mean_ploidy',
                'ploidy_confidence',
                'is_high_dimapd',
                'is_noisy',
                'effective_reads_per_1Mbp'
            ) %in% colnames(tryCatch(
                expr =  read_csv(file,
                                 col_types = cols(),
                                 n_max = 0),
                error = function(x)
                    tibble()
            ))
        )) {
            read_csv(file, col_types = cols()) %>%
                select(
                    cell_id,
                    normalized_dimapd,
                    mean_ploidy,
                    ploidy_confidence,
                    is_high_dimapd,
                    is_noisy,
                    effective_reads_per_1Mbp
                ) %>%
                `colnames<-`(
                    c(
                        'Cell',
                        'normalized_dimapd',
                        'mean_ploidy',
                        'ploidy_confidence',
                        'is_high_dimapd',
                        'is_noisy',
                        'coverage_per_1Mbp'
                    )
                ) %>%
                write_csv(paste0(file.path(opt$out, 'Kronos_format_'), basename(file)))
        } else{
            warning(paste0(file, " is not a per cell stat file"))
        }
    }
}
if ('tracks' %in% names(opt)) {
    for (tracks in opt$tracks) {
        #check if the file exists
        if (!file.exists(tracks)) {
            warning(paste0(tracks, " does not exist"))
            
            # check if it is the right file
        } else if (all(c('id', '#chrom', 'start', 'end', 'copy_number') %in% colnames(tryCatch(
            expr =  read_tsv(
                tracks,
                skip = 2,
                col_types = cols(`#chrom` = 'c'),
                n_max = 0
            ),
            error = function(x)
                tibble()
        )))) {
            #if yes proceed
            read_tsv(tracks,
                     skip = 2,
                     col_types = cols(`#chrom` = 'c')) %>%
                select(id, `#chrom`, start, end, copy_number) %>%
                `colnames<-`(c('Cell', 'chr', 'start', 'end', 'copy_number')) %>%
                mutate(reads = '10X') %>%
                write_tsv(paste0(file.path(opt$out, 'Kronos_format_'), basename(tracks)))
        } else{
            warning(paste0(tracks, " is not a proper track file"))
        }
    }
}

if (!('tracks' %in% names(opt) & 'file' %in% names(opt))) {
    stop('No input')
}

print('done')
