#!/usr/bin/env Rscript

#####################
##### LOAD LIBS #####
#####################

suppressPackageStartupMessages({
    library("here")
    library("glue")
    library("knitr")
    library("tidyverse")
    library("cli")
    library("lmerTest")
    library("easystats")
    library("ggthemes")
    library("svglite")
    library("marginaleffects")
})

# some options
options(width=180)


#####################
##### LOAD FUNS #####
#####################

# function to not print some code to screen
dont_print <- function(expr) {
    invisible(force(expr))
}

# function to fourth power
power4 <- function(x) {
    res <- x^4
    return(res)
}


# function to square the sine
sqsine <- function(x) {
    res <- sin(x)^2
    return(res)
}


# shortcut function to cli report with spaces
report <- function(text,type) {
    if(type=="s") {
        cli::cli_text("")
        cli::cli_text("")
        cli::cli_alert_success(text)
        cli::cli_text("")
    } else if (type=="i") {
        cli::cli_text("")
        cli::cli_text("")
        cli::cli_text("")
        cli::cli_alert_info(text)
        cli::cli_text("")
    } else if (type=="g") {
        cli::cli_text("")
        cli::cli_text("")
        cli::cli_text(text)
    } 
}


# function to summarise qpcr by grouping var ()
summarise_qpcr <- function(df,grouping) {
    df.prop <- df |> 
        dplyr::group_by({{grouping}}) |>
        dplyr::mutate(copiesLitreTotal=sum(copiesLitre)) |>
        dplyr::ungroup() |>
        dplyr::group_by(assay,{{grouping}},copiesLitreTotal) |>
        dplyr::summarise(copiesLitreSpecies=sum(copiesLitre),
            nFilterReps=length(unique(extractionTubeLabel)),
            nFilterRepsPositive=length(unique(extractionTubeLabel[nPcrRepsPositive>0])),
            nPcrReps=n(),
            nPcrRepsPositive=length(which(copiesLitre>0)),
            .groups="drop") |>
        dplyr::mutate(proportion=copiesLitreSpecies/copiesLitreTotal) |> 
        dplyr::arrange({{grouping}},assay) |>
        dplyr::select({{grouping}},assay,copiesLitreSpecies,proportion,nFilterReps,nFilterRepsPositive,nPcrReps,nPcrRepsPositive) |>
        dplyr::rename(group={{grouping}},species=assay,copiesTotal=copiesLitreSpecies) |>
        dplyr::mutate(proportion=dplyr::if_else(is.nan(proportion),0,proportion))
return(df.prop)
} 


# summarise tissues function
summarise_tissues <- function(df,grouping) {
    df.bin <- df |> 
        dplyr::select({{grouping}},specificEpithet) |>
        dplyr::group_by({{grouping}},specificEpithet) |>
        dplyr::count() |> 
        tidyr::pivot_wider(names_from=specificEpithet,values_from=n,values_fill=0) |>
        dplyr::mutate(across(c(niloticus,urolepis,leucostictus),~if_else(.x>0,1,0))) |>
        dplyr::ungroup() |>
        dplyr::rename(group={{grouping}})
return(df.bin)
}
