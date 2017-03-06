## RscriptEcho.R
args <- commandArgs(TRUE)

if (any(grep("Rout",args))){
  srcFile <- args[1]
  outFile <- args[2]
  args <- args[-c(1:2)]
} else {
  args <- commandArgs(TRUE)
  srcFile <- args[1]
  outFile <- paste0(make.names(date()), ".Rout")
  args <- args[-1]
}

sink(outFile, split = TRUE)
source(srcFile, echo = TRUE)
