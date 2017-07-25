library(ypages)
library(rmarkdown)

setwd("src")

files <- list.files(pattern="rmd$")
outfiles <- paste0("../docs/", files)
outfiles <- sub("rmd$", "md", outfiles)

for (i in seq_along(files)) {
    render(files[i], output_file=outfiles[i])
}


