## "kegg_species.txt" is downloaded from "http://rest.kegg.jp/list/organism"
species <- read.table("kegg_species.txt", sep = "\t", header = F, quote = "")
species <- species[, -1]
scientific_name <- gsub(" \\(.*", "", species[,2])
common_name <- gsub(".*\\(", "", species[,2])
common_name <- gsub("\\)", "", common_name)
kegg_species <- data.frame(kegg_code = species[, 1], 
                            scientific_name = scientific_name, 
                            common_name = common_name)
                            
save(kegg_species, file="kegg_species.rda", compress = "xz")                            


