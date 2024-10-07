### Taxonomic Analysis Code
## Load the necessary packages
library(reshape2)
library(readr)
library(plyr)
library(dplyr)
library(tidyr)
## Read txt file
OTU <- read.table(file = 'Kraken_species_count_table.cast.txt', sep = '\t', header = TRUE)

## Calculate the count of each taxon in each sample and only remove taxa not identified
OTU <- OTU %>%
  gather(Species, Count, 2:(ncol(OTU))) %>%
  filter(Count > 0)

## Tabulate read counts by genus
OTU <- OTU %>% group_by(Sample, Species) %>% summarise(Count = sum(Count))

## Convert count to relative abundance and add column (require dplyr)
OTU <- ddply(OTU, .(Sample), mutate, RelativeAbundance = (Count * 100) / sum(Count))

## Save the processed table into txt file
write_tsv(OTU,"kraken_species_abundance.tidy.txt")

## Load the necessary packages for plot
library(wilkoxmisc)
library(ggplot2)

## Read the previously saved txt file
OTU <- read_tsv("kraken_species_abundance.tidy.txt")

## Collapse taxa table to only 10 top phyla, genus, family, etc.
OTUTable <- collapse_taxon_table(OTU, n = 10, Rank = "Species") ## Read the metadata table
Meta <- read_tsv("metadata_skin.tsv")

## Merge relative abundance table and metatable together
OTUTable <- merge(OTUTable, Meta, by = "Sample", all.x = TRUE) ## Save the table in txt file
write_tsv(OTUTable, "Top10Species.txt")

## Plotting the Top 10 Most Abundant Species
OTUTable <- read_tsv("Top10Species.txt")

##Get the unique column names
species <- unique(OTUTable$Species) list(species)
sample <- unique(OTUTable$Sample)

## Set factor for two variables
OTUTable$Species <- factor(OTUTable$Species,
levels = c("Brevundimonas_diminuta",
"Chryseobacterium_taklimakanense", "Cutibacterium_acnes", "Dermacoccus_nishinomiyaensis",
"Gordonia_bronchialis","Janibacter_indicus", "Micrococcus_luteus", "Moraxella_osloensis",
"Xanthomonas_campestris", "Minor/Unclassified")) OTUTable$Sample <- factor(OTUTable$Sample, levels = sample)

## Plot the relative abundance
Plot1 <- ggplot(OTUTable, aes(x = Sample, y = RelativeAbundance, fill = Species)) +
geom_bar(stat="identity") + scale_fill_brewer(palette = "Paired") + theme_classic() +
facet_wrap(~Individual, scales = "free_x") + scale_y_continuous(expand = c(0, 0)) +
ylab("Relative abundance (%)") + xlab("Sample") + theme(axis.text.x = element_text(angle = 90)) +
ggtitle("Ten Most Abundance") + theme(plot.title=element_text(hjust=0.5, size = 18, face = "bold"))

## Save the plot in pdf file
ggsave("top10species_by_individual.pdf")

################################################################
## Prepare data for Indicspecies
## Load the pacakges
library(dplyr)
library(tidyr)
library(indicspecies)

## Read the previously saved txt file
tidy <- read_tsv("kraken_species_abundance.tidy.txt")

## Collapse taxa table to only 30 top phyla, genus, family, etc.
top30 <- collapse_taxon_table(tidy, n = 30, Rank = "Species") ## Read the metadata table
Meta <- read_tsv("metadata_skin.tsv")

## Merge relative abundance table and metatable together
top30 <- merge(top30, Meta, by = "Sample", all.x = TRUE) ## Save the table in txt file
write_tsv(top30, "Top30Species.txt") ## Read the top 30 txt file
top30 <- read_tsv("Top30Species.txt")

#Prepare abundance matrix
abundance_matrix <- top30 %>%
select(Sample, Species, RelativeAbundance) %>% spread(Species, RelativeAbundance, fill = 0)
# Ensure row names are sample IDs
abundance_matrix <- as.data.frame(abundance_matrix) # Convert table to data frame
row.names(abundance_matrix) <- abundance_matrix$Sample abundance_matrix <- abundance_matrix[,-1] # Remove the Sample column

## Indicspecies
# Prepare the grouping vector
grouping_vector <- top30 %>% select(Sample, Individual) %>% distinct() %>%
arrange(Sample)

# Check for distinct individuals
if (length(unique(grouping_vector$Individual)) < 2) { stop("The grouping vector must contain at least two distinct
groups.") }

# Ensure the order of the grouping vector matches the abundance matrix
grouping_vector <- grouping_vector$Individual[match(row.names(abundance_matrix), grouping_vector$Sample)]

# Run the multipatt function
set.seed(123) # For reproducibility 
indicator_species_results <- multipatt(abundance_matrix, grouping_vector, func = "IndVal.g", duleg = TRUE)
# View the summary of the results
summary(indicator_species_results)
# Extract significant results
significant_results <- indicator_species_results$sign significant_results$Species <- row.names(significant_results) significant_results <- significant_results[, c("Species", names(significant_results)[names(significant_results) != "Species"])]
# Reshape `significant_results` to long format

# Initialize an empty data frame for the result
result_df <- data.frame(subject=character(), Species=character(), stat=numeric(), p.value=numeric(), stringsAsFactors=FALSE)

# Loop through each subject and get the top 3 species based on stat
for (subject in c("s.Subject_1", "s.Subject_2", "s.Subject_3", "s.Subject_4")) {
subject_data <- significant_results[significant_results[[subject]] == 1, ]
top_species <- subject_data[order(-subject_data$stat), ][1:3, ] top_species <- top_species[, c("Species", "stat", "p.value")] top_species$subject <- subject
result_df <- rbind(result_df, top_species)
}

# Reorder columns
result_df <- result_df[, c("subject", "Species", "stat", "p.value")]

### Alpha Diversity
## Load the necessary packages
library(readr)
library(reshape2)
library(dplyr)
library(GUniFrac)
library(vegan)
library(stringr)
library(tidyr)
library(wilkoxmisc)
library(fossil)
library(OTUtable)
library(ggplot2)
## Read the table into R
table <- read.table(file = 'Kraken_species_count_table.cast.txt', sep = '\t', header = TRUE)

## Calculate the count of each taxon in each sample and only remove taxa not identified
table <- table %>%
  gather(Species, Count, 2:(ncol(table))) %>%
  filter(Count > 0)

## Cast the table into dataframe and rename their column names
table <- dcast(table, Sample~Species,value.var = "Count", fill = 0) colnames(table) <- c("Sample", paste0("Species",1:4792))

## Save the merged
write_tsv(table,"kraken_report_all_species.cast.txt")
table <- read.csv("kraken_report_all_species.cast.txt", header=TRUE,
sep = "\t", row.names = 1, as.is=TRUE) #rarefy to the minimal depth
table.rarefy <- Rarefy(table,depth = min(rowSums(table)))$otu.tab.rff
table.rarefy <- as.data.frame(as.table(table.rarefy))
names(table.rarefy)[1:2] <- c("Sample","Species") 
names(table.rarefy)[3] <- c("Count")

#for richness
table.rarefy <- table.rarefy %>% filter(!Count== 0) write_tsv(table.rarefy,"kraken_report_all_species.rarefied.tidy.tsv")

# for shannon and evenness
table.rarefy.cast <- dcast(table.rarefy, Species ~Sample, value.var = "Count", fill = 0) write_tsv(table.rarefy.cast,"kraken_report_all_species.rarefied.cast. tsv")

## Calculate richness for rarefied table and save the result
richness <- read_tsv("kraken_report_all_species.rarefied.tidy.tsv") %>%
    select(Sample) %>%
    mutate(richness = 1) %>% 
    group_by(Sample) %>% 
    mutate(richness= sum(richness)) %>% 
    ungroup() %>%
    unique() 
write_tsv(richness,"skin_microbial_richness.txt")

## Calculate shannon diversity for rarefied table and save the result
alpha <- read.table("kraken_report_all_species.rarefied.cast.tsv", header=T, sep="\t", fill =TRUE, quote="", row.names = NULL) alpha[1:6,1:6]
Sample <- colnames(alpha[2:ncol(alpha)])
Species <- alpha[1:nrow(alpha),1] names(alpha) <- NULL
alpha[,1] <- NULL

## Transpose the data frame
alpha_t <- data.frame(t(alpha)) 
alpha_t[is.na(alpha_t)] <- 0 
colnames(alpha_t) <- Species

## Calculate the shannon diversity for global samples
Shannon <- diversity(alpha_t, "shannon")

## Merge the sample info with alpha diversity
data <- data.frame(cbind(Sample,Shannon)) write_tsv(data, "skin_microbial_shannon.txt")

## Calculate pielou's evenness and save the result
table <- read_tsv("kraken_report_all_species.rarefied.cast.tsv") table[,1] <- NULL
evenness <- apply(table,2,pielou)
evenness <- as.data.frame(evenness)
evenness$Sample <- rownames(evenness) write_tsv(evenness,"skin_microbial_pielou_evenness.txt")

## Merge richness, shannon index and evenness as the final alpha diversity result
richness <- read_tsv("skin_microbial_richness.txt")
shannon <- read_tsv("skin_microbial_shannon.txt")
merge <- left_join(richness,shannon)
evenness <- read_tsv("skin_microbial_pielou_evenness.txt") 
merge <- left_join(merge, evenness) 
write_tsv(merge,"skin_microbial_alpha_diversity_rarefied.txt")

## Load the necessary packages for barplot
library(readr)
library(dplyr)
library(reshape2)
library(tidyr)
library(car)
library(viridis)
## Read the alpha diversity table and reformat it
table <- read_tsv("skin_microbial_alpha_diversity_rarefied.txt") 
table <- table %>%
  gather(alpha, value, 2:(ncol(table))) %>%
  filter(value > 0)

## Read metadata information and merge alpha diversity result with meta information
meta <- read_tsv("metadata_skin.tsv") 
table <- left_join(table,meta) View(table)

#Sort data based on Day
day_component <- as.numeric(sub("Day ([0-9]+) .*", "\\1", table$Sampling_day))
sorted_indices <- order(day_component)
sorted_table <- table[sorted_indices, ]

## Factor the day variable in the merged table
level <- unique(sorted_table$Sampling_day)
sorted_table$Sampling_day <- factor(sorted_table$Sampling_day, levels = level)

## Plot the richness and save it in pdf file
install.packages("RColorBrewer") library(RColorBrewer)
colourCount = length(unique(table$Sampling_day)) getPalette = colorRampPalette(brewer.pal(9, "Set1"))
Plot <- ggplot(sorted_table, aes(x= Individual, y=value)) + 
    geom_boxplot(outlier.size = 0.5, outlier.shape = NA) + 
    geom_jitter(size = 3, aes(color = Sampling_day)) + theme_classic() + 
    facet_wrap(~alpha, scales = "free") +
    scale_color_manual(values = getPalette(colourCount)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("Value") +
    xlab("Individual") +
    ggtitle("Alpha Diversity") + theme(plot.title=element_text(hjust=0.5, size = 18, face = "bold"))
    ggsave("alpha_diversity_by_individual.pdf")

## Statistical tests
sorted_table$Day_numeric <- as.numeric(sorted_table$Sampling_day)
subtable <- sorted_table %>% filter(alpha == "richness")

#day as a fixed numeric variable
lm <- lm(subtable$value~Individual+Sample_type+Day_numeric, data=subtable)

#perform type I anova
anova(lm)