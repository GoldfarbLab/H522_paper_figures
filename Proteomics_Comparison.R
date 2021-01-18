library(here)
source(here("common.R"))

################################################################################
# Read data
################################################################################
H522.data <- read_tsv(here("data_processed/proteinsFCvsMockAveraged.txt"))
caco2.data <- read_csv(here("data/Figure 6/MS/Bojkova_Protein_Caco2.csv"))
A549.data <- read_csv(here("data/Figure 6/MS/Mann_proteinDDA_A549ACE2.csv"))
################################################################################

# select columns of interest, extract leading protein and remove isoform suffix
A549.data <- A549.data %>% 
  separate(`protein AC`, c("Leading razor protein"), sep=";", remove=F, extra="drop") %>%
  separate("Leading razor protein", c("Leading canonical"), sep="-", remove=F, extra="drop") %>%
  select("Leading canonical", "A549" = "fold change vs mock (log2)")

caco2.data <- caco2.data %>% 
  separate(`UniProt Accession`, c("Leading canonical"), sep=";", remove=F, extra="drop") %>%
  select("Leading canonical", "Caco2" = "Ratio 24h")

H522.data <- H522.data %>%
  filter(Condition == "24 hr") %>% 
  select(`Leading canonical`, "H522" = FC)

# align the data
data <- inner_join(A549.data, caco2.data, by="Leading canonical") %>%
  inner_join(H522.data, by="Leading canonical")

p <- (ggplot(data=data, aes(x=H522, y=A549))
  + geom_point()
)

print(p)
