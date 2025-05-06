library(ggplot2)
library(ggstatsplot)
library(dplyr)
library(readr)

# get data from GitHub
urlfile="https://raw.githubusercontent.com/roberta-davidson/Davidson_etal_2024-Twist/main/index_hopping_stats2.tsv"
all_data<-read_delim(url(urlfile), delim="\t")

#index hopping 1 round of twist
ggbetweenstats(twist_1_df, Type, rate, type="r") 

#index hopping 1 round of twist
ggbetweenstats(twist_2_df, Type, rate, type="r")

all_data$Type[all_data$Type == "HE"] <- "4-library Pools"
all_data$Type[all_data$Type == "LE"] <- "2-library Pools"
all_data$Type[all_data$Type == "SG"] <- "All Other Pairs of Libraries (Capture & Shotgun)"
all_data$Type[all_data$Type == "TW"] <- "All Other Pairs of Libraries (Capture & Shotgun)"
all_data$Type[all_data$Type == "TW_SG"] <- "All Other Pairs of Libraries (Capture & Shotgun)"

all_data$Round_type <- paste0(all_data$Type," ", "TW", all_data$rounds )
all_data$Round_type[all_data$Round_type == "All Other Pairs of Libraries (Capture & Shotgun) TW1"] <- "All Other Pairs of Libraries (Capture & Shotgun)"
all_data$Round_type[all_data$Round_type == "All Other Pairs of Libraries (Capture & Shotgun) TW2"] <- "All Other Pairs of Libraries (Capture & Shotgun)"

png("/Users/shyamamama/Library/CloudStorage/Box-Box/projects/twist/index_hopping/index_hopping_1.png",
    width     = 8,
    height    = 6,
    units     = "in",
    res       = 1200)
ggbetweenstats(all_data, Type, rate, 
               xlab = "Index Hopping Group", 
               ylab = "Single-Index Hopping Rate", 
               type="r", 
               pairwise.display = "s", 
               results.subtitle=FALSE) +
  theme(axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank()) +
  scale_color_manual(values = c("#019dfe", "#001bcc", "#ff5500"))
dev.off()

png("/Users/shyamamama/Library/CloudStorage/Box-Box/projects/twist/index_hopping/index_hopping_2.png",
    width     = 12,
    height    = 6,
    units     = "in",
    res       = 1200)
ggbetweenstats(all_data, Round_type, rate, 
               xlab = "Index Hopping Group", 
               ylab = "Single-Index Hopping Rate", 
               type="r", 
               pairwise.display = "n",
               results.subtitle=FALSE) +
  theme(axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank()) +
  scale_color_manual(values = c("#00b6c1", "#007f92", "#01b277", "#00998b", "#ff5500"))
dev.off()
