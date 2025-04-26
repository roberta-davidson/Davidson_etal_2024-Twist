library(tidyverse)
library(patchwork)
library(scales)
library(ggtext)
library(ggbreak) 
library(grid)

# Read in f4 output
urlfile="https://raw.githubusercontent.com/roberta-davidson/Davidson_etal_2024-Twist/main/Table_S5.csv"
df<-read_delim(url(urlfile), delim="\t")


# add columns for plotting
f4dat = dplyr::mutate(f4dat, absZ=abs(Z)) %>%
  mutate(BCtype = ifelse(str_detect(PopB, "SG"), "SG",
                         ifelse(str_detect(PopB, "TW1"), "TW1",
                                ifelse(str_detect(PopB, "TW2"), "TW2",              
                                       "NA")))) %>%
  mutate(Dtype = ifelse(str_detect(PopD, "SG"), "SG",
                        ifelse(str_detect(PopD, "TW1"), "TW1",
                               ifelse(str_detect(PopD, "TW2"), "TW2",              
                                      "NA")))) %>%
  mutate(Bcountry = ifelse(str_detect(PopB, "EastIberia"), "EastIberia",
                           ifelse(str_detect(PopB, "WestIberia"), "WestIberia",
                                  ifelse(str_detect(PopB, "CentralAmerica"), "CentralAmerica",              
                                         ifelse(str_detect(PopB, "SouthEastAsia"), "SouthEastAsia",              
                                                "NA"))))) %>%
  mutate(CDcountry = ifelse(str_detect(PopC, "EastIberia"), "EastIberia",
                            ifelse(str_detect(PopC, "WestIberia"), "WestIberia",
                                   ifelse(str_detect(PopC, "CentralAmerica"), "CentralAmerica",              
                                          ifelse(str_detect(PopC, "SouthEastAsia"), "SouthEastAsia",              
                                                 "NA"))))) %>%
  mutate(pair_type = paste(BCtype,"_",Dtype))

textbox_text <- "Population 2"

#plot
plot <- ggplot(f4dat) +
  geom_vline(xintercept = 0, colour="grey40", linewidth=0.5) +
  geom_point(aes(x=F4, y=pair_type), shape=1, size=3) +
  geom_errorbar(aes(xmin=F4-3*StdErr, xmax=F4+3*StdErr, colour=factor(as.numeric(absZ)>2), y=pair_type), linewidth=0.5, width=0) +
  geom_errorbar(aes(xmin=F4-2*StdErr, xmax=F4+2*StdErr, colour=factor(as.numeric(absZ)>2), y=pair_type), linewidth=1, width=0) +
  scale_colour_manual(name= "|Z| score",values = c("blue4", "red"), labels = c("(0,2)", "[2,3)")) +
  scale_shape_manual(values=c(0,1,2,4,5,6)) +
  geom_text(data=subset(f4dat, SNPs<10000), aes(label=paste(SNPs), x=F4, y=pair_type), vjust = -1, size=3.1) +
  theme_bw() + 
  annotation_custom(grob = textGrob("Your Title", hjust = 0, vjust = 1, gp = gpar(colour="red", fontsize = 14, fontface = "bold")), xmin = Inf, xmax = Inf, ymin = Inf, ymax = Inf) +
  theme(legend.position = "bottom", text=element_text(size=12), axis.text.x = element_text(angle=90), plot.margin = margin(0.5, 0.2, 0.5, 0.5, "cm")) +
  facet_grid(CDcountry ~ Bcountry) +
  labs(title= "Population 1", x="f4 (Mbuti, Pop1.MethodA; Pop2.MethodA, Pop2.MethodB)", 
       fill= "Method A _ Method B", y = "Method A _ Method B") 
plot 

textbox <- textGrob("Population 2", gp = gpar(col = "black", fontsize = 14), vjust=1, rot = -90)
combined_plot <- arrangeGrob(plot, right = textbox, widths = c(1, 0.0005))
combined_plot
grid.newpage()
grid.draw(combined_plot)
plot <- as.ggplot(combined_plot)
ggsave("Twist_f4.pdf", width = 12, height = 10, units="in", dpi=600)
ggsave("Twist_f4.png", width = 12, height = 10, units="in", dpi=600)
