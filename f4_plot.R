library(tidyverse)
library(patchwork)
library(scales)
library(ggtext)
library(ggbreak) 
library(grid)

# Read in f4 output
f4dat = read.table("Twist_Mbuti.qpDstat_f4.R.out",
                   col.names=c("PopA", "PopB", "PopC", "PopD","F4", "StdErr", "Z", "ABBA","BABA", "SNPs")) 


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
  mutate(Bcountry = ifelse(str_detect(PopB, "Spain"), "Spain",
                           ifelse(str_detect(PopB, "Portugal"), "Portugal",
                                  ifelse(str_detect(PopB, "Mexico"), "Mexico",              
                                         ifelse(str_detect(PopB, "Indonesia"), "Indonesia",              
                                                "NA"))))) %>%
  mutate(CDcountry = ifelse(str_detect(PopC, "Spain"), "Spain",
                            ifelse(str_detect(PopC, "Portugal"), "Portugal",
                                   ifelse(str_detect(PopC, "Mexico"), "Mexico",              
                                          ifelse(str_detect(PopC, "Indonesia"), "Indonesia",              
                                                 "NA"))))) %>%
  mutate(pair_type = paste(BCtype,"_",Dtype))

#plot
plot <- ggplot(f4dat) +
  geom_vline(xintercept = 0) +
  geom_errorbar(aes(xmin=F4-2*StdErr, xmax=F4+2*StdErr, colour=factor(as.numeric(absZ)>2), y=pair_type), 
                linewidth=1, width=0) +
  geom_point(aes(x=F4, y=pair_type, fill=pair_type), size=4, shape=21) +
  scale_colour_manual(name= "|Z| score",values = c("black", "red"), labels = c("< 2", "> 2")) +
  geom_text(data=subset(f4dat, SNPs<10000), aes(label=paste(SNPs,"SNPs"), x=F4, y=pair_type),
            vjust = -1, size=3.1) +
  theme_bw() + 
  annotation_custom(grob = textGrob("Your Title", hjust = 0, vjust = 1, 
                                    gp = gpar(colour="red", fontsize = 14, fontface = "bold")),
                    xmin = Inf, xmax = Inf, ymin = Inf, ymax = Inf) +
  theme(legend.position = "bottom", text=element_text(size=12), axis.text.x = element_text(angle=90)) +
  facet_grid(CDcountry ~ Bcountry) +
  labs(title= "Population 1", x="f4 (Mbuti, Pop1.MethodA; Pop2.MethodA, Pop2.MethodB)", 
       fill= "Method A _ Method B", y = "Method A _ Method B")

plot 

ggsave("Twist_f4.pdf", width = 12, height = 10, units="in", dpi=600)
ggsave("Twist_f4.png", width = 12, height = 10, units="in", dpi=600)
