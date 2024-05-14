library(tidyverse)
library(patchwork)
library(scales)
library(ggtext)
library(ggbreak) 
library(grid)
library(plotly)
library(cowplot)
library(ggplotify)
library(gtable)
library(readr)
library(gridExtra)

#### get main data frame and set up ####
urlfile="https://raw.githubusercontent.com/roberta-davidson/Davidson_etal_2024-Twist/main/Table_S1.csv"
df<-read_delim(url(urlfile), delim="\t")
colnames(df) <- gsub(" ", "_", colnames(df))
colnames(df) <- gsub("-", "_", colnames(df))
colnames(df) <- gsub("#", "", colnames(df))
colnames(df) <- gsub("_\\(%\\)", "", colnames(df))

#### PLOTTING ####
colours <- c("SG"="#ff0000","TW1"="#01cc54", "TW2"="#066601","1"="#fe7a00", "2"="#019dfe", "4"="#001bcc")
shapes <- c("SG"=21,"TW1"=22, "TW2"=24,"1"=21, "2"=21, "4"=21)

#### Does number of libraries in pool affect SNPs obtained? ####
## Figure S2

plotb2 <- ggplot(subset(df, method!="SG"), aes(x=MappableBioEndoDNA_SG,y=SNPsPer1MSequencedReads, fill=as.factor(n_pool),colour=as.factor(n_pool))) +
  geom_point(size=5, aes(shape=Deidentified_IDs)) +
  scale_shape_manual(guide="none",values=as.character(letters)) +
  theme_bw() +
  scale_y_continuous(labels = scales::comma_format(), limits=c(0,50000)) +
  scale_fill_manual(values=c("1"="#fe7a00", "2"="#019dfe", "4"="#001bcc")) +
  scale_colour_manual(values=c("1"="#fe7a00", "2"="#019dfe", "4"="#001bcc")) +
  scale_x_continuous(limit=c(0,1.5)) +
  theme(text=element_text(size=14), legend.position=c(0.75,0.85)) +
  labs(x="mappable endo%, screening", y="", fill="# libs in pool",colour="# libs in pool", title="") +
  facet_grid(method ~ .)
plotb2

plotb3 <- ggplot(subset(df, method!="SG"), aes(x=MappableBioEndoDNA_SG,y=SNPsPer1MSequencedReads, fill=as.factor(n_pool),colour=as.factor(n_pool))) +
  geom_point(size=5, aes(shape=Deidentified_IDs)) +
  scale_shape_manual(guide="none",values=as.character(letters)) +
  theme_bw() +
  scale_y_continuous(labels = scales::comma_format()) +
  scale_fill_manual(values=c("1"="#fe7a00", "2"="#019dfe", "4"="#001bcc")) +
  scale_colour_manual(values=c("1"="#fe7a00", "2"="#019dfe", "4"="#001bcc")) +
  scale_x_continuous(limit=c(1.5,45)) +
  theme(text=element_text(size=14), legend.position="none") +
  labs(x="mappable endo%, screening", y="SNPs per million sequenced read pairs", fill="# libs in pool",colour="# libs in pool", title="") +
  facet_grid(method ~ .)
plotb3

plot <- plotb3+plotb2 + plot_layout(axes="collect", axis_titles = "collect")
plot

ggsave("Twist_pooling_SNPs.pdf", width = 8, height = 6, units="in", dpi=600)

#### pooling equity ####
### Figure 6
df2 <- select(df, method, De_identified_lib_ID, Deidentified_IDs, Nr._Sequenced_Read_Pairs, Nr._Reads_Into_Mapping, Nr._Mapped_Reads, "Nr._Mapped_Reads_Post_Filter", Nr._Dedup._Mapped_Reads, pool_name) %>%
  filter(str_detect(pool_name,"E")) %>% 
  pivot_longer(cols=c(Nr._Sequenced_Read_Pairs, Nr._Reads_Into_Mapping, Nr._Mapped_Reads, Nr._Mapped_Reads_Post_Filter, Nr._Dedup._Mapped_Reads), names_to = "name", values_to="read_count")

df2$name <- factor(df2$name, levels = c("Nr._Sequenced_Read_Pairs","Nr._Reads_Into_Mapping", "Nr._Mapped_Reads", "Nr._Mapped_Reads_Post_Filter", "Nr._Dedup._Mapped_Reads"))
df2$pool_name <- factor(df2$pool_name, levels = c("LE1", "LE2", "LE3","HE4","HE5","HE6"))

labs <- c("LE1"="LE1: 0.3-0.4 %","LE2"="LE2: 0.5-0.7 %","LE3"="LE3: 0.9-1.4 %","HE4"="HE4: 2.5-3.4 %","HE5"="HE5: 11.0-13.8 %","HE6"="HE6: 38.4-44.1 %")
labs2 <- c("Nr._Sequenced_Read_Pairs"="Sequenced Read Pairs", "Nr._Reads_Into_Mapping"="Reads Into Mapping", "Nr._Mapped_Reads"="Mapped Reads", "Nr._Mapped_Reads_Post_Filter"="Reads MAQ 25", "Nr._Dedup._Mapped_Reads"="Dedup Mapped Reads")

#plots
plot <- ggplot(df2, aes(x=method,y=read_count, fill=Deidentified_IDs)) +
  geom_col(position="fill", colour="black") +
  theme_bw() +
  theme(legend.position = "none", text=element_text(size=12)) +
  facet_grid(pool_name ~ factor(name, levels = c("Nr._Sequenced_Read_Pairs","Nr._Reads_Into_Mapping", "Nr._Mapped_Reads", "Nr._Mapped_Reads_Post_Filter", "Nr._Dedup._Mapped_Reads")), labeller = labeller(pool_name = labs, name = labs2)) + 
  labs(title="Read count", y="Proportion of read count", x="Method")
plot

textbox <- textGrob("Library pool (with range of mappable endo%)", gp = gpar(col = "black", fontsize = 14), vjust=0.8, rot = -90)
combined_plot <- arrangeGrob(plot, right = textbox, widths = c(1, 0.001))
combined_plot
grid.newpage()
grid.draw(combined_plot)
plot <- as.ggplot(combined_plot)
plot

ggsave(plot=plot, "Twist_pooling_equity.pdf", width = 12, height = 10, units="in", dpi=600)

#### capture efficacy #### 
plot1 <- ggplot(df, aes(x=MappableBioEndoDNA_SG, y=as.numeric(SNPs_Covered), fill=method, colour=method )) +
  geom_point(size=5, alpha=1, aes(shape=Deidentified_IDs)) +
  geom_smooth(data=subset(df, method=="SG"), formula = y~x, method="lm", alpha=0.3) +
  geom_smooth(data=subset(df, method!="SG"), formula = y~log(x), method="lm", alpha=0.3) +
  scale_colour_manual(values=colours) + scale_fill_manual(values=colours) + 
  scale_shape_manual(values=as.character(letters)) +
  scale_y_continuous(limits=c(0,1000000), labels=comma_format()) +
   theme_bw() +
  theme(legend.position = "none", text=element_text(size=12)) +
  labs(title="A",x="mappable endo%, screening SG", y="SNPs on Twist aDNA panel")
plot1

plot2 <- ggplot(df, aes(x=MappableBioEndoDNA_SG, y=SNPsPer1MSequencedReads,fill=method, colour=method )) +
  geom_point(size=5, alpha=1, aes(shape=Deidentified_IDs)) +
  geom_smooth(data=subset(df, method=="SG"), formula = y~x, method="lm", alpha=0.3) +
  geom_smooth(data=subset(df, method!="SG"), formula = y~log(x), method="lm", alpha=0.3) +
  scale_colour_manual(values=colours) + scale_fill_manual(values=colours) + 
  scale_shape_manual(values=as.character(letters)) +
  scale_y_continuous(limits=c(0,110000), labels=comma_format()) +
  theme_bw() +
  theme(legend.position = "none", text=element_text(size=12)) +
  labs(title="B",x="mappable endo%, screening SG", y="SNPs per million sequenced read pairs")
plot2

plot3 <- ggplot(df, aes(x=MappableBioEndoDNA_SG, y=SNPsPer1MReadsToMapping, fill=method, colour=method )) +
  geom_point(size=5, alpha=1, aes(shape=Deidentified_IDs)) +
  geom_smooth(data=subset(df, method=="SG"), formula = y~x, method="lm", alpha=0.3) +
  geom_smooth(data=subset(df, method!="SG"), formula = y~log(x), method="lm", alpha=0.3) +
  scale_colour_manual(values=colours) + scale_fill_manual(values=colours) + 
  scale_shape_manual(values=as.character(letters)) +
  scale_y_continuous(limits=c(0,115000), labels=comma_format()) +
  theme_bw() +
  theme(legend.position = "none", text=element_text(size=12)) +
  labs(title="C",x="mappable endo%, screening SG", y="SNPs per million reads into mapping")
plot3

plot4 <- ggplot(df, aes(x=MappableBioEndoDNA_SG, y=SNPsPer1MMappedReads, fill=method, colour=method )) +
  geom_point(size=5, alpha=1, aes(shape=Deidentified_IDs)) +
  geom_smooth(formula = y~x, method="lm", alpha=0.3) +
  scale_colour_manual(values=colours) + scale_fill_manual(values=colours) + 
  scale_shape_manual(values=as.character(letters)) +
  scale_y_continuous(limits=c(0,200000), labels=comma_format()) +
  theme_bw() +
  theme(legend.position = "none",text=element_text(size=12)) +
  labs(title="D", x="mappable endo%, screening SG", y="SNPs per million mapped reads")
plot4

plot5 <- ggplot(df, aes(x=MappableBioEndoDNA_SG, y=SNPsPer1MMappedFilteredReads,fill=method, colour=method )) +
  geom_point(size=5, alpha=1, aes(shape=Deidentified_IDs)) +
  geom_smooth(formula = y~x, method="lm", alpha=0.3) +
  scale_colour_manual(values=colours) + scale_fill_manual(values=colours) + 
  scale_shape_manual(values=as.character(letters)) +
  scale_y_continuous(limits=c(0,300000), labels=comma_format()) +
  theme_bw() +
  theme(legend.position = "none", text=element_text(size=12)) +
  labs(title="E", x="mappable endo%, screening SG", y="SNPs per million filtered reads")
plot5

plot6 <- ggplot(df, aes(x=MappableBioEndoDNA_SG, y=SNPsPer1MUniqueReads,fill=method, colour=method )) +
  geom_point(size=5, alpha=1, aes(shape=Deidentified_IDs)) +
  geom_smooth(formula = y~x, method="lm", alpha=0.3) +
  scale_colour_manual(values=colours) + scale_fill_manual(values=colours) + 
  scale_shape_manual(guide = "none",values=as.character(letters)) +
  scale_y_continuous(limits=c(0,500000), labels=comma_format()) +
  theme_bw() +
  theme(text=element_text(size=12)) +
  labs(title="F",x="mappable endo%, screening SG", y="SNPs per million dedup reads")
plot6

plot7 <- ggplot(df, aes(x=MappableBioEndoDNA_SG, y=SequencedBioEndoDNA, fill=method, colour=method )) +
  geom_point(size=5, alpha=1, aes(shape=Deidentified_IDs)) +
  geom_smooth(data=subset(df, method=="SG"), formula = y~x, method="lm", alpha=0.3) +
  geom_smooth(data=subset(df, method!="SG"), formula = y~log(x), method="lm", alpha=0.3) +
  scale_colour_manual(values=colours) + scale_fill_manual(values=colours) + 
  scale_shape_manual(values=as.character(letters)) +
  geom_abline(slope=1) +
  scale_y_continuous(limits=c(0,100)) +
  theme_bw() +
  theme(legend.position = "none", text=element_text(size=12)) +
  labs(x="mappable endo%, screening SG", y="sequenced endo%",title="G")
plot7

plot8 <- ggplot(df, aes(x=MappableBioEndoDNA_SG, y=MappableBioEndoDNA, fill=method, colour=method )) +
  geom_point(size=5, alpha=1, aes(shape=Deidentified_IDs)) +
  geom_smooth(data=subset(df, method=="SG"), formula = y~x, method="lm", alpha=0.3) +
  geom_smooth(data=subset(df, method!="SG"), formula = y~log(x), method="lm", alpha=0.3) +
  scale_colour_manual(values=colours) + scale_fill_manual(values=colours) + 
  scale_shape_manual(values=as.character(letters)) +
  geom_abline(slope=1) +
  scale_y_continuous(limits=c(0,100)) +
  theme_bw() +
  theme(legend.position = "none", text=element_text(size=12)) +
  labs(x="mappable endo%, screening SG", y="mappable endo%",title="H")
plot8

plot9 <- ggplot(df, aes(x=MappableBioEndoDNA_SG, y=FilteredBioEndoDNA,fill=method, colour=method )) +
  geom_point(size=5, alpha=1, aes(shape=Deidentified_IDs)) +
  geom_smooth(data=subset(df, method=="SG"), formula = y~x, method="lm", alpha=0.3) +
  geom_smooth(data=subset(df, method!="SG"), formula = y~log(x), method="lm", alpha=0.3) +
  scale_colour_manual(values=colours) + scale_fill_manual(values=colours) + 
  scale_shape_manual(values=as.character(letters)) +
  geom_abline(slope=1) +
  scale_y_continuous(limits=c(0,100)) +
  theme_bw() +
  theme(legend.position = "none", text=element_text(size=12)) +
  labs(x="mappable endo%, screening SG", y="filtered endo%",title="I")
plot9

plot10 <- ggplot(df, aes(x=MappableBioEndoDNA_SG, y=UniqueBioEndoDNA,fill=method, colour=method )) +
  geom_point(size=5, alpha=1, aes(shape=Deidentified_IDs)) +
  geom_smooth(data=subset(df, method=="SG"), formula = y~x, method="lm", alpha=0.3) +
  geom_smooth(data=subset(df, method!="SG"), formula = y~log(x), method="lm", alpha=0.3) +
  scale_colour_manual(values=colours) + scale_fill_manual(values=colours) + 
  scale_shape_manual(values=as.character(letters)) +
  #scale_shape_manual(values=shapes) +
  geom_abline(slope=1) +
  scale_y_continuous(limits=c(0,50)) +
  # scale_x_continuous(limits=c(0,4)) +
  theme_bw() +
  theme(legend.position="none",text=element_text(size=12)) +
  labs(x="mappable endo%, screening SG",y="unique endo%", title="J")
plot10

layout <- "
#AABBCC#
#DDEEFF#
GGHHIIJJ
"
plt <- plot1 + plot2 + plot3 + plot4 + plot5 + plot6 + plot7 + plot8 + plot9 + plot10 +
  plot_layout(design = layout, axes="collect", axis_titles = "collect")
plt

ggsave("Twist_endo_SNPs.pdf", width = 14, height = 11.5, units="in", dpi=600)

#### cost effectiveness ####
### Figure 3
#cost
plot9a <- ggplot(df, aes(x=MappableBioEndoDNA_SG, y=CostPerSNP, fill=method, colour=method )) +
  geom_smooth(formula = y~log(x), method="lm", alpha=0.2) +
  geom_point(size=5,alpha=1, aes(shape=Deidentified_IDs)) +
  scale_colour_manual(values=colours) + scale_fill_manual(values=colours) + 
  scale_shape_manual(values=as.character(letters), guide="none") +
  scale_y_log10(labels=comma_format(), n.breaks=12) +
   theme_bw() +
  theme(legend.position=c(0.8,0.8), text=element_text(size=12)) +
  labs(x="mappable endo%, screening SG", y="cost per SNP",title="A")
plot9a

#prepare data frame 
df_2 <- select(df, Deidentified_IDs, CostPerSNP, method, MappableBioEndoDNA_SG) 
df_2 <- pivot_wider(df_2, names_from = "method", values_from="CostPerSNP", names_prefix="costperSNP_")
df_2 <- pivot_longer(df_2, values_to="cost_TW", names_to="method", names_prefix="costperSNP_", cols=c("costperSNP_TW1", "costperSNP_TW2")) %>%
  mutate(CostperSNP_TW = cost_TW/costperSNP_SG) %>%
  mutate(CostperSNP_TW_INV = 1/(cost_TW/costperSNP_SG))

#calculate lm
fit_TW1 <- lm(CostperSNP_TW_INV ~ (MappableBioEndoDNA_SG), data = subset(df_2, method=="TW1"))
rsquared_TW_1 <- summary(fit_TW1)$r.squared
fit_TW2 <- lm(CostperSNP_TW_INV ~ log(MappableBioEndoDNA_SG), data = subset(df_2, method=="TW2"))
rsquared_TW_2 <- summary(fit_TW2)$r.squared

#plot
plot9b <- ggplot(df_2, aes(x=MappableBioEndoDNA_SG, y=CostperSNP_TW_INV, fill=method, colour=method )) +
  geom_smooth(data=subset(df_2, method=="TW2"), formula = y~log(x), method='lm',alpha=0.2) +
  geom_smooth(data=subset(df_2, method=="TW1"), formula = y~x, method="lm", alpha=0.2) +
  geom_point(size=5,alpha=1, aes(shape=Deidentified_IDs)) +
  scale_colour_manual(values=colours) + scale_fill_manual(values=colours) + 
  scale_shape_manual(values=as.character(letters), guide="none") +
  annotate("text", x = 35, y = 80, label = paste("r² TW2 =", round(rsquared_TW_2, 3)), vjust = 0.8, hjust = 0.8, size=5,colour="#066601") +
  annotate("text", x = 35, y = 90, label = paste("r² TW1 =", round(rsquared_TW_1, 3)), vjust = 0.8, hjust = 0.8, size=5,colour="#01cc54") +
  theme_bw() +
  theme(text=element_text(size=12), legend.position="none") +
  labs(x="mappable endo%, screening SG", y="cost saving compared to SG",title="B")
plot9b

layout <- "
AA
BB
"

plt <- plot9a / plot9b + plot_layout(design = layout, axes="collect", axis_titles = "collect")
plt


ggsave("Twist_cost_snp.pdf", width = 8, height = 8, units="in", dpi=600)

### Figure 4
plotA <- ggplot(df, aes(x=as.numeric(SNPs_Covered), y=CostPerSNP)) +
  geom_line(aes(group=Deidentified_IDs, colour=MappableBioEndoDNA_SG), linewidth=1) +
  geom_point(aes(fill=method, colour=MappableBioEndoDNA_SG, shape=method), size=4, stroke=1) +
  scale_fill_manual(values=c("grey90","grey50","black")) +
  scale_shape_manual(values=shapes) +
  scale_color_gradientn(colours = rainbow(3), trans="log10", n.breaks=12) +
  scale_y_log10(labels=comma_format(), n.breaks=20) +
  scale_x_log10(labels=comma_format(), n.breaks=12) + 
  guides(color = guide_colorbar(barheight = 15)) +
  theme_bw() +
  theme(legend.position = "none", text=element_text(size=12), axis.text.x = element_text(angle=90)) +
  labs(title="A", x="SNPs covered on Twist aDNA Panel", y="average AUD per SNP")
plotA

plotB <- ggplot(df, aes(x=SNPsPer1MSequencedReads, y=CostPerSNP)) +
  geom_line(aes(group=Deidentified_IDs, colour=MappableBioEndoDNA_SG), linewidth=1) +
  geom_point(aes(fill=method, colour=MappableBioEndoDNA_SG, shape=method), size=4, stroke=1) +
  scale_fill_manual(values=c("grey90","grey50","black")) +
  scale_shape_manual(values=shapes) +
  scale_color_gradientn(colours = rainbow(3), trans="log10", n.breaks=12) +
  scale_y_log10(labels=comma_format(), n.breaks=20) +
  scale_x_log10(labels=comma_format(), n.breaks=12) +
  guides(color = guide_colorbar(barheight = 15)) +
  theme_bw() +
  theme(text=element_text(size=12), axis.text.x = element_text(angle=90),
        # legend.title = element_text(angle = 0, hjust = 0.5, vjust = 1.5, margin = margin(b = 5))
  ) +
  labs(title="B", x="SNPs per million sequenced read pairs", colour="mappable endo%")
plotB

plotB_f <- plotB + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

patch <- plotA + plotB_f
patch

ggsave("Twist_cost_sample_4.pdf", width = 14, height = 7, units="in", dpi=600)

#### shotgun vs SHOTGUN ####
### Figure S1
df <- df %>%
  filter(method=="SG")

#linear models
#a
lm_model_a <- lm(MappableBioEndoDNA ~ MappableBioEndoDNA_SG, data = df)
slope_a <- coef(lm_model_a)[2]  # Extract the slope
rsquared_a <- summary(lm_model_a)$r.squared  # Extract the R-squared value
#b
lm_model_b <- lm(UniqueBioEndoDNA ~ UniqueBioEndoDNA_SG, data = df)
slope_b <- coef(lm_model_b)[2]  # Extract the slope
rsquared_b <- summary(lm_model_b)$r.squared  # Extract the R-squared value
#c
lm_model_c <- lm(SNPs_Covered ~ SNPs_Covered_screen, data = df)
slope_c <- coef(lm_model_c)[2]  # Extract the slope
rsquared_c <- summary(lm_model_c)$r.squared  # Extract the R-squared value
#d
lm_model_d <- lm(SNPsPer1MReadsToMapping ~ SNPsPer1MReadsToMapping_SG, data = df)
slope_d <- coef(lm_model_d)[2]  # Extract the slope
rsquared_d <- summary(lm_model_d)$r.squared  # Extract the R-squared value

#plots
plota <- ggplot(df, aes(x=MappableBioEndoDNA_SG, y = MappableBioEndoDNA, group=method)) +
  geom_point(size=4, colour="black", shape=21, aes(fill=MappableBioEndoDNA_SG)) +
  geom_smooth(method="lm", alpha=0.2, formula=y~x) +
  geom_abline(slope=1) +
  scale_fill_gradientn(colours = rainbow(3), trans="log10", n.breaks=12) +
  theme_light() +
  annotate("text",  x = 0.1, y = (max(df$MappableBioEndoDNA, na.rm = T))*0.8, label = paste("Slope: ", round(slope_a, 2), "\nR-squared: ", round(rsquared_a, 2)),hjust = 0, vjust = 1,size = 6) +
  theme(legend.position = "none", text=element_text(size=12), axis.text.x = element_text(angle=90, hjust=1)) +
  labs(title="A. mappable biological endogenous DNA %", x="screening shotgun",y="deep shotgun")
plota

plotb <- ggplot(df, aes(x=UniqueBioEndoDNA_SG, y = UniqueBioEndoDNA, group=method)) +
  geom_point(size=4, colour="black", shape=21, aes(fill=MappableBioEndoDNA_SG)) +
  geom_smooth(method="lm", alpha=0.2, formula=y~x) +
  geom_abline(slope=1) +
  scale_fill_gradientn(colours = rainbow(3), trans="log10", n.breaks=12) +
  theme_light() +
  annotate("text",  x = 0.1, y = max(df$UniqueBioEndoDNA, na.rm = T),  label = paste("Slope: ", round(slope_b, 3), "\nR-squared: ", round(rsquared_b, 2)), hjust = 0, vjust = 1,size = 6) +
  theme(legend.position = "none", text=element_text(size=12), axis.text.x = element_text(angle=90, hjust=1)) +
  labs(title="B.unique biological endogenous DNA %", x="screening shotgun",y="deep shotgun")
plotb

plotc <- ggplot(df, aes(x=SNPs_Covered_screen, y = SNPs_Covered, group=method)) +
  geom_point(size=4, colour="black", shape=21, aes(fill=MappableBioEndoDNA_SG)) +
  geom_smooth(method="lm", alpha=0.2, formula=y~x) +
  geom_abline(slope=1) +
  scale_fill_gradientn(colours = rainbow(3), trans="log10", n.breaks=12) +
  theme_light() +
  annotate("text",  x = min(df$SNPs_Covered_screen, na.rm = T), y = max(df$SNPs_Covered, na.rm = T), label = paste("Slope: ", round(slope_c, 2), "\nR-squared: ", round(rsquared_c, 2)), hjust = 0, vjust = 1,size = 6) +
  theme(legend.position = "none", text=element_text(size=12), axis.text.x = element_text(angle=90, hjust=1) ) +
  labs(title="C. 1240k SNPs covered", x="screening shotgun",y="deep shotgun")
plotc

plotd <- ggplot(df, aes(x=SNPsPer1MReadsToMapping_SG, y = SNPsPer1MReadsToMapping, group=method)) +
  geom_point(size=4, colour="black", shape=21, aes(fill=MappableBioEndoDNA_SG)) +
  geom_smooth(method="lm", alpha=0.2, formula=y~x) +
  geom_abline(slope=1) +
  scale_fill_gradientn(colours = rainbow(3), trans="log10", n.breaks=12, name = str_wrap("Mappable Biological Endogenous% (screening)", width = 11)) +
  theme_light() +
  annotate("text",  x = min(df$SNPsPer1MReadsToMapping_SG, na.rm = T), y = max(df$SNPsPer1MReadsToMapping, na.rm = T), label = paste("Slope: ", round(slope_d, 2), "\nR-squared: ", round(rsquared_d, 2)), hjust = 0, vjust = 1,size = 6) +
  theme(text=element_text(size=12), axis.text.x = element_text(angle=90, hjust=1), legend.position = "right", legend.key = element_rect(color = "grey", fill = "white"), legend.box.just = "right", legend.title = element_text(face = "bold", margin = margin(b = 5)) ) +
  guides(fill = guide_colorbar(barheight = 12)) +
  labs(title="D. 1240k SNPs covered per million reads into mapping", x="screening shotgun",y="deep shotgun", fill="")
plotd 

plota + plotb +plotc + plotd

ggsave("Twist_sg_SG.pdf", width = 18, height = 18, units="in", dpi=600)
