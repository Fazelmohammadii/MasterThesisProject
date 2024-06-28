########## Following the preprocessing script ##########

library(maftools)
library(ggridges)
#import mutation list from preprocess step as dataframe
allmutsnodup <- readRDS(allmutsnodup)
# this column is required for data input in Maftools package
allmutsnodup$Tumor_Seq_Allele2 <-""

##### SUB SET DNMT3A, TET2, ASXL1 Mutations ######
DNMT3A_AML <- allmutsnodup[allmutsnodup$Hugo_Symbol=="DNMT3A" & allmutsnodup$status %in% c("AML_only","shared"),]
DNMT3A_CH <- allmutsnodup[allmutsnodup$Hugo_Symbol=="DNMT3A" & allmutsnodup$status %in% c("CH_only","shared"),]
TET2_AML <- allmutsnodup[allmutsnodup$Hugo_Symbol=="TET2" & allmutsnodup$status %in% c("AML_only","shared"),]
TET2_CH <- allmutsnodup[allmutsnodup$Hugo_Symbol=="TET2" & allmutsnodup$status %in% c("CH_only","shared"),]
ASXL1_AML <- allmutsnodup[allmutsnodup$Hugo_Symbol=="ASXL1" & allmutsnodup$status %in% c("AML_only","shared"),]
ASXL1_CH <- allmutsnodup[allmutsnodup$Hugo_Symbol=="ASXL1" & allmutsnodup$status %in% c("CH_only","shared"),]

##### LILLIPOP PLOT #####
#read as MAF file
#DNMT3A
dn_aml <- read.maf(DNMT3A_AML)
dn_ch <-  read.maf(DNMT3A_CH)
#TET2
tt_aml <- read.maf(TET2_AML)
tt_ch <-  read.maf(TET2_CH)
#ASXL1
as_aml <- read.maf(ASXL1_AML)
as_ch <-  read.maf(ASXL1_CH)
#Pairwise Lollipop plot of each gene in CH and AML conditions
dnloll<- lollipopPlot2(m1 = dn_aml, m2 = dn_ch, gene = "DNMT3A", AACol1 = "Protein_Change", AACol2 = "Protein_Change", m1_name = "AML", m2_name = "CH",alpha = 0.4,showDomainLabel = F,domainBorderCol = NA, roundedRect = T,legendTxtSize = 1,pointSize = 2)
ttloll<- lollipopPlot2(m1 = tt_aml, m2 = tt_ch, gene = "TET2", AACol1 = "Protein_Change", AACol2 = "Protein_Change", m1_name = "AML", m2_name = "CH",alpha = 0.4,showDomainLabel = F,domainBorderCol = NA, roundedRect = T,legendTxtSize = 1,pointSize = 2)
asloll<- lollipopPlot2(m1 = as_aml, m2 = as_ch, gene = "ASXL1", AACol1 = "Protein_Change", AACol2 = "Protein_Change", m1_name = "AML", m2_name = "CH",alpha = 0.4,showDomainLabel = F,domainBorderCol = NA, roundedRect = T,legendTxtSize = 1,pointSize = 2)

##### MUT TYPE SPECTRUM #####
# Create the density ridge plot with percentage fill
#DNMT3A_AML
dnamlridge<- ggplot(DNMT3A_AML, aes(y = Variant_Classification, x = AA_position, fill = Variant_Classification)) +
  geom_density_ridges(scale = 3, rel_min_height = 0.01,  alpha = 0.7) +
  scale_fill_manual(values = c("blue","purple","green","red")) +
  theme_ridges() +
  guides(fill = guide_legend(reverse = TRUE)) + xlim(0,1100)+
  theme(axis.title.x=element_blank(),axis.title.y =element_blank())
#DNMT3A_CH
dnchridge<-ggplot(DNMT3A_CH, aes(y = Variant_Classification, x = AA_position, fill = Variant_Classification)) +
  geom_density_ridges(scale = 3, rel_min_height = 0.01,  alpha = 0.7) +
  scale_fill_manual(values = c("blue","purple","green","red")) +
  theme_ridges() +
  guides(fill = guide_legend(reverse = TRUE)) + xlim(0,1100)+
  theme(axis.title.x=element_blank(),axis.title.y =element_blank())
#TET2_AML
ttamlridge<- ggplot(TET2_AML, aes(y = Variant_Classification, x = AA_position, fill = Variant_Classification)) +
  geom_density_ridges(scale = 3, rel_min_height = 0.01,  alpha = 0.7) +
  scale_fill_manual(values = c("blue","purple","green","red")) +
  theme_ridges() +
  guides(fill = guide_legend(reverse = TRUE)) + xlim(0,2100)+
  theme(axis.title.x=element_blank(),axis.title.y =element_blank())
#TET2_CH
ttchridge<-ggplot(TET2_CH, aes(y = Variant_Classification, x = AA_position, fill = Variant_Classification)) +
  geom_density_ridges(scale = 3, rel_min_height = 0.01,  alpha = 0.7) +
  scale_fill_manual(values = c("blue","purple","green","red")) +
  theme_ridges() +
  guides(fill = guide_legend(reverse = TRUE)) + xlim(0,2100)+
  theme(axis.title.x=element_blank(),axis.title.y =element_blank())
#ASXL1_AML
asamlridge<- ggplot(ASXL1_AML, aes(y = Variant_Classification, x = AA_position, fill = Variant_Classification)) +
  geom_density_ridges(scale = 3, rel_min_height = 0.01,  alpha = 0.7) +
  scale_fill_manual(values = c("blue","purple","green","red")) +
  theme_ridges() +
  guides(fill = guide_legend(reverse = TRUE)) + xlim(0,1500)+
  theme(axis.title.x=element_blank(),axis.title.y =element_blank())
#ASXL1_CH
aschridge<-ggplot(ASXL1_CH, aes(y = Variant_Classification, x = AA_position, fill = Variant_Classification)) +
  geom_density_ridges(scale = 3, rel_min_height = 0.01,  alpha = 0.7) +
  scale_fill_manual(values = c("blue","purple","green","red")) +
  theme_ridges() +
  guides(fill = guide_legend(reverse = TRUE)) + xlim(0,1500)+
  theme(axis.title.x=element_blank(),axis.title.y =element_blank())

##### Venn diagram of CH_only, AML_only, and Shared #####
#subset each data frame by Protein Change feature as the most unique value of each variant
ch <- DNMT3A_CH$Protein_Change
aml <- DNMT3A_AML$Protein_Change

tch <- TET2_CH$Protein_Change
taml <- TET2_AML$Protein_Change

aaml<- ASXL1_AML$Protein_Change
ach<- ASXL1_CH$Protein_Change

l <- list(CH=ch,AML=aml)
t <- list(CH=tch,AML=taml)
a <- list(CH=ach,AML=aaml) 


#draw Venn Diagram
library(VennDiagram)
library(grid)

png(filename = "DN_venn.png", width = 3200, height = 3200, res = 300)

DNMT3A_Venn<- venn.diagram(x = l,filename = NULL,
                      output=T,
        category.names = c("CH (683)","AML (270)"),
        print.mode=c("raw","percent"),
       # Circles
        lwd = 1,
        lty = 'solid',
        cex=3,
        fill = c(alpha(ch_color,0.4),alpha(aml_color,0.4)),
        col='black',
        cat.cex=4, cat.fontface = "bold",cat.col= c(ch_color,aml_color),
       cat.pos= c(0.25,-0.15),
       cat.dist= c(0.03),
       cat.default.pos = "outer")
grid.draw(DNMT3A_Venn)

dev.off()


png(filename = "TT_venn.png", width = 3200, height = 3200, res = 300)

TET2_Venn<- venn.diagram(x = t,filename = NULL,
                      output=T,
        category.names = c("CH (451)","AML (263)"),
        print.mode=c("raw","percent"),
       # Circles
        lwd = 1,
        lty = 'solid',
        cex=3,
        fill = c(alpha(ch_color,0.4),alpha(aml_color,0.4)),
        col='black',
        cat.cex=4, cat.fontface = "bold",cat.col= c(ch_color,aml_color),
       cat.pos= c(0.25,-0.15),
       cat.dist= c(0.04),
       cat.default.pos = "outer")

grid.draw(TET2_Venn)

dev.off()


png(filename = "AS_venn.png", width = 3200, height = 3200, res = 300)

ASXL1_Venn<- venn.diagram(x = a,filename = NULL,
                      output=T,
        category.names = c("CH (451)","AML (263)"),
        print.mode=c("raw","percent"),
       # Circles
        lwd = 1,
        lty = 'solid',
        cex=3,
        fill = c(alpha(ch_color,0.4),alpha(aml_color,0.4)),
        col='black',
        cat.cex=4, cat.fontface = "bold",cat.col= c(ch_color,aml_color),
       cat.pos= c(0.25,-0.15),
       cat.dist= c(0.03),
       cat.default.pos = "outer")

grid.draw(ASXL1_Venn)

dev.off()

##### Output the mutation list of each gene #####
library(dplyr)
status_output <- function(named_list, filename) {

  CH <- named_list$CH
  AML <- named_list$AML
  
  
  unique_CH <- setdiff(CH, AML)
  unique_AML <- setdiff(AML, CH)
  shared <- intersect(CH, AML)
  
  df_unique_CH <- data.frame(Element = unique_CH, Category = "ONLY_CH")
  df_unique_AML <- data.frame(Element = unique_AML, Category = "ONLY_AML")
  df_shared <- data.frame(Element = shared, Category = "SHARED")
  
  df <- bind_rows(df_unique_CH, df_unique_AML, df_shared)
  

  write.csv(df, file = filename, row.names = FALSE)
}

status_output(l, "DN_mutlist.csv")
status_output(t, "TT_mutlist.csv")
status_output(a, "AS_mutlist.csv")
