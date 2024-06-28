#The mutation lists from differnt studies of CH are merged and harmonized using the “Lift Genome Annotation” tool of the UCSC Genome Browser and Variant Effect Predictor (VEP). Single VCF file of merged mutation list is created for CH and AML conditions.
#For checking spelling errors and unmatched annotations VCFs are converted to CSV files. 

library(ggplot2)

##### IMPORT DATASETS ######
CHnoUKB<- read.csv("CHold.csv", encoding = "UTF-8")
AMLnoLUE <- read.csv("AMLold.csv", encoding = "UTF-8")

##### FILTERATIONS ######
#Removing mutations with very low VAF
CHnoUKB <- CHnoUKB[CHnoUKB$VAF>=0.05,]
AMLnoLUE <- AMLnoLUE[AMLnoLUE$VAF>=0.05,]
#remove redundency of mutations whitin the 
CHnodup<- unique(CHnoUKB[,c(2:9,11)],fromLast = F)
AMLnodup <- unique(AMLnoLUE[,c(2:9,11)],fromLast = F)

##### VAF FREQUENCY CHECK ######
#select three genes with their VAF in CH population
vafch <- CHnoUKB[CHnoUKB$Hugo_Symbol %in% c("DNMT3A", "TET2", "ASXL1"),c("Hugo_Symbol","VAF")]
vafch$state <- "CH"
#harmonizing VAF values
for (i in 1:nrow(vafch)){
  if(vafch$VAF[i] > 1){
    vafch[[VAF]][i] <- vafch[[VAF]][i] /100 
  }
}
#select three genes with their VAF in AML population
vafaml<- AMLnoLUE[AMLnoLUE$Hugo_Symbol %in% c("DNMT3A", "TET2", "ASXL1"),c("Hugo_Symbol","VAF")]
vafaml$state <- "AML"
vafaml$VAF <- as.numeric(vafaml$VAF)
#harmonizing VAF values
for (i in 1:nrow(vafaml)){
  if(vafaml$VAF[i] > 1){
    vafaml[["VAF"]][i] <- vafaml[["VAF"]][i] / 100 
  }
}
#bind both dataframes
vafbx <- rbind(vafch,vafaml)
vafbx <- vafbx[vafbx$VAF >= 0.05,]

## VAF plot ###
# Assigning colors
ch_color <- "#66C2A5"  # Light blue or teal color for CH
ch_color_high <- "#84ceb7"
ch_color_low <- "#eff8f6"
aml_color <- "#B2182B"  # Dark red or maroon color for AML
aml_color_high <- "#c14655"
aml_color_low <-"#f7e7e9"
# Plot
VAFBXP<-ggplot(vafbx,aes(x=Hugo_Symbol , y=VAF , fill=state))+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~ Hugo_Symbol ,nrow=1, scale="free_x")+
  scale_fill_manual(values = c(aml_color,ch_color)) +
  ylim(0,.7)+
  theme_gray()+
  theme(axis.text.x=element_blank(),axis.ticks.y=element_blank(),axis.title.x=element_blank())
  #    axis.text.y=element_blank(),
  #    axis.ticks.y=element_blank())

##### MUTATION COUNTS SUBSETS #####
allmuts <- merge(CHnodup,AMLnodup, all= TRUE)
allmutsnodup <- unique(allmuts, FromLast=FALSE)
#labeling mutations as their conditions:
allmutsnodup$LABEL <- ifelse(allmutsnodup$status=="CH_only","CH","AML")
#sub setting AML_only, CH_only, and Shared mutations
freqch<- data.frame(table(allmutsnodup[allmutsnodup$status=="CH_only","Hugo_Symbol"]))
freqch$state <- "CH_only"
freqsh<- data.frame(table(allmutsnodup[allmutsnodup$status=="shared","Hugo_Symbol"]))
freqsh$state <- "shared"
freqaml<- data.frame(table(allmutsnodup[allmutsnodup$status=="AML_only","Hugo_Symbol"]))
freqaml$state <- "AML_only"
mutsfreq <- rbind(freqch,freqaml,freqsh)
#count mutations per genes in subsets
genebars <- mutsfreq %>% mutate(Freq = ifelse(state%in%c("AML_only","shared"), Freq,-1*Freq))
names(genebars)[1]<- "genes"

#break to show counts of different subset at once
y_breaks <- pretty(genebars$Freq)
bar_label <- genebars
bar_label$Freq <- abs(bar_label$Freq)
#plot
genefreq<- ggplot(genebars,
       aes(x = reorder(genes, round(abs(Freq)/sum(abs(Freq))*100)),
           y = Freq,
           fill = state)) +
  geom_bar(stat = 'identity', position = 'stack') +
  coord_flip() +
  scale_y_continuous(breaks = y_breaks,
                     labels = abs(y_breaks)) +
  geom_text(aes(label = paste0(abs(Freq))),
            position = position_stack(vjust = 0.4), size = 4) +
  ylab("frequency") +
  xlab("genes") +
  scale_fill_manual(values = c(aml_color_high,ch_color,aml_color_low))+
  theme_minimal()


#####  #####
