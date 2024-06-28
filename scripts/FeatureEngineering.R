#####Following Preprocess Step #####
##### ADD DOMAIN FEATURE #####
#DNMT3A
dout<- allmutsnodup[allmutsnodup$Hugo_Symbol=="DNMT3A",]
dim(dout)
domain <- c()
for(i in dout$AA_position){
  if(i>= 628 & i<=912){
    domain <- c(domain,"MTase")
  }else{
    if(i>=476 & i<=614){
      domain <- c(domain,"ADD")
    }else{
      if(i>=278 & i<=427){
        domain<- c(domain,"PWP")
      }else{
        domain <- c(domain,"out_domain")
      }
    }
  }
}
dout$domain <- domain
#TET2
tout <- allmutsnodup[allmutsnodup$Hugo_Symbol=="TET2",]
domainTT <- c()
for(i in tout$AA_position){
  if(i>= 1047 & i<=1238){
    domainTT <- c(domainTT,"CRD")
  }else{
    if(i>=1238 & i<=1383){
      domainTT <- c(domainTT,"DSBH1")
    }else{
      if(i>=1766 & i<=1856){
        domainTT<- c(domainTT,"DSBH2")
      }else{
        domainTT<- c(domainTT,"out_domain")
      }
    }
  }
}
tout$domain <- domainTT

##### V TYPE FEATURE #####
#to make the Variant_type feature more specific here we change the feature categories into Indels with different steps of frame shift
modify_Vtype <- function(data) {
  for (i in 1:nrow(data)) {
    
      if (abs(data$End_Position[i] - data$Start_Position[i]) == 1 & grepl("Frame", data$Variant_Classification[i], fixed = TRUE)) {
        data$Variant_Type[i] <- paste0(data$Variant_Type[i], "_1bp")
      } else if (abs(data$End_Position[i] - data$Start_Position[i]) == 0 & grepl("Frame", data$Variant_Classification[i], fixed = TRUE)){
       data$Variant_Type[i] <- paste0(data$Variant_Type[i], "_1bp")
      } else if (abs(data$End_Position[i] - data$Start_Position[i]) == 2 & grepl("Frame", data$Variant_Classification[i], fixed = TRUE)) {
        data$Variant_Type[i] <- paste0(data$Variant_Type[i], "_2bp")
      } else if (abs(data$End_Position[i] - data$Start_Position[i]) >= 3 & grepl("Frame", data$Variant_Classification[i], fixed = TRUE)) {
        data$Variant_Type[i] <- paste0(data$Variant_Type[i], "_mbp")
      }
  }
  return(data)
}

dout <- modify_Vtype(dout)
tout <- modify_Vtype_2(tout)

##### SIGNATURE BASED FEATURES #####
#adding neighboring sequence of SNVs 
contexFunc<- function(data){

  library(BSgenome)
  library(rtracklayer)
  library(BSgenome.Hsapiens.UCSC.hg38)

  contextref <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38,
                                                as.character(gsub(" ", "", paste("chr", data$Chromosome))),
                                                as.integer(data$Start_Position)-2, as.integer(data$Start_Position)+2))

  contextalt <- contextref
  substr(contextalt,3,3) <- data$Alternative_Allele

  data <- cbind(data,contextref,contextalt)
  data$Start_Position <- as.integer(data$Start_Position)
  contextl1 <- as.integer(data$Start_Position)-2
  contextl2 <- as.integer(data$Start_Position)-1
  contextl3 <- as.integer(data$Start_Position)+1
  contextl4 <- as.integer(data$Start_Position)+2

  data <- cbind(data,contextl1,contextl2,contextl3,contextl4)

  sbn<- substr(data$contextref,2,2)
  sbp<- substr(data$contextref,4,4)

  data$sbn <- sbn
  data$sbp<- sbp

  return(data)

}
#adding substitution function
subFunc <- function(data){
  data$subs <- gsub(" ","", paste(data$Reference_Allele,">",data$Alternative_Allele))

  data$subs <- gsub("A>C","T>C",data$subs)
  data$subs <- gsub("A>G","T>G",data$subs)
  data$subs <- gsub("A>T","T>A",data$subs)

  data$subs <- gsub("G>C","C>G",data$subs)
  data$subs <- gsub("G>T","C>T",data$subs)
  data$subs <- gsub("G>A","C>A",data$subs)

  return(data)
}

dout <- subFunc(dout)
dout <- contexFunc(dout)
tout <- subFunc(tout)
tout <- contexFunc(tout)

##### TRIMER FEAFUTRES #####
#add trimer feature including neighboring neucleotide and the mutated one in SNVs
process_trimers <- function(df) {
  trimers <- gsub(" ", "", paste(df$sbn, substr(df$subs, 1, 1), df$sbp))
  df$trimers <- trimers
  df$trimers[grep("n", df$trimers)] <- "no_triplet"
  
  return(df)
}

dout <- process_trimers(dout)
tout <- process_trimers(tout)


##### AA FEATURES #####
process_protein_change <- function(df) {
  AAref <- substr(df$Protein_Change, 1, 1)
  AAalt <- gsub(".*\\d(.).*", "\\1", df$Protein_Change)
  
  df <- cbind(df, AAref, AAalt)
  
  AAchange <- paste(df$AAref, df$AAalt, sep = ">")
  df$AAchange <- AAchange
  
  return(df)
}

dout <- process_protein_change(dout)
tout <- process_protein_change(tout)


##### AA PHYSICOCHEMICAL FEATURES #####
# dictrinary for defining each of 22 amino acid physicochemical features imported from references
AAdict <- read.csv("C:\\Users\\fazel\\...\\resources\\AA Dict.csv")
names(AAdict)[1] <- "AA"
AAdict$pI <- as.numeric(AAdict$pI)
AAdict$Hydropobicity <- as.numeric(AAdict$Hydropobicity)
AAdict$Molecular.Weight <- as.numeric(AAdict$Molecular.Weight)

#DNMT3A
#reference AA
dout <- left_join(dout,AAdict, by=c("AAref"= "AA"))
names(dout)[c(11,12,13,14,15,16)] <- c("polarity_ref","charge_ref","sidechain_ref","pI_ref","hydrophobicity_ref","MW_ref")
#Alternative AA
dout <- left_join(dout,AAdict, by=c("AAalt"= "AA"))
names(dout)[c(17,18,19,20,21,22)] <- c("polarity_alt","charge_alt","sidechain_alt","pI_alt","hydrophobicity_alt","MW_alt")

#TET2
#reference AA
tout <- left_join(tout,AAdict, by=c("AAref"= "AA"))
names(tout)[c(11,12,13,14,15,16)] <- c("polarity_ref","charge_ref","sidechain_ref","pI_ref","hydrophobicity_ref","MW_ref")
#Alternative AA
tout <- left_join(tout,AAdict, by=c("AAalt"= "AA"))
names(tout)[c(17,18,19,20,21,22)] <- c("polarity_alt","charge_alt","sidechain_alt","pI_alt","hydrophobicity_alt","MW_alt")


#####ONE HOT ENCODING FUNCTION#####
encFunc <- function(data){
  require(caret)
  df <- data[,-ncol(data)] #split labels
  dummy <- dummyVars(" ~ .", data=df)
  final_df <- data.frame(predict(dummy, newdata=df)) #one hot encoding
  final_df <- cbind(final_df,data[ncol(data)]) #add labels again :)
  return(final_df)
}
