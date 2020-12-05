

sd = function(input){
  control = vector()
  TMEM = vector()
  for(i in 1:nrow(input)){
    for(j in 6:7){
      for(k in 8:9){
        control[i] = sqrt(sum(as.numeric(input[i,j])-as.numeric(input[i,12])**2)/2)
        TMEM[i] = sqrt(sum(as.numeric(input[i,k])-as.numeric(input[i,13])**2)/2)
        #input = add_column(input, GFPsd = control[i], .before = 13)
        #input = add_column(input, TMEMsd = TMEM[i], .after = 14)
        control <<- control
        TMEM <<- TMEM
        return(c(control[1:5],TMEM[1:5])) 
      }
    }
  }
}
sd(input = MS_DF_filter6
)

sqrt(sum(abs(as.numeric(MS_DF_filter6[1,c(8:9)])-as.numeric(MS_DF_filter6[1,13]))**2)/2)

as.numeric(MS_DF_filter6[1,13])/sqrt(sum(abs(as.numeric(MS_DF_filter6[1,c(8:9)])-as.numeric(MS_DF_filter6[1,13]))**2)/2)


MS_DF = read.csv("M:/Erik/Data/Omics/MS/20201012_Samples_View_Report_4441.csv")
## Re-name the column names by the first row
colnames(MS_DF) = MS_DF[1,]
## Delete that row
MS_DF = MS_DF[-1,]
## Re-order the row index
rownames(MS_DF) = NULL

## Replace missing integers with 0s, and question marks with blanks
for (i in 4:8){
  for (j in 1:nrow(MS_DF)){
    if (MS_DF[j,i] == ""){
      MS_DF[j,i] = 0
    } else if (MS_DF[j,i] == "?"){
      MS_DF[j,i] = ""
    }
  }
}

## Co-erce the data to be integers, not characters
for (i in 5:8){
  for (j in 1:nrow(MS_DF)){
    MS_DF[j,i] = as.numeric(MS_DF[j,i])
  }
}

## Remove the last row
MS_DF = MS_DF[-(nrow(MS_DF)),]
## Remove the '#' column
MS_DF = MS_DF[ , -1]
## Sort the data
MS_DF = MS_DF[order(MS_DF$myc_GFP_1, decreasing = FALSE), ]


## Add a column that will be comprised just of human gene names
MS_DF = add_column(MS_DF, HumanGeneID = "", .before = 3)
## Copy the strings
MS_DF$HumanGeneID = MS_DF$HumanGeneID %>% gsub(x = MS_DF$`Identified Proteins (6608)`, replacement = "")


## Add a column that will be comprised just of mouse gene names
MS_DF = add_column(MS_DF, MouseGeneID = "", .before = 4)
## Copy the strings
MS_DF$MouseGeneID = MS_DF$MouseGeneID %>% gsub(x = MS_DF$`Identified Proteins (6608)`, replacement = "")


## Extract the human protein's human gene name
MS_DF$HumanGeneID = str_extract_all(MS_DF$HumanGeneID, "GN=[:alnum:]{1,}")
## Remove the identifier
MS_DF$HumanGeneID = MS_DF$HumanGeneID %>% gsub(x = MS_DF$HumanGeneID, pattern = "GN=", replacement = "")

## Extract the protein's gene name
MS_DF$MouseGeneID = str_extract_all(MS_DF$MouseGeneID, "GN=[:alnum:]{1,}")
## Remove the identifier
MS_DF$MouseGeneID = MS_DF$MouseGeneID %>% gsub(x = MS_DF$MouseGeneID, pattern = "GN=", replacement = "")
## Put the gene name in MM9 EntrezID form (note, however, these were human proteins; consider leaving)
MS_DF$MouseGeneID = paste(substr(MS_DF$MouseGeneID, 1, 1),
                          tolower(substr(MS_DF$MouseGeneID, 2, 10)), sep = "")

## Remove extra strings from the "Identified Proteins" column
MS_DF$`Identified Proteins (6608)` = MS_DF$`Identified Proteins (6608)` %>% gsub(x = MS_DF$`Identified Proteins (6608)`, pattern = " OS.+.?$", replacement = "")

## Fix p53 in mouse
MS_DF$MouseGeneID[which(MS_DF$`Accession Number` == "P53_HUMAN")] = "Trp53"


## Remove the decoy proteins
MS_DF_filter = MS_DF %>% filter(!grepl(MS_DF$`Accession Number`, pattern = "DECOY"))

## Remove proteins under 10 counts in TMEM sample
MS_DF_filter1 = MS_DF_filter[ -c(which(as.numeric(MS_DF_filter$TMEM184B_myc_1) < 5)) , ]
MS_DF_filter2 = MS_DF_filter1[ -c(which(as.numeric(MS_DF_filter1$TMEM184B_myc_2) < 5)) , ]
MS_DF_filter3 = MS_DF_filter2[ -c(which(as.numeric(MS_DF_filter2$TMEM184B_myc_1) < 10)) , ]
#MS_DF_filter4 = MS_DF_filter3[ -c(which(as.numeric(MS_DF_filter3$TMEM184B_myc_2) < 10)) , ]

## Remove Keratin and Tubulin, common contaminants or false positives
MS_DF_filter4 = MS_DF_filter3 %>% filter(!grepl(MS_DF_filter3$`Identified Proteins (6608)` , pattern = "Keratin.+?$|Tubulin.+?$"))

## Remove the "kDa" from the MW column
MS_DF_filter4$`Molecular Weight` = as.numeric(gsub(" kDa","", MS_DF_filter4$`Molecular Weight`))
## Scale all counts by MW
MS_DF_filter4$myc_GFP_1 = as.numeric(MS_DF_filter4$myc_GFP_1)/as.numeric(MS_DF_filter4$`Molecular Weight`)
MS_DF_filter4$myc_GFP_2 = as.numeric(MS_DF_filter4$myc_GFP_2)/as.numeric(MS_DF_filter4$`Molecular Weight`)
MS_DF_filter4$TMEM184B_myc_1 = as.numeric(MS_DF_filter4$TMEM184B_myc_1)/as.numeric(MS_DF_filter4$`Molecular Weight`)
MS_DF_filter4$TMEM184B_myc_2 = as.numeric(MS_DF_filter4$TMEM184B_myc_2)/as.numeric(MS_DF_filter4$`Molecular Weight`)

## Average the counts from both conditions
MS_DF_filter4 = add_column(MS_DF_filter4, GFPavg = (as.numeric(MS_DF_filter4$myc_GFP_1) + as.numeric(MS_DF_filter4$myc_GFP_2))/2)
MS_DF_filter4 = add_column(MS_DF_filter4, TMEMavg = as.numeric(MS_DF_filter4$TMEM184B_myc_1) + as.numeric(MS_DF_filter4$TMEM184B_myc_2)/2)

## Find the Ratio of TMEM1:TMEM2
TMEM_RATIO = as.numeric(MS_DF_filter4[ which(MS_DF_filter4$MouseGeneID == "Tmem184b"), 8 ]) / as.numeric(MS_DF_filter4[ which(MS_DF_filter4$MouseGeneID == "Tmem184b"), 9 ])
## Create a column that stores the given protein's ratio of spectral counts in TMEM samples (scale the counts from TMEM 1: TMEM 2)
MS_DF_filter5 = add_column(MS_DF_filter4, TMEM_RATIO = as.numeric(MS_DF_filter4$TMEM184B_myc_1) / as.numeric(MS_DF_filter4$TMEM184B_myc_2))
## Create a column to include whether the protein's gene was DE in aDRG TMEM mutants
MS_DF_filter5$DEG_in_aDRG = ""


## Remove proteins that have more counts in the 2nd TMEM sample
MS_DF_filter5 = subset(MS_DF_filter5, MS_DF_filter5$TMEM_RAT > 1.0)

## Add the scale that regularizes prey proteins that are farther away from TMEM's ratio
MS_DF_filter6 = add_column(MS_DF_filter5, scale = as.numeric((MS_DF_filter5$TMEMavg*MS_DF_filter5$TMEM_RATIO)-MS_DF_filter5$TMEMavg*sqrt((MS_DF_filter5$TMEM_RATIO-TMEM_RATIO)**2)))






MS_DF_filter6$TMEMavg/(MS_DF_filter6$GFPavg+.1)