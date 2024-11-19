# Import data using the data hub
# Data Hub List
# All datasets are available at https://xenabrowser.net/datapages/.

install.packages("UCSCXenaTools")
library(UCSCXenaTools)

data(XenaData)

head(XenaData)

#-- Get clinical data from 
#----------------------------------------------------------------
#Vie
# The options in XenaFilter function support Regular Expression
XenaGenerate(subset = XenaHostNames=="tcgaHub") %>% 
        XenaFilter(filterDatasets = "HNSC") -> df_todo

%>%
        XenaFilter(filterDatasets = "phenotype")  -> df_todo

XenaGenerate(subset = XenaHostNames=="gdcHub") %>% 
        XenaFilter(filterDatasets = "HNSC") %>%
        XenaFilter(filterDatasets = "phenotype")  -> df_todo

df_todo

XenaQuery(df_todo) %>%
        XenaDownload() -> xe_download
cli = XenaPrepare(xe_download)
class(cli)
dim(cli)
names(cli)[grepl("race",names(cli))]

write.csv(cli,"/")

#-----------------------------------------------------------------
#Vie
# The options in XenaFilter function support Regular Expression
XenaGenerate(subset = XenaHostNames=="gdcHub") %>% 
 XenaFilter(filterDatasets = "HNSC") %>%
 XenaFilter(filterDatasets = "phenotype")  -> df_todo

XenaGenerate(subset = XenaHostNames=="tcgaHub") %>% 
        XenaFilter(filterDatasets = "HNSC") %>%
        XenaFilter(filterDatasets = "clinicalMatrix")  -> df_todo

df_todo

XenaQuery(df_todo) %>%
 XenaDownload() -> xe_download

#Prepare data into R for analysis.
# download 08/20/19
cli = XenaPrepare(xe_download)
class(cli)

 write.table(cli,"Clinical_data_11_25_19.txt")
 write.csv(cli,"Clinical_data_11_25_19.csv")
# CliTp <- read.csv("Clinical_data_08_20_19.csv")
#Browse data
XenaBrowse(df_todo)

XenaBrowse(df_todo, type = "cohort")

#*****************************************************
# Download survival data
# The options in XenaFilter function support Regular Expression
XenaGenerate(subset = XenaHostNames=="gdcHub") %>% 
        XenaFilter(filterDatasets = "HNSC") %>%
        XenaFilter(filterDatasets = "survival")  -> df_todoSurv

df_todoSurv

XenaQuery(df_todoSurv) %>%
        XenaDownload() -> xe_downloadSurv

#Prepare data into R for analysis.
# download 08/20/19
cliSurv = XenaPrepare(xe_downloadSurv)
class(cliSurv)
write.table(cliSurv,"Survival_data_08_20_19.txt")
write.csv(cliSurv,"Survival_data_08_20_19.csv")
#Browse data
XenaBrowse(df_todoSurv)

XenaBrowse(df_todoSurv, type = "cohort")

#*****************************************************
#Download the miRNA data 
#*****************************************************
XenaGenerate(subset = XenaHostNames=="gdcHub") %>% 
 XenaFilter(filterDatasets = "HNSC") %>%
 XenaFilter(filterDatasets = "miRNA")  -> df_miRNA

XenaGenerate(subset = XenaHostNames=="tcgaHub") %>% 
        XenaFilter(filterDatasets = "HNSC") %>%
        XenaFilter(filterDatasets = "miRNA_HiSeq_gene")  -> df_miRNA

df_miRNA

XenaQuery(df_miRNA) %>%
 XenaDownload() -> xemiRNA_download

#---  


#Prepare data into R for analysis.
# download 08/20/19
gdcHub_miRNA = XenaPrepare(xemiRNA_download)
class(gdcHub_miRNA)
write.table(gdcHub_miRNA,"miRNA_data_11_25_19.txt")
write.csv(gdcHub_miRNA,"miRNA_data_11_25_19.csv")

#Browse data
XenaBrowse(df_miRNA)

XenaBrowse(df_miRNA, type = "cohort")


#*****************************************************
#Download the mRNA data 
#*****************************************************
XenaGenerate(subset = XenaHostNames=="gdcHub") %>% 
        XenaFilter(filterDatasets = "HNSC") %>%
        XenaFilter(filterDatasets = "htseq_counts")  -> df_mRNACounts

XenaGenerate(subset = XenaHostNames=="tcgaHub") %>% 
        XenaFilter(filterDatasets = "HNSC") %>%
        XenaFilter(filterDatasets = "miRNA_HiSeq_gene")  -> df_mRNACounts


df_mRNACounts

XenaQuery(df_mRNACounts) %>%
        XenaDownload() -> df_mRNACounts_download

#Prepare data into R for analysis.
# download 08/20/19
gdcHub_mRNACounts = XenaPrepare(df_mRNACounts_download)
class(gdcHub_mRNACounts)
write.table(gdcHub_mRNACounts,"mRNA_datanormalized_log2RSEM_11_25_19.txt")
write.csv(gdcHub_mRNACounts,"mRNA_datanormalized_log2RSEM_11_25_19.csv")


#*****************************************************
XenaGenerate(subset = XenaHostNames=="gdcHub") %>% 
        XenaFilter(filterDatasets = "HNSC") %>%
        XenaFilter(filterDatasets = "htseq_fpkm.")  -> df_mRNAfpkm

df_mRNAfpkm
df_mRNAfpkm@datasets<-df_mRNAfpkm@datasets[1]
df_mRNAfpkm@datasets

XenaQuery(df_mRNAfpkm) %>%
        XenaDownload() -> df_mRNAfpkm_download

#Prepare data into R for analysis.
# download 08/20/19
gdcHub_mRNAfpkm = XenaPrepare(df_mRNAfpkm_download)
class(gdcHub_mRNAfpkm)
write.table(gdcHub_mRNAfpkm,"mRNA_datafpkm_08_20_19.txt")
write.csv(gdcHub_mRNAfpkm,"miRNA_datafpkm_08_20_19.csv")






#Browse data
XenaBrowse(df_miRNA)

XenaBrowse(df_miRNA, type = "cohort")




# The options in XenaFilter function support Regular Expression
XenaGenerate(subset = XenaHostNames=="tcgaHub") %>% 
 XenaFilter(filterDatasets = "clinical") %>% 
 XenaFilter(filterDatasets = "LUAD|LUSC|LUNG") -> df_todo

df_todo
