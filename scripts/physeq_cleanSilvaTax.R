library(glue)
library(stringi)
library(varhandle)

physeq_cleanSilvaTax=function(physeq, append_prefix=TRUE) {

tax.tbl=as.data.frame(as(tax_table(physeq),"matrix"))
tax.tbl$ASV=rownames(tax.tbl)

#concatenate Genus-Species assignments
for (i in seq_along(rownames(tax.tbl))) {

  if (is.na(tax.tbl$Species[i])) {
  tax.tbl$Genussp[i]=paste0(tax.tbl$Genus[i])
  } else {
    tax.tbl$Genussp[i]=paste0(tax.tbl$Genus[i]," ",tax.tbl$Species[i])
    #remove duplicate genus names
    tax.tbl$Genussp[i]=paste0(glue_collapse(intersect(unlist(str_split(tax.tbl$Genussp[i]," ")),unlist(str_split(tax.tbl$Genus[i]," "))),sep=" ")," ",glue_collapse(setdiff(unlist(str_split(tax.tbl$Species[i]," ")),unlist(str_split(tax.tbl$Genus[i]," "))),sep=" "))

  }

}

#for (i in 1:ncol(tax.tbl)) {
    
#    tax.tbl[,i]=unfactor(as.data.frame(tax.tbl[,i]))
    
#}

tax.tbl$Genussp[which(is.na(tax.tbl$Genussp))]="NA"
tax.tbl$Genussp[which(tax.tbl$Genussp=="NA")]="NAN"

#fill Genus-level NAs with next higher tax (Family)
tax.tbl$Family=unfactor(tax.tbl$Family)

tax.tbl$Genussp[which(tax.tbl$Genussp=="NAN")]=tax.tbl$Family[which(tax.tbl$Genussp=="NAN")]
#tax.tbl$Genussp[tax.tbl$Genussp=="NA"]=tax.tbl$Family[tax.tbl$Genussp=="NA"]	

#fill Genus-level ambigs with next higher tax (Family)
tax.tbl$Genussp[which(tax.tbl$Genussp=="uncultured")]=tax.tbl$Family[which(tax.tbl$Genussp=="uncultured")]
tax.tbl$Genussp[which(tax.tbl$Genussp=="uncultured bacterium")]=tax.tbl$Family[which(tax.tbl$Genussp=="uncultured bacterium")]
tax.tbl$Genussp[which(tax.tbl$Genussp=="uncultured organism")]=tax.tbl$Family[which(tax.tbl$Genussp=="uncultured organism")]
tax.tbl$Genussp[which(tax.tbl$Genussp=="uncultured human gut metagenome")]=tax.tbl$Family[which(tax.tbl$Genussp=="uncultured human gut metagenome")]
tax.tbl$Genussp[which(tax.tbl$Genussp=="uncultured gut metagenome")]=tax.tbl$Family[which(tax.tbl$Genussp=="uncultured gut metagenome")]

tax.tbl$Genussp[which(tax.tbl$Genussp=="gut metagenome")]=tax.tbl$Family[which(tax.tbl$Genussp=="gut metagenome")]

tax.tbl$Genussp[which(tax.tbl$Genussp=="unidentified")]=tax.tbl$Family[which(tax.tbl$Genussp=="unidentified")]

#fill Family-level NAs with next higher tax (Order)
tax.tbl$Genussp[which(is.na(tax.tbl$Genussp))]="NA"
tax.tbl$Genussp[which(tax.tbl$Genussp=="NA")]="NAN"

tax.tbl$Genussp[which(tax.tbl$Genussp=="NAN")]=tax.tbl$Order[which(tax.tbl$Genussp=="NAN")]

#fill Family-level ambigs with next higher tax (Order)
tax.tbl$Genussp[which(tax.tbl$Genussp=="uncultured")]=tax.tbl$Order[which(tax.tbl$Genussp=="uncultured")]
tax.tbl$Genussp[which(tax.tbl$Genussp=="uncultured bacterium")]=tax.tbl$Order[which(tax.tbl$Genussp=="uncultured bacterium")]
tax.tbl$Genussp[which(tax.tbl$Genussp=="uncultured organism")]=tax.tbl$Order[which(tax.tbl$Genussp=="uncultured organism")]
tax.tbl$Genussp[which(tax.tbl$Genussp=="uncultured human gut metagenome")]=tax.tbl$Order[which(tax.tbl$Genussp=="uncultured human gut metagenome")]
tax.tbl$Genussp[which(tax.tbl$Genussp=="uncultured gut metagenome")]=tax.tbl$Order[which(tax.tbl$Genussp=="uncultured gut metagenome")]

tax.tbl$Genussp[which(tax.tbl$Genussp=="gut metagenome")]=tax.tbl$Order[which(tax.tbl$Genussp=="gut metagenome")]

tax.tbl$Genussp[which(tax.tbl$Genussp=="unidentified")]=tax.tbl$Order[which(tax.tbl$Genussp=="unidentified")]

#fill out all succeeding NAs with next higher tax (Class upwards)
tax.tbl$Genussp[which(is.na(tax.tbl$Genussp))]="NA"
tax.tbl$Genussp[which(tax.tbl$Genussp=="NA")]="NAN"

tax.tbl$Genussp[which(tax.tbl$Genussp=="NAN")]=tax.tbl$Class[which(tax.tbl$Genussp=="NAN")]

tax.tbl$Genussp[which(is.na(tax.tbl$Genussp))]="NA"
tax.tbl$Genussp[which(tax.tbl$Genussp=="NA")]="NAN"

tax.tbl$Genussp[which(tax.tbl$Genussp=="NAN")]=tax.tbl$Phylum[which(tax.tbl$Genussp=="NAN")]

tax.tbl$Genussp[which(is.na(tax.tbl$Genussp))]="NA"
tax.tbl$Genussp[which(tax.tbl$Genussp=="NA")]="NAN"

tax.tbl$Genussp[which(tax.tbl$Genussp=="NAN")]=tax.tbl$Kingdom[which(tax.tbl$Genussp=="NAN")]

#clean names
tax.tbl$Genussp=gsub("\\[","",tax.tbl$Genussp)
tax.tbl$Genussp=gsub("\\]","",tax.tbl$Genussp)
tax.tbl$Genussp=gsub(" bacterium","",tax.tbl$Genussp)# 
tax.tbl$Genussp=gsub(" organism","",tax.tbl$Genussp)
tax.tbl$Genussp=gsub(" group","",tax.tbl$Genussp)

#fix long entries (taxon-specific)
tax.tbl$Genussp=gsub("Bacteroides massiliensis B84634 = Timone 84634 DSM 17679 JCM 13223","Bacteroides massiliensis B84634",tax.tbl$Genussp)
tax.tbl$Genussp=gsub("Eubacterium coprostanoligenes human gut metagenome","Eubacterium coprostanoligenes",tax.tbl$Genussp)
tax.tbl$Genussp=gsub("Actinomyces Chlamydia trachomatis","Chlamydia trachomatis",tax.tbl$Genussp)
tax.tbl$Genussp=gsub("Cutibacterium Propionibacterium sp.","Cutibacterium sp.",tax.tbl$Genussp)

if (append_prefix==TRUE) {
#tax.tbl$Genussp=make.unique(tax.tbl$Genussp)
#tax.tbl$Genussp=if_else(duplicated(tax.tbl$Genussp),paste0(tax.tbl$Genussp,".ASV",stri_sub(tax.tbl$ASV,-4,-1)),tax.tbl$Genussp) 
tax.tbl$Genussp=paste0(tax.tbl$Genussp,".ASV",stri_sub(tax.tbl$ASV,-4,-1))

}

else if (append_prefix==FALSE) {
  
}

#re-order tables
tax.tbl$Lowest_taxon=tax.tbl$Genussp
tax.tbl=tax.tbl[,c("Kingdom","Phylum","Class","Order","Family","Genus","Species","ASV","Lowest_taxon")]

}
