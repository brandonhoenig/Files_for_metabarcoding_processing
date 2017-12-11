#filtering reference datafile
#referencefile.txt is the outputfile from downloading the TSV data of your reference sequences from Barcode of Life Database#

#install necessary packages
install.packages("data.table")
install.packages("stringr")
install.packages("seqinr")

#load necessary libraries
library(data.table)
library(stringr)
library(seqinr)

#reading in txt file & converting it to a .csv file
tab = read.delim("referencefile.txt")
write.table(tab, file="referencefile.csv",sep=",",col.names=TRUE,row.names=FALSE)

#reading in only the recordID and taxonomy columns into taxdat.csv &
#removing all rows without genus-level ID
taxdat<- fread('referencefile.csv', select = c(3,10,12,14,16,20,22))
is.na(taxdat) <- taxdat==''
taxdat <- as.data.frame(taxdat)
taxdat <- taxdat[!is.na(taxdat[,6]), ]
  #6 above is for genus-level, 7 would remove all rows without a species-level ID, 
    #5 removes all without family-level ID, etc.

#make taxonomy file
#removing genus name from species name column
justspec<-str_split_fixed(taxdat$species_name, " ", 2)
justspec<- as.data.frame(justspec)
taxdat<-taxdat[,-7]
taxdat<-cbind(taxdat, justspec$V2)
names(taxdat)[7]<- "species_name"

#add taxonomy-level ID to each level of taxonomy
taxdat$phylum_name <- sub("^", "p__", taxdat$phylum_name)
taxdat$class_name <- sub("^", "c__", taxdat$class_name)
taxdat$order_name <- sub("^", "o__", taxdat$order_name)
taxdat$family_name <- sub("^", "f__", taxdat$family_name)
taxdat$genus_name <- sub("^", "g__", taxdat$genus_name)
taxdat$species_name <- sub("^", "s__", taxdat$species_name)

#combine all taxonomy into one row and remove all other rows
taxdat$taxonomy<- paste(taxdat$phylum_name, taxdat$class_name, taxdat$order_name, 
                        taxdat$family_name, taxdat$genus_name, taxdat$species_name, sep =";")
taxdat<- taxdat[,c(-2:-7)]
#convert data into matrix & write it as a text file named taxadata.txt
taxadata<- as.matrix(taxdat)
taxadata<- matrix(taxadata, ncol = ncol(taxdat), dimnames = NULL)
write.table(taxadata, "taxadata.txt", sep="\t", col.names = F, row.names = F)
#taxadata.txt can now be used as the .txt file in function 'assign_taxonomy.py' in QIIME

#make sequence file
#read in sequence data
seqdat<- fread('referencefile.csv', select = c(3,72))
is.na(seqdat) <- seqdat==''
seqdat <- as.data.frame(seqdat)
View(seqdat)
seqdat <- seqdat[!is.na(seqdat[,2]),]

#write the reference sequence fasta file as referencesequence.fasta
write.fasta(as.list(seqdat$nucleotides), seqdat$recordID, 'referencesequences.fasta')
#referencesequences.fasta can now be used as the .fasta file in function 'assign_taxonomy.py' in QIIME
