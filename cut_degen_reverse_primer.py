#Put this file in the same folder as your .fasta file and open this file from that folder
#Cut Primers for Q30 sequences 
#requires Python2 & Biopython

from Bio import SeqUtils
from Bio import SeqIO
import sys



#for records and adapter, should be sys.argv[1 and 2]
fastafile =  "seqs_w_for_removed.fasta" #This is the input .fasta file from 'cut_degen_forward_primer.py' 
adapter =  "ATGATGATGATGATGATGATGATG" #This is the reverse primer as a string (degenerate primers can be used)
keepreads = False #True or False, this will determine whether or not reads are kept. (Keep to False unless you want to keep potentially errenous sequences)
removeadapters = True #True or False. If this is True, the adapters will be removed. (Keep True as adapters need removed) 
end_defn = 3 #If 5, the primer is removed from the 5' end of the sequence (forward primer). If 3, then it is removed from the 3' end of the sequence (reverse compliment of reverse primer).
adapter_name = 'Reverse Primer Name' #This is the name of the adapter that you can put into the output text file.

#make changes above as needed for specific primer and preferences

keepreads = str(keepreads)
removeadapters = str(removeadapters)
fastafile=str(fastafile)
end_defn = str(end_defn)

fh = open(fastafile, mode='r+')
len_adapter = len(adapter)
count_adapter_found = 0
count_adapter_not_found = 0
total_seq_count = 0

parsed = SeqIO.parse(fh, format="fasta")

output_fh_name = "seqs_w_primers_removed"

if fastafile=="test3prime.fasta":
    output_fh_name="seqs_w_primers_removed.fasta"

output_fh = open(output_fh_name, mode='w+')

output_text_name = "info_w_primers_removed.txt"
if fastafile=="test3prime.fasta":
    output_text_name="info_w_primers_removed.txt"
    
output_text_fh = open(output_text_name, mode='w+')


for record in parsed:
    try:
        sequence = str(record.seq)
        search = SeqUtils.nt_search(sequence, adapter) #This will search the
        index = int(search[1]) #If it finds the adapter, is the starting index from which it was found.
        adapter_start = index
        adapter_end = index+len_adapter
        count_adapter_found +=1
        total_seq_count+=1
        if removeadapters == "True": #if the value is true, it removes the adapters from the sequences.
            if end_defn=="5":
                record = record[adapter_end:] #If a 5' adapter, you remove adapter from beginning
            elif end_defn=="3":
                record = record[:adapter_start] #If it is a 3' adapter, you remove the adapter at the end
        elif removeadapters == "False": #if the value is false, it does not remove the adapters from the sequences.
            record = record
        SeqIO.write(record, output_fh, format="fasta") #No matter what, write the reads.
    except IndexError:
        count_adapter_not_found+=1
        total_seq_count+=1
        record = record
        if keepreads=="True":
            SeqIO.write(record, output_fh, format="fasta")
        elif keepreads=="False":
            pass
        else:
            pass

output_fh.close()

percent_cut = 100*(float(count_adapter_found)/float(total_seq_count))


output_text_fh.write("The total number of sequences that were analyzed was %i.\n\n"%total_seq_count)
output_text_fh.write("Adapter was found and removed for %i sequences (%i%% of total).\n\n"%(count_adapter_found, percent_cut))

if keepreads =="True":
    output_text_fh.write("Sequences that did not contain the adapter were kept.\n\n")
elif keepreads=="False":
    output_text_fh.write("Sequences that did not contain the adapter were removed from the dataset.\n\n")
if removeadapters=="True":
    output_text_fh.write("The adapters were removed from the dataset.\n\n")
elif removeadapters=="False":
    output_text_fh.write("The adapters were not removed from the dataset.\n\n")
if end_defn=="5":
    output_text_fh.write("Adapters were removed from the 5\' end.\n\n")
elif end_defn=="3":
    output_text_fh.write("Adapters were removed from the 3\'end.\n\n")

output_text_fh.write("The name of the adapter that was removed was named %s, and had the sequence %s.\n\n"%(adapter_name,adapter))
output_text_fh.close()
