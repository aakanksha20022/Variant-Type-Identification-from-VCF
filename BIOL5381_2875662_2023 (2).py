#Importing all the relevant modules
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
import gffutils
import vcf
import os
import matplotlib.pyplot as plt
import logging
import argparse

#Creating a logger instance
logger = logging.getLogger()

#setting the logging level to INFO, which will capture messages with severity INFO and above
logger.setLevel(logging.INFO)

#Creating a Filehandler to log messages to a file named 'Log.txt'
ch= logging.FileHandler('2875662_log.txt')

#Setting the log message format using a Formatter
ch.setFormatter(logging.Formatter('%(levelname)s-%(asctime)s-%(message)s'))

#Adding the fileHandler to the logger to start logging to the specified file
logger.addHandler(ch)

#Creating an ArgumentParser for handling command line arguments
parser = argparse.ArgumentParser(description='Process VCF, GFF, and FASTA files.')

#Defining the command line arguments for VCF, GFF and FASTA files
parser.add_argument('--vcf', required=True, help='Path to the VCF file.')
parser.add_argument('--gff', required=True, help='Path to the GFF file.')
parser.add_argument('--fasta', required=True, help='Path to the FASTA file.')

#Parsing the command line arguments with the parse_args method
args = parser.parse_args()

#Extracting the file paths from command line arguments
vcf_file = args.vcf
gff_file = args.gff
fasta_file = args.fasta

#Generating the default name for the database file based on the GFF file
db_file= os.path.splitext(args.gff)[0] + ".db"

#Logging file paths for referance
logger.info(f'VCF file: {args.vcf}\n')
logger.info(f'GFF file: {args.gff}\n')
logger.info(f'FASTA file: {args.fasta}\n')


# Initializing counters for each variant type
non_coding_count = 0
synonymous_count = 0
non_synonymous_count = 0
qual_count = 0

#Checking if database file already exists
if not os.path.isfile(db_file):
    try:
    #If the database file does not exist, creating a new database
        print('Creating database...\n')
        db = gffutils.create_db(args.gff, dbfn= db_file, force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
    except ValueError:
        logger.error(f'Cannot create database {db_file}.\n')
        raise SystemExit(1)

else:
   #If the database file exists, connecting to the existing database
    try:
        print('Connecting to existing database..\n')
        db = gffutils.FeatureDB(db_file, keep_order=True)
    except ValueError:
        logger.error(f'Cannot read existing database {db_file}.\n')
        raise SystemExit(1)

#Creating a reader object of the vcf file
try:
    vcfReader= vcf.Reader(filename= args.vcf)
except FileNotFoundError:
    logger.error(f'File {vcf_file} cannot be opened for reading. Please check and try again.\n')
    raise SystemExit(1)

#Creating a function to extract the concatenated cds sequence within a transcript
def process_cds_sequence(cds, fasta_file, seq, loc, feature):
    #Concatenating the sequence of the cds to the empty sequence string
    seq += cds.sequence(fasta_file, use_strand=True)
    #Checking if the current cds is equal to the feature containing the variant
    if cds == feature:
        #If equal, returning the updated sequence, location and True to indicate the loop should break
        return seq, loc, True  # Indicates the loop should break
    else:
        #If not eual, updating the location with the sequence of the coding region
        loc += cds.sequence(fasta_file, use_strand=True)
        return seq, loc, False  # Indicates the loop should continue

#Creating a function to calculate the protein location in the translated sequences
#Getting the ref and alt amino acids  
def calc_protein_location(var_loc, prot_seq_ref, prot_seq_mut):
    
    #Calculate the protein location by diving the variant location by (3+1)
    #Adding 1 ensures that the protein location is reported in a one-based index
    pro_loc = int((var_loc) // 3 + 1)

    #Checking if the calculated protein location is within the reference and mutant protein sequences
    if 0 <= pro_loc <= len(prot_seq_ref) and 0 <= pro_loc <= len(prot_seq_mut):

        #If it is, extracting the ref and alt amino acids at the calculated position
        #Subtracting 1 to revert back to zero based indexing to access the correct amino acid
        ref_aa = prot_seq_ref[pro_loc - 1]
        alt_aa = prot_seq_mut[pro_loc - 1]

        #Returning the protein location, ref and alt amino acids
        return pro_loc, ref_aa, alt_aa
    else:
        #Returning None for pro_loc, ref_aa, and alt_aa if conditions are not met
        return None, None, None  
    

    
try:
    # Opening a file for writing the table
    with open('2875662_output.tsv', 'w') as file:
        # Writing the header
        file.write("Chrom\tPos\tRef\tAlt\tType\tTranscript\tProtein_Location\tRef_AA\tAlt_AA\n")
        
        #Iterating through records in vcfReader
        for record in vcfReader:

            #Initializing the protein location, transcript, ref and alt amino acid values to "NA"
            #For variants which are not present in the coding regions
            trans= 'NA'
            pro_loc= 'NA'
            ref_aa= 'NA'
            alt_aa= 'NA'

            #Checking if the quality of the variants is less than or equal to 20 and if the condition is true, incrementing the qual_count variable by 1
            if record.QUAL <= 20:
                qual_count +=1
            
            #Checking if the quality of the variants is greater than 20
            if record.QUAL > 20:
                
                #Retrieving a list of coding regions for the current VCF record position
                coding_regions = list(db.region(seqid=record.CHROM, start=record.POS, end=record.POS, featuretype='CDS'))

                #Setting the default variant type to Non-coding
                variant_type= 'Non-coding'

                #Incrementing the count of non-coding variables
                non_coding_count +=1

                #Checking if there are coding regions for the current position
                if coding_regions:
                    for feature in coding_regions:
                        for transcript in db.parents(feature, featuretype='mRNA'):
                            trans= transcript.id
                            
                            seq= ''  #Initialize an empty sequence
                            loc= ''  #Initialize an empty location
                            if transcript.strand == '+':
                                #Processing coding regions for the positive strand by calling the get_cds_seq function
                                for cds in db.children(transcript, order_by= 'start', featuretype= 'CDS'):
                                    seq, loc, should_break = process_cds_sequence(cds, args.fasta, seq, loc, feature)

                                    if should_break:
                                        break
                                #Calculating the variant location from cds start position
                                start_pos = feature.start
                                var_loc = len(loc) + (record.POS - start_pos)
                                
                                
                                #Creating Seq and MutableSeq objects for reference and mutated sequences
                                ref_seq= Seq(seq)
                                mut_seq= MutableSeq(seq)

                                #Checking if the variant location is within the processed sequence
                                if 0 <= var_loc < len(mut_seq):
                                    for alt in record.ALT:
                                        mut_seq[var_loc] = str(alt)
                                

                                    # Translating both sequences
                                    prot_seq_ref = ref_seq.translate()
                                    prot_seq_mut = mut_seq.translate()

                                    #Calculating protein location and reference/alternate amino acids by calling the calc_protein_location function
                                    pro_loc, ref_aa, alt_aa = calc_protein_location(var_loc, prot_seq_ref, prot_seq_mut)
                        
                                    
                            elif transcript.strand == '-':

                                #Processing the coding regions for the negative strand
                                for cds in db.children(transcript, order_by='start', featuretype='CDS', reverse= True):
                                    seq, loc, should_break = process_cds_sequence(cds, args.fasta, seq, loc, feature)

                                    if should_break:
                                        break
                                #Calculating the variant location from cds end position    
                                end_pos = feature.end
                                var_loc = len(loc) + (end_pos-record.POS)
                                
                                ref_seq= Seq(seq)
                                mut_seq= MutableSeq(seq)

                                #Checking if the variant location is within the processed sequence
                                if 0 <= var_loc < len(mut_seq):

                                    #Converting ALT nucleotide based on complementary strand
                                    for alt in record.ALT:
                                        if alt == 'G':
                                            alt = 'c'
                                        elif alt == 'A':
                                            alt  = 't'
                                        elif alt == 'C':
                                            alt = 'g'
                                        elif alt == 'T':
                                            alt = 'a'
                                        mut_seq[var_loc] = str(alt).upper()
                                        

                                        # Translating both sequences
                                        prot_seq_ref = ref_seq.translate()
                                        
                                        prot_seq_mut = mut_seq.translate()

                                        #Calculating protein location and ref/alt amino acids by calling the calc_protein_location function
                                        pro_loc, ref_aa, alt_aa = calc_protein_location(var_loc, prot_seq_ref, prot_seq_mut)

                                        
                            #Checking if the reference and mutated protein sequences are equal
                            #Categorizing them as synonymous if they are equal and non-synonymous if they are unequal
                            if prot_seq_ref == prot_seq_mut:
                                variant_type = 'Synonymous'
                                synonymous_count +=1
                            else:
                                variant_type = 'Non-synonymous'
                                non_synonymous_count +=1

                        #Writing the results to the output file
                        file.write(f"{record.CHROM}\t{record.POS}\t{record.REF}\t{record.ALT[0]}\t{variant_type}\t{trans}\t{pro_loc}\t{ref_aa}\t{alt_aa}\n")
                else:
                    #Writing default values to the output file if no coding regions are present
                    file.write(f"{record.CHROM}\t{record.POS}\t{record.REF}\t{record.ALT[0]}\t{variant_type}\tNA\tNA\tNA\tNA\n")
except FileNotFoundError:
    logger.error(f'File {file} cannot be opened for writing. Please check and try again.\n')
    raise SystemExit(1)

#Logging information about number of variants with quality <=20
logger.info(f'Number of variants with quality <= 20 is {qual_count}\n')
logger.info(f'Output table file created: 2875662_output.tsv\n')

#Counts for different variant types           
non_coding_proportion = non_coding_count 
synonymous_proportion = synonymous_count 
non_synonymous_proportion = non_synonymous_count 

# Bar plot
labels = ['Non-coding', 'Synonymous', 'Non-synonymous']
counts = [non_coding_proportion, synonymous_proportion, non_synonymous_proportion]

#Creating a bar plot using Matplotlib
plt.bar(labels, counts, color=['blue', 'green', 'red'])
plt.title('count of Variants by Type')
plt.xlabel('Variant Type')
plt.ylabel('count')
plt.show()

#Saving the plot as an image file names '2875662_barplot.png'
plt.savefig('2875662_barplot.png')
logger.info('Output plot file created: 2875662_barplot.png\n')



                       