Variant-Type-Identification-from-VCF
---
The provided script is designed to analyze genetic variants from VCF, GFF, and FASTA files, generating detailed output and visualizations. The script utilizes a combination of logging and 
command-line argument handling for efficient execution.

The script leverages the argparse module for handling command-line arguments. Users are required to provide paths to the VCF, GFF, and FASTA files using the following command-line arguments:
1. 	--vcf: Path to the VCF file.
2. 	--gff: Path to the GFF file.
3. 	--fasta: Path to the FASTA file.

Example usage:
python3 BIOL5381_2875662_2023.py --vcf filename.vcf --gff filename.gff --fasta filename.fasta

Output Generation
1. The '2875662_output.tsv' serves as a comprehensive summary of processed genetic variants, providing a structured tabular representation of key attributes. Each row in the table corresponds to a 
   distinct genetic variant and includes essential information such as the chromosome identifier ('Chrom'), position of the variant on the chromosome ('Pos'), reference allele ('Ref'), alternative 
   allele(s) ('Alt'), variant type ('Type') categorized as non-coding, synonymous, or non-synonymous, associated transcript identifier ('TranscriptID'), protein location within the translated 
   amino acid sequence ('Protein_Location'), reference amino acid ('Ref_AA'), and alternative amino acid(s) resulting from the variant ('Alt_AA'). 
2. File-related errors and the count of variants with quality scores less than or equal to 20 is accounted for in the log file names ‘2875662_log.txt’. Total counts of non-coding, synonymous, and 
   non-synonymous variants are calculated.
3. A bar plot ('2875662_barplot.png') is generated to visually represent the variant counts by type.

