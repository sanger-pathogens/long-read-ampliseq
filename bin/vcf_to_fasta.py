#!/usr/bin/env python
from pysam import VariantFile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
import argparse, sys, os
import pandas as pd

#base script taken from nf-core bactmap and adjusted for long-read-ampliseq pipeline

class ParserWithErrors(argparse.ArgumentParser):
    def error(self, message):
        print('{0}\n\n'.format(message))
        self.print_help()
        sys.exit(2)

    def is_valid_file(self, parser, arg):
        if not os.path.isfile(arg):
            parser.error("The file %s does not exist!" % arg)
        else:
            return arg


def argparser():
    description = """
    A script to parse a filtered VCF and
    """
    parser = ParserWithErrors(description = description)
    parser.add_argument("-r", "--reference_file", required=True,
                        help="reference fasta file path",
                        type=lambda x: parser.is_valid_file(parser, x))
    parser.add_argument("-v", "--filtered_vcf_file", required=True,
                        help="filtered bcf file path",
                        type=lambda x: parser.is_valid_file(parser, x))
    parser.add_argument("-o", "--output_fasta_file", required=True,
                    help="file path to output fasta file")
    parser.add_argument("-i", "--fasta_id",
                    default="auto", help="fasta header ID")
    parser.add_argument("-m", "--multifasta",
                    default=False, help="multifasta output")
    parser.add_argument("-b", "--bed_file", type=lambda x: parser.is_valid_file(parser, x), required=True,
                        help="BED file (TSV) defining regions (<name>\t<start>\t<end>)" )
    parser.add_argument("-rr", "--replace_reference",
                    default=False, help="replace reference with gap")
    parser.add_argument("-g", "--gap_character", default= 'N', help="character to use for between amplicon gaps default - leave blank to only include variants (bit useless)")
    return parser

def calculate_reference_lengths(reference_file):
    reference_lengths = OrderedDict()
    for record in SeqIO.parse(reference_file, format = 'fasta'):
        reference_lengths[record.id] = len(record.seq)
    return reference_lengths

def get_variant_info(vcf_file):
    variant_info = OrderedDict()
    with VariantFile(vcf_file) as vcf:
        for record in vcf:
            # Extract position, reference allele, and alternate allele from the record
            position = record.pos
            ref_allele = record.ref
            alt_allele = record.alts[0]  # Assuming only one alternate allele is present
            variant_info[position] = (ref_allele, alt_allele)
    return variant_info

def variants_in_range(bed_range, variants):
    variants_list = []
    for key, value in variants.items():
        if key in bed_range:
            variants_list.append((key, value))
    return variants_list

def change_base_with_checks(sequence, position, base_change):
    """
    Change the base at a specified position in the sequence to a new base.

    Parameters:
    sequence (str): The biological sequence (e.g., DNA or RNA).
    position (int): The position (0-based index) where the base should be changed.
    new_base (str): The new base to replace the one at the specified position.

    Returns:
    str: The modified sequence with the base changed.
    """
    old_base, new_base = base_change
    if position < 0 or position >= len(sequence):
        raise ValueError(f"Position {position} is out of range: {len(sequence)}")
    if old_base != sequence[position] and sequence[position] != "N":
        raise ValueError(f"{sequence[position-5:position+6]} center base was expected to be {old_base}")
    if new_base not in 'ACGTUacgtu':
        raise ValueError("Variant not an accepted base.")
    
    
    # Change the base at the specified position
    new_sequence = sequence[:position] + new_base + sequence[position+1:]
    return new_sequence

def calculate_gaps_to_add(gap_start_position, gap_end_position, gapcharacter):
    return [gapcharacter] * (gap_end_position - gap_start_position)

def extract_sequences_from_bed_and_include_variants(reference_file, bed_file, variant_info, replace_reference, gap_character):
    #Reference genome as seqrecord
    ref_dict = SeqIO.to_dict(SeqIO.parse(reference_file, "fasta"))
    
    # Read the bed file in
    bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=['chromosome', 'start', 'end'])

    # Extract sequences
    extracted_sequences = []
    for index, row in bed_df.iterrows():
        chromosome = row['chromosome'] #not needed for this but good to ensure we are looking in the correct place? Multibed/multifasta ref later?
        start = row['start']
        end = row['end']
        variants = variants_in_range(range(start, end), variant_info)
        if variants:
            if replace_reference:
                gaps = calculate_gaps_to_add(start, end, gap_character)
                sequence = Seq("".join(gaps))
            else:
                sequence = ref_dict[chromosome].seq[start:end]
            for variant in variants:
                adjusted_start = variant[0] - start-1 #adjust to VCF not starting at 1
                sequence = change_base_with_checks(sequence, adjusted_start, variant[1])
            sequence_record = SeqRecord(sequence, id=f"{chromosome}_{start}_{end}", description="")
            extracted_sequences.append(sequence_record)
        else:
            if replace_reference:
                gaps = calculate_gaps_to_add(start, end, gap_character)
                sequence = Seq("".join(gaps))
            else:
                sequence = ref_dict[chromosome].seq[start:end]
            sequence_record = SeqRecord(sequence, id=f"{chromosome}_{start}_{end}", description="")
            extracted_sequences.append(sequence_record)
    
    return extracted_sequences

def write_sequence(filepath, multifasta, fasta_id, sequence_list):
    if multifasta:
        counter = 1
        with open(filepath, 'w') as output:
            for sequence in sequence_list:
                record = SeqRecord(Seq(sequence.seq), id = f"{fasta_id}_{counter}", description = '')
                SeqIO.write(record, output, "fasta")
                counter +=1
    else:
        with open(filepath, 'w') as output:
            master_record = ""
            for sequence in sequence_list:
                master_record += str(sequence.seq)
            record = SeqRecord(Seq(master_record), id = f"{fasta_id}_multi_locus", description = 'joined ref bed file regions with variants') 
            SeqIO.write(record, output, "fasta")
                



if __name__ == '__main__':
    parser = argparser()
    args = parser.parse_args()

    reference_lengths = calculate_reference_lengths(args.reference_file)
    variants = get_variant_info(args.filtered_vcf_file)

    extracted_sequences = extract_sequences_from_bed_and_include_variants(args.reference_file, args.bed_file, variants, args.replace_reference, args.gap_character)

    if args.fasta_id == "auto":
        fasta_id = os.path.basename(args.filtered_vcf_file).split('.')[0]
    else:
        fasta_id = args.fasta_id

    write_sequence(args.output_fasta_file, args.multifasta, fasta_id, extracted_sequences)