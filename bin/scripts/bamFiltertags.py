#########################################################
# title: process_bam_with_barcodes
# author: Dawn
# date: 20250424
#########################################################

import argparse
import pandas as pd
import json
import pysam

def load_data(barcode_path, beads_path, best_barcode_path):
    """
    Load barcode and beads data from files and create a filtered barcode dictionary.
    """
    # Load barcode and beads data
    barcode_data = pd.read_csv(barcode_path, sep="\t", header=None)
    beads_data = pd.read_csv(beads_path, sep="\t", header=None)

    # Load best barcode mapping
    with open(best_barcode_path, "r") as f:
        best_barcode = json.load(f)

    # Filter barcodes
    filter_barcode = barcode_data[barcode_data[1].isin(beads_data[0].to_list())]
    filter_barcode = filter_barcode.set_index(0)
    filter_barcode_dict = filter_barcode.to_dict()[1]

    # Create best barcode dictionary
    best_barcode_dict = {}
    for k, v in best_barcode.items():
        if v in filter_barcode_dict:
            best_barcode_dict[k] = [v, filter_barcode_dict[v]]

    return best_barcode_dict

def process_bam(input_bam_path, output_bam_path, best_barcode_dict, sample_name):
    """
    Process the BAM file and update the barcode tags.
    """
    # Open input and output BAM files
    samfile = pysam.AlignmentFile(input_bam_path, "rb")
    outsam = pysam.AlignmentFile(output_bam_path, "wb", header=samfile.header)

    # Process each record in the BAM file
    for i in samfile:
        try:
            old_barcode = i.get_tag('CR')
            old_barcode = old_barcode[:10] + old_barcode[16:]
            if old_barcode in best_barcode_dict:
                CB = best_barcode_dict[old_barcode][0]
                DB = best_barcode_dict[old_barcode][1]
                i.set_tag("CB", CB)
                i.set_tag("DB", DB)
                i.set_tag("NB", "{}_{}".format(sample_name,DB))
                outsam.write(i)
        except KeyError:
            continue

    # Close files
    samfile.close()
    outsam.close()

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description="Process BAM file and update barcode tags.")
    parser.add_argument("--barcode_path", help="Path to barcodeTranslate.txt file")
    parser.add_argument("--beads_path", help="Path to beads_barcodes.txt file")
    parser.add_argument("--best_barcode_path", help="Path to best_barcode.json file")
    parser.add_argument("--input_bam_path", help="Path to input BAM file")
    parser.add_argument("--output_bam_path", help="Path to output BAM file")
    parser.add_argument("--sample_name", help="Sample name for the NB tag")

    # Parse arguments
    args = parser.parse_args()

    # Load data and create barcode dictionary
    best_barcode_dict = load_data(args.barcode_path, args.beads_path, args.best_barcode_path)

    # Process the BAM file
    process_bam(args.input_bam_path, args.output_bam_path, best_barcode_dict, args.sample_name)

if __name__ == "__main__":
    main()