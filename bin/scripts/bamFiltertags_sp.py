#########################################################
# title: process_bam_with_barcodes
# author: Dawn
# date: 20250507
#########################################################

import argparse
import json
import pysam
import logging
import os

# 配置日志记录
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_best_barcode_dict(best_barcode_path):
    """
    Load the best barcode dictionary from a JSON file.
    :param best_barcode_path: Path to the best_barcode.json file
    :return: Dictionary of best barcodes
    """
    try:
        with open(best_barcode_path, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        logging.error(f"File not found: {best_barcode_path}")
        raise
    except json.JSONDecodeError:
        logging.error(f"Invalid JSON format in file: {best_barcode_path}")
        raise

def process_bam(input_bam_path, output_bam_path, best_barcode_dict, sample_name):
    """
    Process the BAM file and update the barcode tags.
    :param input_bam_path: Path to the input BAM file
    :param output_bam_path: Path to the output BAM file
    :param best_barcode_dict: Dictionary of best barcodes
    :param sample_name: Sample name for the NB tag
    """
    if not os.path.exists(input_bam_path):
        logging.error(f"Input BAM file not found: {input_bam_path}")
        return

    try:
        # Open input and output BAM files
        samfile = pysam.AlignmentFile(input_bam_path, "rb")
        outsam = pysam.AlignmentFile(output_bam_path, "wb", header=samfile.header)

        # Process each record in the BAM file
        for i in samfile:
            try:
                old_barcode = i.get_tag('CR')
                if old_barcode in best_barcode_dict:
                    CB = best_barcode_dict[old_barcode][0]
                    DB = best_barcode_dict[old_barcode][1]
                    i.set_tag("CB", CB)
                    i.set_tag("DB", DB)
                    i.set_tag("NB", "{}_{}".format(sample_name,DB))
                    outsam.write(i)
            except KeyError:
                logging.warning(f"Tag 'CR' not found in record: {i}")
            except Exception as e:
                logging.error(f"Error processing record: {i}. Error: {e}")

        # Close files
        samfile.close()
        outsam.close()
    except Exception as e:
        logging.error(f"Error processing BAM file: {e}")

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description="Process BAM file and update barcode tags.")
    parser.add_argument("--best_barcode_path", required=True, help="Path to best_barcode.json file")
    parser.add_argument("--input_bam_path", required=True, help="Path to input BAM file")
    parser.add_argument("--output_bam_path", required=True, help="Path to output BAM file")
    parser.add_argument("--sample_name", required=True, help="Sample name for the NB tag")

    # Parse arguments
    args = parser.parse_args()

    # Load best barcode dictionary
    best_barcode_dict = load_best_barcode_dict(args.best_barcode_path)

    # Process the BAM file
    process_bam(args.input_bam_path, args.output_bam_path, best_barcode_dict, args.sample_name)

if __name__ == "__main__":
    main()