#########################################################
# title:filter cell bam
# author:dawang
# date:20241023
#########################################################


import sys
import pysam
import json
import pandas as pd
def main():
    if len(sys.argv) < 4:
        print("Usage: python bamFilterNB.py <input.bam> <beads_barcodes.txt> <sample><output.bam>")
        sys.exit(1)

    input_bam_path = sys.argv[1]
    beads_barcodes_path = sys.argv[2]
    sample_name = sys.argv[3]
    output_bam_path = sys.argv[4]

    cell_list = pd.read_csv(beads_barcodes_path,header=None)[0].to_list()
    cell_list=['{}_{}'.format(sample_name,i) for i in cell_list]

    try:
        samfile = pysam.AlignmentFile(input_bam_path, "rb")
        outsam = pysam.AlignmentFile(output_bam_path, "wb", header=samfile.header)
        
        for i in samfile:
            try:
                if i.has_tag('NB'):
                    old_barcode = i.get_tag('NB')
                    if old_barcode in cell_list:                 
                        outsam.write(i)
            except KeyError:
                continue
    except Exception as e:
        print(f"Error processing SAM/BAM file: {e}")
        sys.exit(1)
    finally:
        outsam.close()

if __name__ == "__main__":
    main()