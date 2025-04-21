#########################################################
# title:addCB
# author:dawang
# date:20241023
#########################################################


import sys
import pysam
import json

def main():
    if len(sys.argv) < 4:
        print("Usage: python script.py <input.bam> <best_barcode.json> <output.bam>")
        sys.exit(1)

    input_bam_path = sys.argv[1]
    best_barcode_path = sys.argv[2]
    output_bam_path = sys.argv[3]

    try:
        with open(best_barcode_path, "r") as f:
            best_barcode = json.load(f)
    except Exception as e:
        print(f"Error reading best barcode file: {e}")
        sys.exit(1)

    barcode_list = list(best_barcode.keys())
    try:
        samfile = pysam.AlignmentFile(input_bam_path, "rb")
        outsam = pysam.AlignmentFile(output_bam_path, "wb", header=samfile.header)
        
        for i in samfile:
            try:
                if i.has_tag('CR'):
                    old_barcode = i.get_tag('CR')
                    old_barcode =old_barcode[:10]+old_barcode[16:]
                    if old_barcode in barcode_list:
                        new_barcode = best_barcode[old_barcode]
                        i.set_tag("CB", new_barcode)
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