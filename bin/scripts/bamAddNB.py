import sys
import pysam
import json
import pandas as pd

def read_barcode_file(file_path, samplename):
    try:
        barcode = pd.read_csv(file_path, sep="\t", header=None)
        barcode.columns = ['barcode', 'cell']
        barcode['NB'] = ["{}_{}".format(samplename, i) for i in barcode['cell']]
        del barcode['cell']
        barcode = barcode.set_index("barcode")
        nb_tag_mapping = barcode.to_dict()['NB']
        return nb_tag_mapping
    except FileNotFoundError:
        print(f"Error: Barcode file '{file_path}' not found.")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Error: Barcode file '{file_path}' is empty.")
        sys.exit(1)
    except pd.errors.ParserError as e:
        print(f"Error: Invalid format in barcode file '{file_path}': {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading barcode file: {e}")
        sys.exit(1)

def process_bam_file(input_bam_path, nb_tag_mapping, output_bam_path):
    try:
        samfile = pysam.AlignmentFile(input_bam_path, "rb")
        outsam = pysam.AlignmentFile(output_bam_path, "wb", header=samfile.header)
        
        barcode_list = list(nb_tag_mapping.keys())
        for i in samfile:
            try:
                if i.has_tag('CB'):
                    old_barcode = i.get_tag('CB')
                    if old_barcode in barcode_list:
                        new_barcode = nb_tag_mapping[old_barcode]
                        i.set_tag("NB", new_barcode)
                outsam.write(i)
            except KeyError:
                continue
    except pysam.SamtoolsError as e:
        print(f"Error processing SAM/BAM file: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error processing SAM/BAM file: {e}")
        sys.exit(1)
    finally:
        outsam.close()

def main():
    if len(sys.argv) < 5:
        print("Usage: python script.py <samplename> <barcode_valid> <input.bam> <output.bam>")
        sys.exit(1)

    samplename = sys.argv[1]
    barcode_valid = sys.argv[2]
    input_bam_path = sys.argv[3]
    output_bam_path = sys.argv[4]

    nb_tag_mapping = read_barcode_file(barcode_valid, samplename)
    process_bam_file(input_bam_path, nb_tag_mapping, output_bam_path)

if __name__ == "__main__":
    main()