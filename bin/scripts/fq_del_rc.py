import argparse

## CR: 16:41 
## UR: 0:10

def modify_fastq(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        while True:
            # Read four lines at a time (FASTQ format)
            identifier = infile.readline().strip()
            sequence = infile.readline().strip()
            plus = infile.readline().strip()
            quality = infile.readline().strip()

            # Break the loop if we've reached the end of the file
            if not identifier:
                break

            # Modify the identifier
            identifier = identifier.split(" ")[0]
            new_identifier = f"{identifier}"
            
            # Write the modified lines to the output file
            outfile.write(f"{new_identifier}\n{sequence}\n{plus}\n{quality}\n")

def main():
    parser = argparse.ArgumentParser(description="Modify FASTQ files by extracting a subsequence and modifying the identifier.")
    parser.add_argument('--input_file', type=str, required=True, help='Path to the input FASTQ file')
    parser.add_argument('--output_file', type=str, required=True, help='Path to the output FASTQ file')

    args = parser.parse_args()

    modify_fastq(args.input_file, args.output_file)

if __name__ == "__main__":
    main()