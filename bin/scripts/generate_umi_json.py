import pandas as pd
import json
import argparse
import os
import logging



def _numerical_to_dna(numerical_str, length):
    """
    Convert a hexadecimal numerical string to a DNA sequence.
    
    Parameters:
    numerical_str (str): The hexadecimal numerical string.
    length (int): The length of the DNA sequence to generate.
    
    Returns:
    str: The DNA sequence.
    """
    NT_COMP = {'0': 'A', '1': 'C', '2': 'G', '3': 'T'}
    numerical = int(numerical_str, 16)  # Convert hex string to integer
    dna = []

    for _ in range(length):
        numerical, remainder = divmod(numerical, 4)
        dna.append(NT_COMP[str(remainder)])
    dna = dna[::-1]  # Reverse the list to get the correct order
    return ''.join(dna)


def process_data(path, data_type):
    """
    Read and process the input CSV file based on data type (LRS or SRS).
    
    Parameters:
    path (str): Path to the input CSV file.
    data_type (str): Type of data ('lrs' or 'srs').
    
    Returns:
    pd.DataFrame: Processed DataFrame.
    """
    # Read the CSV file, assuming columns are separated by spaces and there is no header
    data = pd.read_csv(path, sep=" ", header=None)
    data = data.rename(columns={0: "id", 1: "umi", 2: "gene"})
    data = data.drop_duplicates()
    
    if data_type == 'srs':
        data['barcode'] = [_numerical_to_dna(i, 20) for i in data['id']]
    if data_type == 'lrs':
        data['barcode'] = [i[:10]+i[-10:] for i in data['id']]
        
    data['gb'] = data['gene'] +"_"+ data['barcode']
        
    return data


def save_data(data, output_dir, data_type):
    """
    Save aggregated data to CSV and JSON files.
    
    Parameters:
    data_info (pd.DataFrame): Aggregated DataFrame.
    output_dir (str): Output directory.
    data_type (str): Type of data ('lrs' or 'srs').
    """
    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    data_info = data.groupby("gb").agg({"umi": list})
    # Save aggregated data to JSON
    json_file = f"{output_dir}/{data_type}_gene_barcode_umi.json"
    d_dict = data_info.to_dict()['umi']
    with open(json_file, 'w') as f:
        json.dump(d_dict, f)
    logging.info(f"Output files saved to {output_dir} for {data_type.upper()} data.")


def main(lrs_path, srs_path, output_dir):
    """
    Main function to process LRS and SRS data and generate output files.
    
    Parameters:
    lrs_path (str): Path to the LRS input CSV file.
    srs_path (str): Path to the SRS input CSV file.
    output_dir (str): Output directory for generated files.
    """
    # Process and save LRS data
    lrs_data = process_data(lrs_path, 'lrs')
    save_data(lrs_data, output_dir, 'lrs')
    
    # Process and save SRS data
    srs_data = process_data(srs_path, 'srs')
    save_data(lrs_data, output_dir, 'srs')


if __name__ == "__main__":
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description="Generate JSON and CSV files from LRS and SRS CSV inputs.")
    parser.add_argument("--lrs_path", help="The path to the LRS input CSV file.")
    parser.add_argument("--srs_path", help="The path to the SRS input CSV file.")
    parser.add_argument("--outdir", help="The directory where the output files will be saved.")
    
    # Parse command line arguments
    args = parser.parse_args()
    
    # Call the main function
    main(args.lrs_path, args.srs_path, args.outdir)