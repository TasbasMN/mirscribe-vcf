import os
import csv
import json

def generate_csv_and_json_from_vcf(vcf_file_path, output_dir):
    try:
        # Create the output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)

        # Construct the output file paths
        vcf_filename = os.path.basename(vcf_file_path)
        vcf_name, _ = os.path.splitext(vcf_filename)
        vcf_path = os.path.join(output_dir, f"{vcf_name}.csv")
        json_path = os.path.join(output_dir, f"{vcf_name}.json")

        # Open the VCF file for reading
        with open(vcf_file_path, 'r') as vcf_file:
            # Prepare to write to CSV and collect data for JSON
            with open(vcf_path, 'w', newline='') as file1:
                csv_writer = csv.writer(file1, delimiter='\t')  # Use tab delimiter
                variant_dict = {}

                # Process each line in the VCF file
                for line in vcf_file:
                    if line.startswith('#'):
                        continue  # Skip header lines

                    fields = line.strip().split('\t')
                    chr, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]

                    # Write data to CSV file with file name as unique ID
                    csv_writer.writerow([chr, pos, vcf_name, ref, alt])

                    # Process INFO field for JSON file
                    info = fields[7]
                    variant_classification = next(
                        (
                            item.split('=')[1]
                            for item in info.split(';')
                            if item.startswith('Variant_Classification=')
                        ),
                        'Unknown',
                    )
                    unique_id = '_'.join([chr, pos, ref, alt])
                    if variant_classification in variant_dict:
                        variant_dict[variant_classification].append(unique_id)
                    else:
                        variant_dict[variant_classification] = [unique_id]

                # Export the dictionary to a JSON file
                with open(json_path, 'w') as file2:
                    json.dump(variant_dict, file2, indent=4)

    except Exception as e:
        print(f"An error occurred: {e}")

    return vcf_path