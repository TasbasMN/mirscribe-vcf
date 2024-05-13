import csv
import json

def generate_csv_and_json_from_vcf(vcf_file_path, vcf_path, json_path):
    try:
        # Open the VCF file for reading
        with open(vcf_file_path, 'r') as vcf_file:
            # Prepare to write to CSV and collect data for JSON
            with open(vcf_path, 'w', newline='') as file1:
                csv_writer = csv.writer(file1)
                variant_dict = {}

                # Process each line in the VCF file
                for line in vcf_file:
                    if line.startswith('#'):
                        continue  # Skip header lines

                    fields = line.strip().split('\t')
                    chr, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]

                    # Write data to CSV file
                    csv_writer.writerow([chr, pos, ref, alt])

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
                    # Create unique_id and update variant_dict
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