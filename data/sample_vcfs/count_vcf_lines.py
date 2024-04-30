import os
import sys

def count_lines_in_vcf_files(directory):
    report_file_path = os.path.join(directory, "vcf_file_line_counts.txt")
    vcf_files = [f for f in os.listdir(directory) if f.endswith(".vcf")]
    vcf_files.sort()  # Sort the list of VCF files alphabetically

    with open(report_file_path, "w") as report_file:
        for filename in vcf_files:
            file_path = os.path.join(directory, filename)
            with open(file_path, "r") as file:
                line_count = sum(1 for _ in file)
            report_file.write(f"{filename} line count: {line_count}\n")
            print(f"Line count for {filename} written to {report_file_path}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <directory_path>")
        sys.exit(1)

    directory_path = sys.argv[1]
    count_lines_in_vcf_files(directory_path)
