import argparse
import os


parser = argparse.ArgumentParser(
    description='Process a VCF file in chunks using concurrent futures.')
parser.add_argument('file_path', type=str, help='Path to the VCF file')
parser.add_argument("-c", '--chunksize', default=100,
                    type=int, help='Number of lines to process per chunk')
parser.add_argument("-o", '--output_dir', type=str,
                    default='./results', help='Path to the output directory')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='Enable verbose logging')
parser.add_argument('-w', '--workers', default=os.cpu_count(),
                    type=int, help='Number of concurrent workers')
parser.add_argument('--skip-rnaduplex', action='store_true',
                    help='Skip RNAduplex analysis')
args = parser.parse_args()

VCF_FULL_PATH = args.file_path
VCF_ID = os.path.basename(VCF_FULL_PATH).split(".")[0]
VERBOSE = args.verbose
CHUNKSIZE = args.chunksize

WORKERS = args.workers

SKIP_RNADUPLEX = args.skip_rnaduplex

OUTPUT_DIR = os.path.join(args.output_dir, f"{VCF_ID}_{CHUNKSIZE}")

# Set other configuration variables based on command-line arguments
