import argparse
import os

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Process a VCF file in chunks using concurrent futures.')
    parser.add_argument('file_path', default="data/sample_vcfs/sample.vcf",
                        type=str, help='Path to the VCF file')
    parser.add_argument("-c", '--chunksize', default=200,
                        type=int, help='Number of lines to process per chunk')
    parser.add_argument("-o", '--output_dir', type=str,
                        default='./results', help='Path to the output directory')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose logging')
    parser.add_argument('-w', '--workers', default=os.cpu_count(),
                        type=int, help='Number of concurrent workers')
    parser.add_argument('--skip-rnaduplex', action='store_true',
                        help='Skip RNAduplex analysis')
    parser.add_argument('-t', '--threshold', default=0.2, type=float, 
                        help='Threshold for filtering out pairs that have less prediction difference than the threshold')
                        
                        
    return parser.parse_args()


args = parse_arguments()

VCF_FULL_PATH = args.file_path
CHUNKSIZE = args.chunksize
VERBOSE = args.verbose
WORKERS = args.workers
SKIP_RNADUPLEX = args.skip_rnaduplex
FILTER_THRESHOLD = args.threshold


VCF_ID = os.path.basename(VCF_FULL_PATH).split(".")[0]
OUTPUT_DIR = os.path.join(args.output_dir, f"{VCF_ID}_{CHUNKSIZE}")
