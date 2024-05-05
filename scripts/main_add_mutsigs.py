from scripts.pipeline_steps.step5 import *
import argparse
from scripts.globals import WT_THRESHOLD, MUT_THRESHOLD, MUTSIG_PROBABILITIES
import os
from pyensembl import EnsemblRelease

def main():
    
    
    parser = argparse.ArgumentParser()
    parser.add_argument("file_path", type=str, help="path to the csv file")
    
    args = parser.parse_args()

    VCF_FULL_PATH = args.file_path
    

    
    assembly = EnsemblRelease(75)

    df = apply_step_5(VCF_FULL_PATH, MUT_THRESHOLD, WT_THRESHOLD, MUTSIG_PROBABILITIES, assembly)


    csv_dir = os.path.dirname(VCF_FULL_PATH)
    VCF_ID = os.path.basename(VCF_FULL_PATH).split(".")[0]
    
    step5_csv = os.path.join(csv_dir, f"{VCF_ID}_step5.csv")
    df.to_csv(step5_csv, index=False)


if __name__ == '__main__':
    main()


