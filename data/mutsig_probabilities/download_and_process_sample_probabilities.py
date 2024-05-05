# this function is taken directly from https://github.com/AlexandrovLab/SigProfilerTopography/blob/master/SigProfilerTopography/Topography.py#L434

import os
import pandas as pd


def check_download_sample_probability_files():
    current_path = os.getcwd()
    os.makedirs(os.path.join(current_path, 'sample_probabilities'), exist_ok=True)
    sample_probability_files_path = os.path.join(current_path, 'sample_probabilities')

    if os.path.isabs(sample_probability_files_path):
        os.chdir(sample_probability_files_path)

        probability_files = ['COSMIC_DBS78_Decomposed_Mutation_Probabilities.txt',
                            'COSMIC_SBS96_Decomposed_Mutation_Probabilities.txt']

        for probability_filename in probability_files:
            probability_file_path = os.path.join(sample_probability_files_path, probability_filename)
            if not os.path.exists(probability_file_path):
                print(f'Does not exists: {probability_file_path}')
            try:
                print(
                    f'Downloading {probability_filename} under {sample_probability_files_path}'
                )

                # wget -c Continue getting a partially-downloaded file
                # wget -nc  If a file is downloaded more than once in the same directory, the local file will be clobbered, or overwritten
                # cmd="bash -c '" + 'wget -r -l1 -c -nc --no-parent -nd -P ' + chrombased_npy_path + ' ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerTopography/lib/nucleosome/chrbased/' + filename + "'"

                # -r When included, the wget will recursively traverse subdirectories in order to obtain all content.
                # -l1 Limit recursion depth to a specific number of levels, by setting the <#> variable to the desired number.
                # -c option to resume a download
                # -nc, --no-clobber If a file is downloaded more than once in the same directory, Wget's behavior depends on a few options, including -nc.  In certain cases, the local file will be clobbered, or overwritten, upon repeated download.  In other cases it will be preserved.
                # -np, --no-parent Do not ever ascend to the parent directory when retrieving recursively.  This is a useful option, since it guarantees that only the files below a certain hierarchy will be downloaded.
                # -nd, --no-directories When included, directories will not be created. All files captured in the wget will be copied directly in to the active directory
                cmd = "bash -c '" + 'wget -r -l1 -c -nc --no-parent -nd ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerTopography/sample_probability_files/' + probability_filename + "'"
                print(f"cmd: {cmd}")
                os.system(cmd)
            except Exception:
                print("The UCSD ftp site is not responding...")
    else:
        # It has to be an absolute path
        print(f'{sample_probability_files_path} is not an absolute path.')

def main():
    # check_download_sample_probability_files()
    
    df = pd.read_csv('sample_probabilities/COSMIC_SBS96_Decomposed_Mutation_Probabilities.txt', sep='\t')
    mutsig = df.iloc[:, 2:].idxmax(axis=1)

    # Append the mutsig column to the DataFrame
    df['mutsig'] = mutsig
    
    df = df.iloc[:, [0, 1, -1]]
    
    df.to_csv("sample_probabilities.csv", index=False)

if __name__ == '__main__':
    main()