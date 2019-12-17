#!/usr/bin/env python3
# Run Kraken2 on multiple files, choose database.


import argparse
import os
import subprocess
from tqdm import tqdm
from datetime import date
from time import time 


reports_folder = "/home/ubuntu/Data/Reports"
path_database = "/home/ubuntu/database/kraken2"

kraken2_db_choices = {
    "2015-10c"     : f"{path_database}/2015-10c",
    "2015-standard": f"{path_database}/2015-standard",
}

path_isolates = "/home/ubuntu/Disks/HDD500/ONT_isolates"
path_mocks = "/home/ubuntu/Disks/HDD500/ONT_Mock_Communities"
default_fastq_files = {
    "100000-SyntReads_20-BacGut": "/home/ubuntu/Disks/HDD1000/Segmentation/Test-Data/Synthetic_from_Genomes/2019-11-26_100000-SyntReads_20-BacGut.fastq"
    #"isolate_clostridioides_difficile": f"{path_isolates}/Clostridioides_difficile_SRR5457531.fastq",
    #"isolate_dengue_virus": f"{path_isolates}/DRR048282.fastq",
    #"isolate_e_coli": f"{path_isolates}/Escherichia_coli_SRR6118140.fastq",
    #"mock_10k" : f"{path_mocks}/Mock_10000-uniform-bacteria-l1000-q8.fastq",
#     "mock_100k": f"{path_mocks}/Mock_100000-bacteria-l1000-q10.fastq",
}

def run_kraken(path_file, db_choice):
    time_split = time()
        
    if "10c" in db_choice:
        for db_split in range(10):
            cmd = [
                "kraken2", "--threads", f"{cores}",
                "--db", f"{kraken2_db_choices[db_choice]}/{db_split}", 
                fastq_files[fastq], 
                "--output", f"{path_reports}/{db_choice}_{db_split}_{fastq}.kraken2_out",
                "--report", f"{path_reports}/{db_choice}_{db_split}_{fastq}.kraken2_report",
            ]
            print(f"report: {path_reports}/{db_choice}_{db_split}_{fastq}.kraken2_report")
            results = subprocess.check_output(cmd)
        
    else:
        cmd = [
            "kraken2", "--threads", f"{cores}",
            "--db", kraken2_db_choices[db_choice], 
            fastq_files[fastq], 
            "--output", f"{path_reports}/{db_choice}_{fastq}.kraken2_out",
            "--report", f"{path_reports}/{db_choice}_{fastq}.kraken2_report",
        ]
        print(cmd)
        results = subprocess.check_output(cmd)
    
    
    print(f"{fastq} done, {time()-time_split:.1f}s \n")


def scan_files(db_choice, cores, fastq_files=default_fastq_files):
    total = len(fastq_files)
    started = time()
    
    for i, fastq in enumerate(fastq_files):
        print(f"Processing {fastq} with database {db_choice}. File {i}/{total}.")
        path_reports = f"{reports_folder}/{date.today()}/{fastq}"
        if not os.path.exists(path_reports):
            os.makedirs(path_reports)
            
        time_split = time()
        
        if "10c" in db_choice:
            for db_split in range(10):
                cmd = [
                    "kraken2", "--threads", f"{cores}",
                    "--db", f"{kraken2_db_choices[db_choice]}/{db_split}", 
                    fastq_files[fastq], 
                    "--output", f"{path_reports}/{db_choice}_{db_split}_{fastq}.kraken2_out",
                    "--report", f"{path_reports}/{db_choice}_{db_split}_{fastq}.kraken2_report",
                ]
                print(f"report: {path_reports}/{db_choice}_{db_split}_{fastq}.kraken2_report")
                results = subprocess.check_output(cmd)
            
        else:
            cmd = [
                "kraken2", "--threads", f"{cores}",
                "--db", kraken2_db_choices[db_choice], 
                fastq_files[fastq], 
                "--output", f"{path_reports}/{db_choice}_{fastq}.kraken2_out",
                "--report", f"{path_reports}/{db_choice}_{fastq}.kraken2_report",
            ]
            print(cmd)
            results = subprocess.check_output(cmd)
        
        
        print(f"{fastq} done, {time()-time_split:.1f}s \n")
        
    print(f"Finished, total time = {time()-started:.1f}s")
    return


if __name__ == '__main__':
    # add defaults for -h  (formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser = argparse.ArgumentParser(description='Launch Kraken',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)  
    parser.add_argument('-d', '--database', help='Choose the Database to use', 
                        choices=['2019-mini', '2019-standard', '2015-standard', '2015-10c'], default='2019-mini') 
    parser.add_argument('-p', '--processes', help='Number of threads', type=int, default=2)
    args = parser.parse_args()

    scan_files(args.database, args.processes)
    
    
    
    