# User input. Choose synbio file, choose domain of synbio
import subprocess
import os
import shutil

print("Welcome to EcoGeno! Enter your pathway for the synbio FASTA file: ")
sb_path= input()
print("Enter where you want results saved: ")
desired_location=input()
finding_name=sb_path.split("/")
group_num=len(finding_name)
sb_name=finding_name[group_num]
matches=sb_name+"_matches"

#DIAMOND run
blastp=['diamond', 'blastp', '-d', '/home/anna/PycharmProjects/pythonProject/DIAMOND_matches/reference.dmnd', '-q', sb_path,
        '-o', matches, '--max-target-seqs', '1', '--outfmt', '6']  # completes DIAMOND search by using full path
subprocess.run(blastp)
shutil.move(os.path.abspath(matches), desired_location)

#Matrix Comparison 