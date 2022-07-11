import subprocess
def to_r_data_processing(doc_name_1, doc_name_2, doc_name_3, location):
    command = 'Rscript'
    path2script = '/home/anna/rstudio-2022.02 (1).3-492-amd64-debian/completion_of_PCA_and_Dendro.R'
    args= [location+'/'+doc_name_1, location+'/'+doc_name_2, location+'/'+doc_name_3]
    cmd =[command, path2script] +args
    subprocess.check_output(cmd, universal_newlines=True)
    #subprocess.call(['RStudio.sh', '2022_07_11_pca_dendro.R', doc_name_1, doc_name_2, doc_name_3])
