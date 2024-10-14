import streamlit as st
import numpy as np
import pandas as pd
import os
import shutil
import os.path
import numpy as np 
import pandas as pd 
import numpy as np 
red = "\033[91m"
reset_color = "\033[0m"
import seaborn as sns
import matplotlib.pyplot as plt
import tempfile
from EnCen_Functions import EC_extract, tsv_to_fasta, diamond_impl, genome_extractor_ref, genome_extractor_syn, genome_to_genome_diffcomp, read_in_binary_matrix, calculating_distance, pass_to_distance, upload_file, upload_file2, move_files_to_folder
import altair as alt


#Input Block__________________________________________________________________________________________________________________________________________________________________________________________
st.markdown("<h1 style='text-align: center;'>Environmental Census</h1>", unsafe_allow_html=True)
st.header('A Bioinformatics Tool for Synthetic Biology Risk Assessments', divider='rainbow')
st.write(':blue[Developed by John Docter, University of Colorado Boulder]')
st.write(':blue[For Troubleshooting and Inquiries, Please Contact john.docter@colorado.edu]')
st.header('Select Metagenome(s) to Analyze')
intake = st.multiselect('Please choose which metagenomes to analyze', 
                        ['Industrial Wastewater', 'Wastewater Treatment Plant', 'River'])
'You selected: ', str(intake)

choices = [choice.strip().lower() for choice in intake]

#Analysis______________________________________________________________________________________________
home_dir = '/home/anna/Documents/JGI_soil_genomes' 


for mg_to_analyze in choices:
    if mg_to_analyze == 'industrial wastewater':
        metagenome_name = 'reference_diamond_analysis_output' #-> folder
        home_dir = '/home/anna/Documents/JGI_soil_genomes' 
        # IW = '/home/anna/Documents/JGI_soil_genomes/IW_Metagenome'
        IW = '/home/anna/Documents/JGI_soil_genomes/' + mg_to_analyze + '_metagenome_bins'
        # abspath = os.path.abspath()
    #Automatically Makes a folder for the metagenomic bins to upload into, then asks the user for files and moves those files into the bins folder_____________________________________________________###
        os.chdir(home_dir)
        if os.path.exists(IW):
            shutil.rmtree(IW)
            os.mkdir(IW)
        else:
            os.mkdir(IW)
        


        #Below is the original file uplaod using tkinter
        # file_paths = upload_file(home_dir, mg_to_analyze)
          # upload_location = IW
        # move_files_to_folder(file_paths, upload_location)
        st.header(mg_to_analyze.title() + ' Analysis', divider='gray')
        uploaded_file = st.file_uploader("Please upload the Biome .faa file(s) you would like to analyze against", type='.faa', accept_multiple_files=True, key='IW')
        if not uploaded_file:
            st.stop()
        if uploaded_file:
           for f in uploaded_file:
                temp_dir = tempfile.mkdtemp()
                path = os.path.join(temp_dir, f.name)
                with open(path, "wb") as file:
                        file.write(f.getvalue())
                shutil.move(path, IW)
                shutil.rmtree(temp_dir)
        else:
            st.write('Waiting on File Upload')
        
    ##__________________________________________________________________________________#Diamond Analysis
        os.chdir(home_dir)
        if os.path.exists(metagenome_name):
            shutil.rmtree(metagenome_name)
            os.mkdir(metagenome_name)
        else:
            os.mkdir(metagenome_name)
        # home_dir = home_dir + "/" + metagenome_name 

        os.chdir(IW)
        with st.spinner('Diamond Sequence Aligner Matching Biome.faa files to Unitprot Reference'):
            diamond = diamond_impl(IW, '') #-> Takes in the path and directory
        st.success('Biome Sequences Aligned and Present Sequences Identified')

    # #________________________________________________________________________________# Creating reference functional profile

        dmnd_folder = '/home/anna/Documents/JGI_soil_genomes/reference_diamond_analysis_output'
        ff_name = 'functional_profiles'
        functional_folder = '/home/anna/Documents/JGI_soil_genomes/functional_profiles'
        name = mg_to_analyze + '_metagenome'

        for item in os.listdir(IW):
            if item.endswith(('_matches.tsv', '.dmnd')):
                source = os.path.join(IW, item)
                destination = os.path.join(dmnd_folder, item)
                shutil.move(source, destination)
        

        output = genome_extractor_ref(dmnd_folder, name, home_dir, functional_folder, ff_name)
        st.success('Enzymatic Profile Created')

        # os.chdir(home_dir)
        # if os.path.exists(functional_folder):
        #     shutil.rmtree(functional_folder)
        #     os.mkdir(ff_name)
        # else:#makes a new directory called metagenome_name
        #     os.mkdir(ff_name)

        for item in os.listdir(dmnd_folder):
            if item.endswith('_profile'):
                source = os.path.join(dmnd_folder, item)
                destination = os.path.join(functional_folder, item)
                shutil.move(source, destination)
            
    #-----------------------------------------------------------------------------------# Organism to analyze matches and functional profile creation 
        synbio = '/home/anna/Documents/JGI_soil_genomes/' + mg_to_analyze + '_synbio_inputs_and_outputs'
        name = 'Synbio'
        syn_folder_name = mg_to_analyze + '_synbio_inputs_and_outputs'
        desired_location2 = '/home/anna/Documents/JGI_soil_genomes'

        os.chdir(home_dir)
        if os.path.exists(synbio):
            print('Synbio directory already exists')
        else:#makes a new directory called metagenome_name
            os.mkdir(syn_folder_name)

       
 
        uploaded_file = st.file_uploader("Please upload the synbio .faa file you would like to analyze", type='.faa', key = 'IW_syn')
        if not uploaded_file: 
            st.stop()
        if uploaded_file:
            temp_dir = tempfile.mkdtemp()
            path = os.path.join(temp_dir, uploaded_file.name)
            with open(path, "wb") as f:
                f.write(uploaded_file.getvalue())
        

        path2 = os.path.join(synbio, uploaded_file.name)
        if os.path.exists(path2):
            pass
        else:
            shutil.move(path, synbio)
            shutil.rmtree(temp_dir)

        os.chdir(synbio) 
        with st.spinner('Diamond Aligner Matching Synbio.faa Sequences to Unitpro Reference'):
            diamond_syn = diamond_impl(synbio, name) #diamond_syn = synbio
        st.success('Synbio Sequences Aligned')
        output2 = genome_extractor_syn(diamond_syn, name, home_dir)
        st.success('Synbio Enzymatic Profile Created')


        for item in os.listdir(synbio):
            if item.endswith('_profile'):
                source = os.path.join(synbio, item)
                destination = os.path.join(functional_folder, item)
                shutil.move(source, destination)
    #____________________________________________________________________________________#Distance Scoring
        synbio_binary = '/home/anna/Documents/JGI_soil_genomes/functional_profiles/Synbio_functional_profile'
        [distance_list_for_synbio, new_loc ]= pass_to_distance(synbio_binary, name, desired_location2, mg_to_analyze)
        # st.balloons()

        

    elif mg_to_analyze == 'wastewater treatment plant':
        metagenome_name = 'reference_diamond_analysis_output' #-> folder
        home_dir = '/home/anna/Documents/JGI_soil_genomes' 
        # IW = '/home/anna/Documents/JGI_soil_genomes/IW_Metagenome'
        IW = '/home/anna/Documents/JGI_soil_genomes/' + mg_to_analyze + '_metagenome_bins'
        # abspath = os.path.abspath()
    #Automatically Makes a folder for the metagenomic bins to upload into, then asks the user for files and moves those files into the bins folder_____________________________________________________###
        os.chdir(home_dir)
        if os.path.exists(IW):
            shutil.rmtree(IW)
            os.mkdir(IW)
        else:
            os.mkdir(IW)

        st.header(mg_to_analyze.title() + ' Analysis', divider='gray')
        uploaded_file = st.file_uploader("Please upload the Biome .faa file(s) you would like to analyze against", type='.faa', accept_multiple_files=True, key='wwtp')
        if not uploaded_file:
            st.stop()
        if uploaded_file:
           for f in uploaded_file:
                temp_dir = tempfile.mkdtemp()
                path = os.path.join(temp_dir, f.name)
                with open(path, "wb") as file:
                        file.write(f.getvalue())
                shutil.move(path, IW)
                shutil.rmtree(temp_dir)
        else:
            st.write('Waiting on File Upload')

    ##__________________________________________________________________________________#Diamond Analysis
        os.chdir(home_dir)
        if os.path.exists(metagenome_name):
            shutil.rmtree(metagenome_name)
            os.mkdir(metagenome_name)
        else:
            os.mkdir(metagenome_name)
        # home_dir = home_dir + "/" + metagenome_name 

        os.chdir(IW)
        with st.spinner('Diamond Sequence Aligner Matching Biome.faa files to Unitprot Reference'):
            diamond = diamond_impl(IW, '') #-> Takes in the path and directory
        st.success('Biome Sequences Aligned and Present Sequences Identified')
    # #________________________________________________________________________________# Creating reference functional profile

        dmnd_folder = '/home/anna/Documents/JGI_soil_genomes/reference_diamond_analysis_output'
        ff_name = 'functional_profiles'
        functional_folder = '/home/anna/Documents/JGI_soil_genomes/functional_profiles'
        name = mg_to_analyze + '_metagenome'

        for item in os.listdir(IW):
            if item.endswith(('_matches.tsv', '.dmnd')):
                source = os.path.join(IW, item)
                destination = os.path.join(dmnd_folder, item)
                shutil.move(source, destination)
        

        output = genome_extractor_ref(dmnd_folder, name, home_dir, functional_folder, ff_name)
        st.success('Enzymatic Profile Created')
        # os.chdir(home_dir)
        # if os.path.exists(functional_folder):
        #     shutil.rmtree(functional_folder)
        #     os.mkdir(ff_name)
        # else:#makes a new directory called metagenome_name
        #     os.mkdir(ff_name)

        for item in os.listdir(dmnd_folder):
            if item.endswith('_profile'):
                source = os.path.join(dmnd_folder, item)
                destination = os.path.join(functional_folder, item)
                shutil.move(source, destination)
            
    #-----------------------------------------------------------------------------------# Organism to analyze matches and functional profile creation 
        synbio = '/home/anna/Documents/JGI_soil_genomes/' + mg_to_analyze + '_synbio_inputs_and_outputs'
        name = 'Synbio'
        syn_folder_name = mg_to_analyze + '_synbio_inputs_and_outputs'
        desired_location2 = '/home/anna/Documents/JGI_soil_genomes'

        os.chdir(home_dir)
        if os.path.exists(synbio):
            print('Synbio directory already exists')
        else:#makes a new directory called metagenome_name
            os.mkdir(syn_folder_name)

        uploaded_file = st.file_uploader("Please upload the synbio .faa file you would like to analyze", type='.faa', key='wwtp_syn')
        if not uploaded_file: 
            st.stop()
        if uploaded_file:
            temp_dir = tempfile.mkdtemp()
            path = os.path.join(temp_dir, uploaded_file.name)
            with open(path, "wb") as f:
                f.write(uploaded_file.getvalue())
        
        path2 = os.path.join(synbio, uploaded_file.name)
        if os.path.exists(path2):
            pass
        else:
            shutil.move(path, synbio)
            shutil.rmtree(temp_dir)

        os.chdir(synbio) 
        with st.spinner('Diamond Aligner Matching Synbio.faa Sequences to Unitpro Reference'):
            diamond_syn = diamond_impl(synbio, name) #diamond_syn = synbio
        st.success('Synbio Sequences Aligned')
        output2 = genome_extractor_syn(diamond_syn, name, home_dir)
        st.success('Synbio Enzymatic Profile Created')

        for item in os.listdir(synbio):
            if item.endswith('_profile'):
                source = os.path.join(synbio, item)
                destination = os.path.join(functional_folder, item)
                shutil.move(source, destination)
                
        synbio_binary = '/home/anna/Documents/JGI_soil_genomes/functional_profiles/Synbio_functional_profile'
        [distance_list_for_synbio, new_loc ]= pass_to_distance(synbio_binary, name, desired_location2, mg_to_analyze)
        # st.success(mg_to_analyze + ' Synbio Analysis Complete')


    elif mg_to_analyze == 'river':

        metagenome_name = 'reference_diamond_analysis_output' #-> folder
        home_dir = '/home/anna/Documents/JGI_soil_genomes' 
        # IW = '/home/anna/Documents/JGI_soil_genomes/IW_Metagenome'
        IW = '/home/anna/Documents/JGI_soil_genomes/' + mg_to_analyze + '_metagenome_bins'
        # abspath = os.path.abspath()
    #Automatically Makes a folder for the metagenomic bins to upload into, then asks the user for files and moves those files into the bins folder_____________________________________________________###
        os.chdir(home_dir)
        if os.path.exists(IW):
            shutil.rmtree(IW)
            os.mkdir(IW)
        else:
            os.mkdir(IW)

        # # Example usage
        st.header(mg_to_analyze.title() + ' Analysis', divider='gray')
        uploaded_file = st.file_uploader("Please upload the Biome .faa file(s) you would like to analyze against", type='.faa', accept_multiple_files=True, key='river')
        if not uploaded_file:
            st.stop()
        if uploaded_file:
           for f in uploaded_file:
                temp_dir = tempfile.mkdtemp()
                path = os.path.join(temp_dir, f.name)
                with open(path, "wb") as file:
                        file.write(f.getvalue())
                shutil.move(path, IW)
                shutil.rmtree(temp_dir)
        else:
            st.write('Waiting on File Upload')

    ##__________________________________________________________________________________#Diamond Analysis
        os.chdir(home_dir)
        if os.path.exists(metagenome_name):
            shutil.rmtree(metagenome_name)
            os.mkdir(metagenome_name)
        else:
            os.mkdir(metagenome_name)
        # home_dir = home_dir + "/" + metagenome_name 

        os.chdir(IW)
        with st.spinner('Diamond Sequence Aligner Matching Biome.faa files to Unitprot Reference'):
            diamond = diamond_impl(IW, '') #-> Takes in the path and directory
        st.success('Biome Sequences Aligned and Present Sequences Identified')
    # #________________________________________________________________________________# Creating reference functional profile

        dmnd_folder = '/home/anna/Documents/JGI_soil_genomes/reference_diamond_analysis_output'
        ff_name = 'functional_profiles'
        functional_folder = '/home/anna/Documents/JGI_soil_genomes/functional_profiles'
        name = mg_to_analyze + '_metagenome'

        for item in os.listdir(IW):
            if item.endswith(('_matches.tsv', '.dmnd')):
                source = os.path.join(IW, item)
                destination = os.path.join(dmnd_folder, item)
                shutil.move(source, destination)
        

        output = genome_extractor_ref(dmnd_folder, name, home_dir, functional_folder, ff_name)
        st.success('Enzymatic Profile Created')

        # os.chdir(home_dir)
        # if os.path.exists(functional_folder):
        #     shutil.rmtree(functional_folder)
        #     os.mkdir(ff_name)
        # else:#makes a new directory called metagenome_name
        #     os.mkdir(ff_name)

        for item in os.listdir(dmnd_folder):
            if item.endswith('_profile'):
                source = os.path.join(dmnd_folder, item)
                destination = os.path.join(functional_folder, item)
                shutil.move(source, destination)
            
    #-----------------------------------------------------------------------------------# Organism to analyze matches and functional profile creation 
        synbio = '/home/anna/Documents/JGI_soil_genomes/' + mg_to_analyze + '_synbio_inputs_and_outputs'
        name = 'Synbio'
        syn_folder_name = mg_to_analyze + '_synbio_inputs_and_outputs'
        desired_location2 = '/home/anna/Documents/JGI_soil_genomes'

        os.chdir(home_dir)
        if os.path.exists(synbio):
            print('Synbio directory already exists')
        else:#makes a new directory called metagenome_name
            os.mkdir(syn_folder_name)

        
        uploaded_file = st.file_uploader("Please upload the synbio .faa file you would like to analyze", type='.faa', key='river_syn')
        if not uploaded_file: 
            st.stop()
        if uploaded_file:
            temp_dir = tempfile.mkdtemp()
            path = os.path.join(temp_dir, uploaded_file.name)
            with open(path, "wb") as f:
                f.write(uploaded_file.getvalue())
        
        path2 = os.path.join(synbio, uploaded_file.name)
        if os.path.exists(path2):
            pass
        else:
            shutil.move(path, synbio)
            shutil.rmtree(temp_dir)

        os.chdir(synbio) 
        with st.spinner('Diamond Aligner Matching Synbio.faa Sequences to Unitpro Reference'):
            diamond_syn = diamond_impl(synbio, name) #diamond_syn = synbio
        st.success('Synbio Sequences Aligned')
        output2 = genome_extractor_syn(diamond_syn, name, home_dir)
        st.success('Synbio Enzymatic Profile Created')
        for item in os.listdir(synbio):
            if item.endswith('_profile'):
                source = os.path.join(synbio, item)
                destination = os.path.join(functional_folder, item)
                shutil.move(source, destination)
                
        synbio_binary = '/home/anna/Documents/JGI_soil_genomes/functional_profiles/Synbio_functional_profile'
        [distance_list_for_synbio, new_loc ]= pass_to_distance(synbio_binary, name, desired_location2, mg_to_analyze)
        # st.success(mg_to_analyze.title() + ' Synbio Analysis Complete')

    else:
        pass

##Vizualization______________________________________________________________________________________________________________________________________________________###
st.header('Difference Scoring', divider = 'gray')
first_choice = None
second_choice = None
third_choice = None
try:
    first_choice = choices[0]
    # first_choice
    second_choice = choices[1]
    # second_choice
    third_choice = choices[2]
    # third_choice

except IndexError:
    pass
    
diff_score1 = None
diff_score2 = None
diff_score3 = None 
try:
    diff_score1 = pd.read_csv('/home/anna/Documents/JGI_soil_genomes/' + first_choice + '_synbio_inputs_and_outputs/Absolute_Difference_Comparison_Score.txt', sep='\t')
    diff_score1['Biome'] = first_choice
    # diff_score1
    diff_score2 = pd.read_csv('/home/anna/Documents/JGI_soil_genomes/' + second_choice + '_synbio_inputs_and_outputs/Absolute_Difference_Comparison_Score.txt', sep='\t')
    diff_score2['Biome'] = second_choice
    # diff_score2
    diff_score3 = pd.read_csv('/home/anna/Documents/JGI_soil_genomes/' + third_choice + '_synbio_inputs_and_outputs/Absolute_Difference_Comparison_Score.txt', sep='\t')
    diff_score3['Biome'] = third_choice
    # diff_score3
except TypeError:
    pass






try:
    frames = [diff_score1, diff_score2, diff_score3]
    # frames
    combined= pd.concat(frames)
    combined_sorted = combined.sort_values(by='Metagenome Bin ID')
    combined_sorted = combined_sorted.drop('Unnamed: 0', axis=1)
    combined_sorted['Metagenome Bin ID'] = combined_sorted['Metagenome Bin ID'].str.replace('_matches.tsv', '')
    combined_sorted

    # st.bar_chart(combined_sorted, x = 'Metagenome Bin ID', y= 'Difference Score', color='Biome') 
    chart = alt.Chart(combined_sorted).mark_bar().encode(x='Metagenome Bin ID', y='Difference Score', color=alt.Color('Biome',  scale=alt.Scale(scheme='blues')))
    st.altair_chart(chart, use_container_width=True)
    st.toast('Analysis Complete!', icon='🧬')
except ValueError:
    st.write('Waiting on input')

##Deletes all the files for the next run. Turn off if you want files to analyze, turn on if you are going to do multiple runs. 
#Will have to delete files manually if you turn it off, else program gets confused
 
os.chdir(home_dir)
for file in os.listdir(home_dir):
     if file.endswith(('_output', 'outputs', '_bins', '_profiles')):
        shutil.rmtree(os.path.join(home_dir, file))

