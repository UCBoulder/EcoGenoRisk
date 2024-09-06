<span style="color:green;">EnCen</span>
=====
 Encen_HPC.py and EnCen_Local.ipynb are the scripts for generating large scale and small scale functional binary matrices and difference scores. 
 
 Analysis_scripts folder, specifically _1_EnCen_Topmatch_profile_creation.py and _2_Dataframe_Creation.py contains the processing scripts for EnCen outputs. 

 faa_extractor is a small script for cleaning JGI downloads and getting only the .faa files. 

1) EnCen Streamlit
   - Proof of concept web application with UI for local analyses
   - Open "EnCen_Streamlit.py" and in the terminal, after navigating to the dir. with the script, type "steamlit run Encen_Streamlit.py"
2) Analysis_Scripts
   - For running after Encen_HPC.py or EnCen_Local.ipynb, depending on your computation resources
   - Takes in the outputs of the above scripts and analyzes the data
   - The graphs folder contains the scripts for turning the output of _2_Dataframe_Creation.py into visualizations. 

 
