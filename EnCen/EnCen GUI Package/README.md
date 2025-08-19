<img src="87ed0cbc-842e-45e6-8ffd-a7d8254f337d.png" alt="Eco" width="400"/>

Please reach out to john.docter@colorado.edu with any questions 

Introduction
============
GUI_Package is the software for running analyses on your local machine. 

Data heavy and custom analyses can be done using the full code (EnCen HPC). 

Instructions
============
1. Download or git pull GUI_Package from the repository. 
2. Ensure all python libraries and packages in EcoGenoRisk.py (Main script) and Functions.py (Custom functions) are installed
3. Navigate to working directory, or where GUI_Package is downloaded
4. In the terminal, run 'streamlit run EcoGenoRisk.py'
5. A webpage should open and the program should be running. 
6. Follow all prompts to output a threat assessment. 

Useful Tips
===========
"EC_List_For-Parsing.txt" needs to stay where it is relative to the main script. 

"Post Processing Scripts" are for going outside the bounds of the GUI. There are scripts for cleaning data downloaded from the Joint Genome Institute as well as analyzing the data that comes out of the program and visualizing it. The program will run fine without these, they are there in case of need. 
