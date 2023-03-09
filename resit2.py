import os
import pandas as pd
import numpy as np
import subprocess
import re
#got  my AP from NCBI, should now work
os.environ['NCBI_API_KEY'] = 'cb6b860872dfccb5f04396d6f056d4b36c08'

###Introducing the program name, purpose
while True:
    print("Welcome to OQSIL-Online protein Quality aSsesIng tooL")
    #and taking values from the user
    oqsil = input("Please enter the name of the protein you would like to analyze: ")
    takson = input("Please enter the Taxon ID(txid): ")
    ##Winning some time
    print("Connecting to NCBI")
    print(f"Searching for {oqsil} fastas for txid{takson}")
    # sending request to ncbi for the with protein name and takson
    query = f'{oqsil}[Protein Name] AND txid{takson}[Organism]'
    search_cmd = f'esearch -db protein -query "{query}" | efetch -format fasta'
    ##calculating the number of found fastas
    result_soni = os.popen(f'{search_cmd} | grep ">" | wc -l').read().strip()
    # clarifying whether user they wants to download all fastas
    print("Generating results")
    print(f'Overall:{result_soni} {oqsil} fastas found')
    #Downloading all what user wants if he/she wants
    skachay_vsyo = input(f'Do you want to download all of the {oqsil} fastas for for txid {takson} ? (y/n): ')

    if skachay_vsyo.lower() == 'y':
        print("Downloading all the fastas")
        # downloading all fastas
        skachivaem = f'{search_cmd} > {oqsil}.fa'
        os.system(skachivaem)
        ##Saving the  output
        print(f'All {result_soni} fastas for txid {takson} have been saved to {oqsil}.fa')
        with open(f'{oqsil}.fa', 'r') as oqsillar:
            fasta_tarkibi = oqsillar.read()
	    # assessing the the length and completeness of each sequence
        while True:
            keyingi_etap = input('What would you like to do next? (a) Look for another protein, (b) Continue for similarity analysis, (c) Exit: ')
            if keyingi_etap.lower() == 'a':
                break
            elif keyingi_etap.lower() == 'b':


            # Taking the downloaded file as an input and naming the output file name
                fasta_fayli = f"{oqsil}.fa"
                clust_fayli = fasta_fayli[:-3] + ".msf"

            # Another, hopefully working clustalo command
                clustalo_ishla = f"clustalo -v -i {fasta_fayli} -o {clust_fayli} --outfmt=msf --threads=20 --force"
                ##So people say this code should now work!!!! I really hope so!
                try:
                    subprocess.run(clustalo_ishla, shell=True, check=True)
                    print("Clustal Omega alignment is finished. Generating the results")
                except subprocess.CalledProcessError as e:
                    print(f"Error running Clustal Omega: {e}")
                # Print the alignment result to screen
                with open(clust_fayli, "r") as f:
                    print(f.read())
                # Does the user wants to save the result?
                natijani_saqla = input("Would You like to save the output? (y/n)")
                if natijani_saqla.lower() == "y":
                    # Saving the output file with the input protein {oqsil} name
                    yangi_clust_fayli = f"{oqsil}_aligned.msf"
                    os.rename(clust_fayli, yangi_clust_fayli)
                    # Printing confirmation message
                    print(f"Alignment result saved as {yangi_clust_fayli}.")

            elif keyingi_etap.lower() == 'c':
                break


    else:
        # identifying whether the user wants to filter the data for specific organism
          organism_kirit = input('Do you want to specify the organism? (y/n): ')
	       ##asking the user whether or not they want spefici fastas for specific organism within this taxon
          if organism_kirit.lower() == 'y':
             ##if they say yes, then asking them to enter the name o f that organism
             organism = input('Enter organism name: ')
	         ###Sending query for this specific organims
             query += f' AND {organism}[Organism]'
             search_cmd = f'esearch -db protein -query "{query}" | efetch -format fasta'
             #counting the number of found results
             result_soni = os.popen(f'{search_cmd} | grep ">" | wc -l').read().strip()
             print(f'Number of results found for {organism}: {result_soni}')
             #again giving a word to user()
             download_spec = input(f'Do you want to download all {result_soni} {oqsil} fastas for {organism}? (y/n): ')
             if download_spec.lower() == 'y': # downloading filteres fastas
                 print(f"Downloading {oqsil} fasta files for {organism} ")
                 skachivaem = f'{search_cmd} > {oqsil}.fa'
                 os.system(skachivaem)
                 print(f'All {result_soni} {oqsil} fastas for {organism} have been downloaded to {oqsil}.fa')
                 with open(f'{oqsil}.fa', 'r') as oqsillar:
                     fasta_tarkibi = oqsillar.read()
		             
                     # assessing the the length and completeness of each sequence
                 while True:
                     keyingi_etap = input('What would you like to do next? (a) Look for another protein, (b) Continue for similarity analysis, (c) Exit: ')
                     if keyingi_etap.lower() == 'a':
                         break
                     elif keyingi_etap.lower() == 'b':


                     # Taking the downloaded file as an input and naming the output file name
                         fasta_fayli = f"{oqsil}.fa"
                         clust_fayli = fasta_fayli[:-3] + ".msf"

                     # Another, hopefully working clustalo command
                         clustalo_ishla = f"clustalo -v -i {fasta_fayli} -o {clust_fayli} --outfmt=msf --threads=20 --force"
                         ##So people say this code should now work!!!! I really hope so!
                         try:
                             subprocess.run(clustalo_ishla, shell=True, check=True)
                             print("Clustal Omega alignment is finished. Generating the results")
                         except subprocess.CalledProcessError as e:
                             print(f"Error running Clustal Omega: {e}")
                         # Print the alignment result to screen
                         with open(clust_fayli, "r") as f:
                             print(f.read())
                         # Does the user wants to save the result?
                         natijani_saqla = input("Would You like to save the output? (y/n)")
                         if natijani_saqla.lower() == "y":
                             # Saving the output file with the input protein {oqsil} name
                             yangi_clust_fayli = f"{oqsil}_aligned.msf"
                             os.rename(clust_fayli, yangi_clust_fayli)
                             # Printing confirmation message
                             print(f"Alignment result saved as {yangi_clust_fayli}.")
 
                     elif keyingi_etap.lower() == 'c':
                         break

