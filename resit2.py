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


            # Again specifying the input and output files for clustalo
            ###So this code will be taking the downloaded file and creating two files
            ##one multiple sequence aligned and one distance matrix. Really dont know
            ##how else i can do it to work like show first names of two prot pairs AND
            ##then how similar of different are they
            ##but could generate only the distance matrix, which shows the similarity with each
            #fasta but not separately, Still, understandable
                fasta_fayli = f"{oqsil}.fa"
                clust_fayli = fasta_fayli[:-3] + ".clustal"
                dist_fayli = f"{oqsil}.dist"

                # Generating dist matrix with the similarity scores for sequence
                clustalo_ishla = f"clustalo -v -i {fasta_fayli} -o {clust_fayli} --outfmt=msf --threads=20 --force --full --distmat-out={dist_fayli}"

                try:
                    subprocess.run(clustalo_ishla, shell=True, check=True, capture_output=True, text=True)
                    print("Clustal Omega alignment is finished. Generating the results. Please wait")
                    # Taking out  similarity scores from the generated results
                    similarity_scores = []
                    with open(clust_fayli, "r") as f:
                        for line in f:
                            if line.startswith('#'):
                                score = float(line.split(':')[-1].strip())
                                similarity_scores.append(f"{score:.2%}")
                    print(f"Similarity scores between fastas: {similarity_scores}")
                    # Distance matrix going to screen and might be saved as file at a later stage
                    with open(dist_fayli, "r") as f:
                        print(f.read())
                except subprocess.CalledProcessError as e:
                    print(f"Unfortunately, there was an error running Clustal Omega: {e}")

                # Visualising the aligned sequences to the screen
                with open(clust_fayli, "r") as f:
                    print(f.read())

                #Clarifying whether the user wants the result to be saved?
                natijani_saqla = input("Would you like to save the output? (y/n)")
                #if the user wants to save then downloading results into two different files
                #the aligned pairs and the distance matrix
                if natijani_saqla.lower() == "y":
                    # At this point we are saving the output files with the protein {oqsil} name
                    ##clustal one
                    yangi_clust_fayli = f"{oqsil}_aligned.msf"
                    os.rename(clust_fayli, yangi_clust_fayli)
                    print(f"Alignment result saved as {yangi_clust_fayli}.")
                    ###and the dist matrix
                    yangi_dist_fayli = f"{oqsil}_dist.matrix"
                    os.rename(dist_fayli, yangi_dist_fayli)
                    print(f"Distance matrix saved as {yangi_dist_fayli}.")

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
                        fasta_fayli = f"{oqsil}.fa"
                        clust_fayli = fasta_fayli[:-3] + ".clustal"
                        dist_fayli = f"{oqsil}.dist"

                        # Generating dist matrix with the similarity scores for sequence
                        clustalo_ishla = f"clustalo -v -i {fasta_fayli} -o {clust_fayli} --outfmt=msf --threads=20 --force --full --distmat-out={dist_fayli}"

                        try:
                            subprocess.run(clustalo_ishla, shell=True, check=True, capture_output=True, text=True)
                            print("Clustal Omega alignment is finished. Generating the results. Please wait")
                            # Taking out  similarity scores from the generated results
                            similarity_scores = []
                            with open(clust_fayli, "r") as f:
                                for line in f:
                                    if line.startswith('#'):
                                        score = float(line.split(':')[-1].strip())
                                        similarity_scores.append(f"{score:.2%}")
                            print(f"Similarity scores between fastas: {similarity_scores}")
                            # Distance matrix going to screen and might be saved as file at a later stage
                            with open(dist_fayli, "r") as f:
                                print(f.read())
                        except subprocess.CalledProcessError as e:
                            print(f"Unfortunately, there was an error running Clustal Omega: {e}")

                        # Visualising the aligned sequences to the screen
                        with open(clust_fayli, "r") as f:
                            print(f.read())

                        #Clarifying whether the user wants the result to be saved?
                        natijani_saqla = input("Would you like to save the output? (y/n)")
                        #if the user wants to save then downloading results into two different files
                        #the aligned pairs and the distance matrix
                        if natijani_saqla.lower() == "y":
                            # At this point we are saving the output files with the protein {oqsil} name
                            ##clustal one
                            yangi_clust_fayli = f"{oqsil}_aligned.msf"
                            os.rename(clust_fayli, yangi_clust_fayli)
                            print(f"Alignment result saved as {yangi_clust_fayli}.")
                            ###and the dist matrix
                            yangi_dist_fayli = f"{oqsil}_dist.matrix"
                            os.rename(dist_fayli, yangi_dist_fayli)
                            print(f"Distance matrix saved as {yangi_dist_fayli}.")
                     elif keyingi_etap.lower() == 'c':
                         break
