import os
import pandas as pd
import subprocess
import matplotlib.pyplot as plt

#got  my AP from NCBI, should now work
os.environ['NCBI_API_KEY'] = 'cb6b860872dfccb5f04396d6f056d4b36c08'

###Introducing the program name, purpose
while True:
    print("Welcome to OQSIL-Online protein Quality aSsesIng tooL")
    #and taking values from the user
    oqsil = input("Please enter the name of the protein you would like to analyze: ")
    takson = input("Please enter the Taxon ID(txid): ")
    ##Winning some time
    print(f">Connecting to NCBI\n\n>Searching for {oqsil} fastas for txid{takson}")
    # sending request to ncbi for the with protein name and takson
    query = f'{oqsil}[Protein Name] AND txid{takson}[Organism]'
    search_cmd = f'esearch -db protein -query "{query}" | efetch -format fasta'
    ##calculating the number of found fastas
    result_soni = os.popen(f'{search_cmd} | grep ">" | wc -l').read().strip()
    # clarifying whether user they wants to download all fastas
    if result_soni == "0":
        print("Error: No results found for the specified organism or txid. \nEither there is no fastas presented for this specific protein and/or txid. \nPlease check https://www.ncbi.nlm.nih.gov/, specifying protein {oqsil} and txid{txid}\n \nIf these data are presented in NCBI, then Please check your input and try again. \nLets start from the beginning!!!")
        continue
    else:
        print(f">Generating results \nOverall:{result_soni} {oqsil} fastas found. \nPLEASE KEEP IN MIND. If the number of proteins is too high, further generated results, might be not that much beautiful and informative as You expect")
    
    #Downloading all what user wants if he/she wants
    skachay_vsyo = input(f'Do you want to download all of the {oqsil} fastas for for txid {takson} ? (y/n): ')

    if skachay_vsyo.lower() == 'y':
        print("\n>Downloading all the fastas")
        # downloading all fastas
        skachivaem = f'{search_cmd} > {oqsil}.fa'
        os.system(skachivaem)
        ##Saving the  output
        print(f'\n\n>All {result_soni} fastas for txid {takson} have been saved to {oqsil}.fa')
        
        oqsillar = f"{oqsil}.fa"                              
        with open(oqsillar, "r") as f:
            lines = f.readlines()  
        # Creating a dictionary to hold the protein names and their sequences
        proteins = {}
        for line in lines:
            if line.startswith(">"):
                protein_name = line.strip()[1:]
                proteins[protein_name] = ""
            else:
                proteins[protein_name] += line.strip()
        
	    #MIN/MAX/MEAN
                
                
        # Creating a pandas DataFrame from the dictionary
        oqsilcha = pd.DataFrame.from_dict(proteins, orient="index", columns=["Sequence"])
        # Counting the number of nucleotides for each protein
        oqsilcha["Nucleotide Count"] = oqsilcha["Sequence"].apply(len)
        smallest_one = oqsilcha["Nucleotide Count"].min()
        longest_one = oqsilcha["Nucleotide Count"].max()
        oqsil_mean_length = oqsilcha["Nucleotide Count"].mean()
        print(f"The min. length of these fastas is {smallest_one} bp and the longest sequence comprises of {longest_one} bp. \n Mean length of {oqsil} sequences is {oqsil_mean_length} bp")
                           
        with open(f'{oqsil}.fa', 'r') as oqsillar:
            fasta_tarkibi = oqsillar.read()     

        while True:
            keyingi_etap = input('\n\n>>>What would you like to do next? (a) Look for another protein, (b) Continue for similarity analysis, (c) Exit: ')
            if keyingi_etap.lower() == 'a':
                continue
            elif keyingi_etap.lower() == 'b':
            # Again specifying the input and output files for clustalo
            ###So this code will be taking the downloaded file and creating two files
            ##one multiple sequence aligned and one distance matrix. 
                fasta_fayli = f"{oqsil}.fa"
                clust_fayli = fasta_fayli[:-3] + ".msf"
                dist_fayli = f"{oqsil}.dist"
                # Generating Clustalo matrix with the similarity scores for sequence
                clustalo_path = "/usr/bin/clustalo"
                os.environ["PATH"] += os.pathsep + clustalo_path
                clustalo_ishla = f"clustalo -v -i {fasta_fayli} -o {clust_fayli} --outfmt=msf --threads=20 --force --full --distmat-out={dist_fayli}"
                try:
                    subprocess.run(clustalo_ishla, shell=True, check=True, capture_output=True, text=True)
                    print("\n\n>>>Clustal Omega alignment is finished. Generating the results. Please wait")
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
                    yangi_clust_fayli = f"{oqsil}.msf"
                    os.rename(clust_fayli, yangi_clust_fayli)
                    print(f"Alignment result saved as {yangi_clust_fayli}.")
                    ###and the dist matrix
                    yangi_dist_fayli = f"{oqsil}.matrix"
                    os.rename(dist_fayli, yangi_dist_fayli)
                    print(f"Similarity scores saved as {yangi_dist_fayli}.")
                    # Prompt the user to plot the level of conservation between protein sequences
                    plot_con = input("Do you want to plot the level of conservation between protein sequences? (y/n): ")

                    if plot_con.lower() == "y":
                        # Run the command to plot conservation level
                        conservation_plotting = f"plotcon -sformat msf {yangi_clust_fayli} -winsize 16 -graph pdf"
                        
                        os.system(conservation_plotting)
                        # view the generated PDF file on screen
                        view_pdf = "gs plotcon.pdf"
                        os.system(view_pdf)
                        # Display the plot
                        plt.show()
                        print("Conservation plot saved as plotcon.pdf")
                       
                        ### PATMATMOTIFS 
                       
                        
                        # asking the  user if they want to compare protein motifs
                        compare_motifs = input("Do you want to identify protein motifs in your set? (y/n): ")
                        if compare_motifs.lower() == "y":
                            # Setting the path to the patmatmotifs command
                            
                            patmatmotifs_path = "/usr/bin/patmatmotifs"
                            #The output will look like this
                            oqsil_patmat=f"{oqsil}.patmatmotifs"
                            #Tho patmat command itself
                            patmatmotifs_cmd = f"patmatmotifs -sequence {yangi_clust_fayli} -outfile {oqsil_patmat} -full"
                            print(patmatmotifs_cmd)  # Check the command before running it
                            #Running the command, fingers crossed
                            try:
                                subprocess.run(patmatmotifs_cmd, shell=True, check=True)
                                print("patmatmotifs finished successfully.")
                                print(f"The results saved into {oqsil}.patmatmotifs")
                            except subprocess.CalledProcessError as e:
                                print(f"Error running patmatmotifs: {e}")
                    
                    
                    
                    
                    # Ask the user if they want to analyze another protein/the same protein for another txid or exit
                    analyze_another = input("Do you want to analyze another protein/the same protein for another txid or exit? (a/n/e): ")

                    if analyze_another.lower() == "a":
                        continue
                    elif analyze_another.lower() == "n":
                        break
                    elif analyze_another.lower() == "e":
                        print("Thank you for using OQSIL.")
                        break
                    else:
                        print("Invalid input. Please enter 'a' to analyze another protein/the same protein for another txid, 'n' to exit, or 'e' to exit.")


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
             if result_soni == "0":
                 print("Error: No results found for the specified organism or txid. \nEither there are no fastas presented for this specific protein and/or txid. \nPlease check https://www.ncbi.nlm.nih.gov/ specifying protein {oqsil}, species {organism}\nIf these data are presented in NCBI, then Please check your input and try again. \nLets start again!")
                 continue
             else:
                 print(f'Number of results found for {organism}: {result_soni}\n PLEASE KEEP IN MIND. If the number of proteins is too high, further generated results, might be not that much beautiful and informative as You expect')
             #again giving a word to user()
             download_spec = input(f'Do you want to download all {result_soni} {oqsil} fastas for {organism}? (y/n): ')
             if download_spec.lower() == 'y': # downloading filteres fastas
                 print(f"Downloading {oqsil} fasta files for {organism} ")
                 skachivaem = f'{search_cmd} > {oqsil}.fa'
                 os.system(skachivaem)
                 print(f'All {result_soni} {oqsil} fastas for {organism} have been downloaded to {oqsil}.fa')
                 
                 
                 oqsillar = f"{oqsil}.fa"                              
                 with open(oqsillar, "r") as f:
                     lines = f.readlines()  
                 # Creating a dictionary to hold the protein names and their sequences
                 proteins = {}
                 for line in lines:
                     if line.startswith(">"):
                         protein_name = line.strip()[1:]
                         proteins[protein_name] = ""
                     else:
                         proteins[protein_name] += line.strip()

                 #MIN/MAX/MEAN
                         
                         
                 # Creating a pandas DataFrame from the dictionary
                 oqsilcha = pd.DataFrame.from_dict(proteins, orient="index", columns=["Sequence"])
                 # Counting the number of nucleotides for each protein
                 oqsilcha["Nucleotide Count"] = oqsilcha["Sequence"].apply(len)
                 smallest_one = oqsilcha["Nucleotide Count"].min()
                 longest_one = oqsilcha["Nucleotide Count"].max()
                 oqsil_mean_length = oqsilcha["Nucleotide Count"].mean()
                 print(f"The min. length of these fastas is {smallest_one} bp and the longest sequence comprises of {longest_one} bp. \n Mean length of {oqsil} sequences is {oqsil_mean_length} bp")
                                    
                 with open(f'{oqsil}.fa', 'r') as oqsillar:
                     fasta_tarkibi = oqsillar.read()

                     # assessing the the length and completeness of each sequence
                 while True:
                     keyingi_etap = input('What would you like to do next? (a) Look for another protein, (b) Continue for similarity analysis, (c) Exit: ')
                     if keyingi_etap.lower() == 'a':
                         continue
                     elif keyingi_etap.lower() == 'b':
                         fasta_fayli = f"{oqsil}.fa"
                         clust_fayli = fasta_fayli[:-3] + ".msf"
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
                             yangi_clust_fayli = f"{oqsil}.msf"
                             os.rename(clust_fayli, yangi_clust_fayli)
                             print(f"\nAlignment result saved as {yangi_clust_fayli}.")
                             ###and the dist matrix
                             yangi_dist_fayli = f"{oqsil}.matrix"
                             os.rename(dist_fayli, yangi_dist_fayli)
                             print(f"Similarity scores saved as {yangi_dist_fayli}.")
                             # Prompt the user to plot the level of conservation between protein sequences
                             
                             plot_con = input("\nDo you want to plot the level of conservation between protein sequences? (y/n): ")
                        
                             if plot_con.lower() == "y":
                                # Visualising the plot of conservation level on screen
                                 conservation_plotting = f"plotcon -sformat msf {yangi_clust_fayli} -winsize 16 -graph pdf" 
                                 os.system(conservation_plotting)
                                 # view the generated PDF file on screen
                                 view_pdf = "gs plotcon.pdf"
                                 os.system(view_pdf)
                                 # Display the plot
                                 plt.show()
                                 print("Conservation plot saved as plotcon.pdf")
                               
                                 
                              
                                
                              #PATMATMOTIFS
                                
                               
                                
                               
                                # asking the  user if they want to compare protein motifs
                                 compare_motifs = input("Do you want to identify protein motifs in your set? (y/n): ")
                                 if compare_motifs.lower() == "y":
                                     # Set the path to the patmatmotifs command
                                     patmatmotifs_path = "/usr/bin/patmatmotifs"
                                     oqsil_patmat=f"{oqsil}.patmatmotifs"
                                     patmatmotifs_cmd = f"patmatmotifs -sequence {yangi_clust_fayli} -outfile {oqsil_patmat}"
                                     print(patmatmotifs_cmd)  # Check the command before running it
                                     try:
                                         subprocess.run(patmatmotifs_cmd, shell=True, check=True)
                                         print("patmatmotifs finished successfully.")
                                         print(f"The results saved into {oqsil}.patmatmotifs")
                                     except subprocess.CalledProcessError as e:
                                         print(f"Error running patmatmotifs: {e}")

                           
                            
                           # Asking the user, if they want to analyze another protein/the same protein for another txid or exit
                             analyze_another = input("Do you want to analyze another protein/the same protein for another txid or exit? (a/n/e): ")
                        
                             if analyze_another.lower() == "a":
                                 continue
                             elif analyze_another.lower() == "n":
                                 break
                             elif analyze_another.lower() == "e":
                                 print("Thank you for using OQSIL.")
                                 break
                             else:
                                 print("Invalid input. Please enter 'a' to analyze another protein/the same protein for another txid, 'n' to exit, or 'e' to exit.")
                     elif keyingi_etap.lower() == 'c':
                         break
