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
                #if the user wants to do some cleaning we will run the analysis through this code, taking the downloaded fastas and cleaning it from non aa elements
                yuklangan_fasta = f"{oqsil}.fa"
                # clean fasta sequeneces will be stored here
                tozalangan_fasta = f"cleaned_{oqsil}.fa"
                # through this dictionary.
                toza_sekv = {}

                # reading downloaded fasta file
                with open(yuklangan_fasta, "r") as f:
                    lines = f.readlines()
                    seq_id = None
                    seq = ""
                    for line in lines:
                        line = line.strip()
                        if line.startswith(">"):
                            if seq_id:
                                # removin non-aa characters from the sequence
                                toza_seq = ''.join(filter(lambda x: x in 'ACDEFGHIKLMNPQRSTVWY-', seq.upper()))
                                toza_sekv[seq_id] = toza_seq
                            seq_id = line[1:]
                            seq = ""
                        else:
                            seq += line
                    # processing the final sequence in the file
                    toza_seq = ''.join(filter(lambda x: x in 'ACDEFGHIKLMNPQRSTVWY-', seq.upper()))
                    toza_sekv[seq_id] = toza_seq

                # writing the cleaned sequences to file
                with open(tozalangan_fasta, "w") as f:
                    for seq_id, seq in toza_sekv.items():
                        f.write(f">{seq_id}\n{seq}\n")
                print(f"All fastas within {oqsil}.fa cleaned and stored in cleaned_{oqsil}.fa")
                # running Clustal Omega(at least trying but i dont know why its not working, i suspect its becoz of the clening fail, but how to correct it??????)
                ###i was specifying here the path to clustalo, ("usr/bin/clustalo") but its now gone.  Where? will sort it out later, by now its somehow working |
                clustalo_ishla = f'clustalo -i {tozalangan_fasta} --outfmt=phylip --full > yakuniy.matrica'
                ##running clustal o
                subprocess.run(clustalo_ishla, shell=True, check=True)
                valid_lines = []
                #trying to make this stuff work with all the fastas just skipping over the non standard elements, but its finding a problem in eahc line. Errrrrrr
                with open("yakuniy.matrica", "r") as f:
                    matriks_chiziq = f.readlines()[1:]  
                for line in matriks_chiziq:
                    if "--------------------------------------------------" in line:
                        continue
                    if any(c.isalpha() for c in line):
                        print(f"Skipping line in distance matrix due to invalid format: {line.strip()}")
                        continue
                    valid_lines.append(line)

                # transforming all this list into an array
                n = len(valid_lines)
                matrica = np.zeros((n, n))
                for i, line in enumerate(matriks_chiziq):
                    if "--------------------------------------------------" in line:
                        continue
                    if any(c.isalpha() for c in line):
                        continue
                    matrica[i] = np.array(line.strip().split(), dtype=float)
                print(matrica)###Pleeeeeeeease work.

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
		             ##assessing the completeness of the protein(Do i really need this?)
                     check_quality = input(f"Do you want to conduct quality check on the downloaded {oqsil} fasta file? (y/n): ")
                     if check_quality.lower() == 'y':
                         # checking the mean length and completeness
                         fasta_list = fasta_tarkibi.split('>')[1:]
                         fasta_lens = [len(fasta.split('\n', 1)[1].replace('\n', '')) for fasta in fasta_list]
                         mean_len = np.mean(fasta_lens)
                         print(f"Mean length of the {oqsil} fasta sequences: {mean_len} bp.")
                         completeness = len(fasta_lens) / int(result_soni)
                         print(f"Completeness of the {oqsil} fasta sequences: {completeness*100:.2f}%.")
                         keyingi_etap = input('What would you like to do next? (a) Look for another protein, (b) Continue for similarity analysis, (c) Exit: ')
                         if keyingi_etap.lower() == 'a':
                             continue
                         elif keyingi_etap.lower() == 'b':
                             # path to input FASTA file
                             yuklangan_fasta = f'{oqsil}.fa'
                             # path to output distance matrica file
                             yakuniy_matriks_path = 'oqsil_sim.matrica'
                             # run Clustal Omega
                             clustalo_ishla = f'clustalo -i {yuklangan_fasta} --outfmt=phylip --full > {yakuniy_matriks_path}'
                             subprocess.run(clustalo_ishla, shell=True, check=True)
                             # lets now load the distance matrica into a NumPy array
                             with open(yakuniy_matriks_path, 'r') as f:
                                 matriks_chiziq = f.readlines()[1:]  # skip first line
                                 print(matriks_chiziq)  # add this line to check the contents of matriks_chiziq
                             n = len(matriks_chiziq)
                             matrica = np.zeros((n, n))
                             for i, line in enumerate(matriks_chiziq):
                                 matrica[i] = np.array(line.strip().split(), dtype=float)
                                 # print the distance matrica
                                 print(matrica)
                         elif keyingi_etap.lower() == 'c':
                             break

                    #clarifying what exactly does user want for the next step
                     else:
                         keyingi_etap = input('What would you like to do next? (a) Look for another protein, (b) Continue for similarity analysis, (c) Exit: ')
                         if keyingi_etap.lower() == 'a':
                             continue
                         elif keyingi_etap.lower() == 'b':
                            ##loading the input fasta file
                            yuklangan_fasta = f'{oqsil}.fa'
                            #final cleaned fasta_file
                            yakuniy_matriks_path = f'{oqsil}.matrica'
                            # running Clustal Omega (it doesnt like me)
                            clustalo_ishla = f'clustalo -i {yuklangan_fasta} --outfmt=phylip --full > {yakuniy_matriks_path}'
                            subprocess.run(clustalo_ishla, shell=True, check=True)
                            # lets now try to load the distance matrix into a NumPy array
                            with open(yakuniy_matriks_path, 'r') as f:
                                matriks_chiziq = f.readlines()[1:]  # skippin the first line
                                print(matriks_chiziq)  # adding this line to check the contents of matrix_lines(chiziq is line in my lang)
                            n = len(matriks_chiziq)
                            matrica = np.zeros((n, n))
                            for i, line in enumerate(matriks_chiziq):
                                matrica[i] = np.array(line.strip().split(), dtype=float)
                                # print the distance matrica
                                print(matrica)
                         elif keyingi_etap.lower() == 'c':
                             break
                             ###please woooooork
