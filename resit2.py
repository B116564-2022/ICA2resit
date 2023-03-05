import os
import pandas as pd
import numpy as np

#got  my AP from NCBI, should now work
os.environ['NCBI_API_KEY'] = 'cb6b860872dfccb5f04396d6f056d4b36c08'

###Introducing the program name, purpose and taking values from the user

while True:
    print("Welcome to OQSIL-Online protein Quality aSsesIng tooL")
    oqsil = input("Please enter the name of the protein you would like to analyze: ")
    takson = input("Please enter the Taxon ID(txid): ")
    print("Connecting to NCBI")
    print(f"Searching for {oqsil} fastas for txid{takson}")
    print("One more second")
    # sending request to ncbi for the with protein name and takson
    query = f'{oqsil}[Protein Name] AND txid{takson}[Organism]'
    search_cmd = f'esearch -db protein -query "{query}" | efetch -format fasta'
    num_results = os.popen(f'{search_cmd} | grep ">" | wc -l').read().strip()

    # clarifying whether user they wants to download all fastas
    print("Generating results")
    print(f'Number of results found: {num_results}')
    skachay_vsyo = input(f'Do you want to download all of the {oqsil} fastas for for txid {takson} ? (y/n): ')

    if skachay_vsyo.lower() == 'y':
        # downloading all fastas
        print("Downloading all the fastas")
        skachivaem = f'{search_cmd} > {oqsil}.fa'
        os.system(skachivaem)
        print(f'All {num_results} fastas for txid {takson} have been saved to {oqsil}.fa.')
        with open(f'{oqsil}.fa', 'r') as oqsillar:
            fasta_tarkibi = oqsillar.read()
	    # assessing the the length and completeness of each sequence
        check_quality = input(f"Do you want to conduct quality check on the downloaded {oqsil} fasta file? (y/n): ")
        if check_quality.lower() == 'y':
            ###Calculaig the mean lengths of the fastas
            fasta_list = fasta_tarkibi.split('>')[1:]
            fasta_lens = [len(fasta.split('\n', 1)[1].replace('\n', '')) for fasta in fasta_list]
            mean_len = np.mean(fasta_lens)
            print(f"Mean length of the {oqsil} fasta sequences: {mean_len} bp.")
	      ##assessing the completeness of the protein
            completeness = len(fasta_lens) / int(num_results)
            print(f"Completeness of the {oqsil} fasta sequences: {completeness*100:.2f}%.")
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
             num_results = os.popen(f'{search_cmd} | grep ">" | wc -l').read().strip()
             print(f'Number of results found for {organism}: {num_results}')
             download_spec = input(f'Do you want to download all {num_results} {oqsil} fastas for {organism}? (y/n): ')
             if download_spec.lower() == 'y': # downloading filteres fastas
                 print(f"Downloading {oqsil} fasta files for {organism} ")
                 skachivaem = f'{search_cmd} > {oqsil}.fa'
                 os.system(skachivaem)
                 print(f'All {num_results} {oqsil} fastas for {organism} have been downloaded to {oqsil}.fa.')
                 with open(f'{oqsil}.fa', 'r') as oqsillar:
                     fasta_tarkibi = oqsillar.read()
		##assessing the completeness of the protein
                     check_quality = input(f"Do you want to conduct quality check on the downloaded {oqsil} fasta file? (y/n): ")
                     if check_quality.lower() == 'y':
                         # checking the mean length and completeness
                         fasta_list = fasta_tarkibi.split('>')[1:]
                         fasta_lens = [len(fasta.split('\n', 1)[1].replace('\n', '')) for fasta in fasta_list]
                         mean_len = np.mean(fasta_lens)
                         print(f"Mean length of the {oqsil} fasta sequences: {mean_len} bp.")
                         completeness = len(fasta_lens) / int(num_results)
                         print(f"Completeness of the {oqsil} fasta sequences: {completeness*100:.2f}%.")
                         next_step = input('What would you like to do next? (a) Look for another protein, (b) Exit: ')
                         if next_step.lower() == 'a':
                             continue
                         elif next_step.lower() == 'b':
                             break

                    #clarifying what exactly does user want for the next step
                     else:
                         next_step = input('What would you like to do next? (a) Look for another protein, (b) Exit: ')
                         if next_step.lower() == 'a':
                             continue
                         elif next_step.lower() == 'b':
                             break
