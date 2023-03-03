import os

while True:
    print("Welcome to OQSIL-Online protein Quality aSsesIng tooL")
    оqsil = input("Pleasе  enter the name of the protein you would like to analyze: ")
    takson = input("Please enter the Taxon ID(txid): ")
    print(f"Searching {oqsil} fastas for txid{takson}")
    print("Counting the outputs")

    # search for dehydrogenase fastas for given txid
    query = f'{oqsil}[Protein Name] AND txid{takson}[Organism]'
    search_cmd = f'esearch -db protein -query "{query}" | efetch -format fasta'
    num_results = os.popen(f'{search_cmd} | grep ">" | wc -l').read().strip()

    # ask user if they want to download all fastas
    print(f'Number of results found: {num_results}')
    download_all = input('Do you want to download all of the fastas? (y/n): ')

    if download_all.lower() == 'y':
        # download all fastas
        download_cmd = f'{search_cmd} > oqsil.fa'
        os.system(download_cmd)
        print(f'All {num_results} fastas have been downloaded to {oqsil}.fa.')
    else:
        # ask user if they want to specify the organism
        organism_input = input('Do you want to specify the organism? (y/n): ')

        if organism_input.lower() == 'y':
            organism = input('Enter organism name or taxid: ')
            query += f' AND {organism}[Organism]'
            search_cmd = f'esearch -db protein -query "{query}" | efetch -format fasta'
            num_results = os.popen(f'{search_cmd} | grep ">" | wc -l').read().strip()
            print(f'Number of results found for {organism}: {num_results}')
            download_spec = input(f'Do you want to download all {num_results} fastas for {organism}? (y/n): ')
            if download_spec.lower() == 'y':
                # download specified fastas
                download_cmd = f'{search_cmd} > oqsilspec.fa'
                os.system(download_cmd)
                print(f'All {num_results} fastas for {organism} have been downloaded to oqsilspec.fa.')
        else:
            # ask user if they want to search for another protein
            search_again = input('Do you want to search for another protein? (y/n): ')
            if search_again.lower() == 'n':
                break
