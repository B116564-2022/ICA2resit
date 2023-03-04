
import os

#got  my AP from NCBI, should now work 
os.environ['NCBI_API_KEY'] = 'cb6b860872dfccb5f04396d6f056d4b36c08'

###Introducing the program name, purpose and taking values from the users
import os

while True:
    print("Welcome to OQSIL-Online protein Quality aSsesIng tooL")
    oqsil = input("Please enter the name of the protein you would like to analyze: ")
    takson = input("Please enter the Taxon ID(txid): ")
    print(f"Searching for {oqsil} fastas for txid{takson}")

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
        print(f'All {num_results} fastas for txid {takson} have been saved to {oqsil}.fa')
        with open(f'{oqsil}.fa', 'r') as oqsillar:
            fasta_tarkibi = oqsillar.read()
    else:
        # identifying whether the user wants to filter the data for specific organism
        organism_kirit = input('Do you want to specify the organism? (y/n): ')

        if organism_kirit.lower() == 'y':
            organism = input('Enter organism name or taxid: ')
            query += f' AND {organism}[Organism]'
            search_cmd = f'esearch -db protein -query "{query}" | efetch -format fasta'
            num_results = os.popen(f'{search_cmd} | grep ">" | wc -l').read().strip()
            print(f'Number of results found for {organism}: {num_results}')
            download_spec = input(f'Do you want to download all {num_results} {oqsil} fastas for {organism}? (y/n): ')
            if download_spec.lower() == 'y':
                # downloading filteres fastas
                print(f"Downloading {oqsil} fasta files for {organism} ")
                skachivaem = f'{search_cmd} > {oqsil}.fa'
                os.system(skachivaem)
                print(f'All {num_results} {oqsil} fastas for {organism} have been downloaded to {oqsil}.fa')
                with open(f'{oqsil}.fa', 'r') as oqsillar:
                    fasta_tarkibi = oqsillar.read()
        else:
            # does user want to search for another protein?
            yana_qidir = input('Would You like to continue your research on another protein? (y/n): ')
            if yana_qidir.lower() == 'n':
                break

