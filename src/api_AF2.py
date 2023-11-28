import requests



def requesting_structure_(accession_name: str) -> tuple:
    """
    Request structure and confidence data for a given UniProt accession.
    Args:
        accession_name (str): UniProt accession.

    Returns:
        tuple: A tuple containing a list of confidence metrics and the URL to download the structure.
    """
    # Construct URLs for structure and residue data
    url_structure = f"https://alphafold.ebi.ac.uk/api/prediction/{accession_name}?key=AIzaSyCeurAJz7ZGjPQUtEaerUkBZ3TaBkXrY94"
    url_residue = f"https://alphafold.ebi.ac.uk/api/uniprot/summary/{accession_name}.json?key=AIzaSyCeurAJz7ZGjPQUtEaerUkBZ3TaBkXrY94"

    # Send HTTP requests for structure and residue data
    response_structure = requests.get(url_structure)
    response_residue = requests.get(url_residue)
    
    if response_structure.status_code == 200 and response_residue.status_code == 200:
        data_structure = response_structure.json()
        data_residue = response_residue.json()
        
        # Get confidence metrics
        confidence_type = data_residue['structures'][0]['summary'].get('confidence_type')
        confidence_avg_local_score = data_residue['structures'][0]['summary'].get('confidence_avg_local_score')
        confidence_list = [confidence_type, confidence_avg_local_score]

        # Get URL to download structure
        data_url = data_structure[0]['cifUrl']

        return confidence_list, data_url
    else:
        print(f"Request failed with status code (Structure): {response_structure.status_code}")
        print(f"Request failed with status code (Residue): {response_residue.status_code}")
        return [], None
