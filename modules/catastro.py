import os
import requests
import urllib.parse

BASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def get_building_gml(ref_cat, save_folder):
    """
    Get the GML file of a building from the parcel reference
    and save it in the specified folder.

    Parameters
    ----------
    ref_cat : str
        Reference of the parcel.
    save_folder : str
        Folder where the GML file will be saved.

    Returns
    -------
    None
    """

    # Define the query to get the building parts
    base_url = "http://ovc.catastro.meh.es/INSPIRE/wfsBU.aspx"
    params = {
        "service": "wfs",
        "version": "2",
        "request": "getfeature",
        "STOREDQUERIE_ID": "GETBUILDINGPARTBYPARCEL",
        "refcat": ref_cat,
        "srsname": "EPSG::25829"
    }
    query = f"{base_url}?{urllib.parse.urlencode(params)}"

    try:
        response = requests.get(query, timeout=10)  # Define a timeout of 10 seconds

        if response.status_code == 200:
            print(f"The request was successful for the parcel {ref_cat}")

            # Check if the folder exists
            if not os.path.exists(save_folder):
                raise FileNotFoundError(f"The folder {save_folder} does not exist")

            file_name = os.path.join(save_folder, f"Building_{ref_cat}.gml")

            # Check if the file already exists and ask the user to overwrite it or not
            if os.path.exists(file_name):
                overwrite = input(f"The file {file_name} already exists. Do you want to overwrite it? (y/n): ")
                if overwrite.lower() != "y":
                    print("The file was not overwritten")
                    return

            # Save the file
            with open(file_name, "wb") as file:
                file.write(response.content)
                print(f"The file was saved in {file_name}")
        
        else:
            print(f"Request error: {response.status_code}\n{response.text}")

    except requests.RequestException as e:
        print(f"Connection error: {e}")

    
if __name__ == "__main__":
    # Define the reference of the parcel
    ref_cat = "3530001TG3433S"
    
    # Get the building GML
    out_folder = os.path.join(BASE_DIR, "data", "gml", "buildings")
    get_building_gml(ref_cat, out_folder)