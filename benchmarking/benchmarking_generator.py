import random
from app.app import process_request
from tqdm import trange
import os

def random_time(h_start=9, h_end=20):
    """
    Generate a random time (HH:MM) between two hours (inclusive).
    Returns a string in the format 'HH:MM'.
    """
    hour = random.randint(h_start, h_end)
    minute = random.randint(0, 59)
    return f"{hour:02d}:{minute:02d}"

def random_location():
    LOCATIONS = [
        "4, Calle Sófocles, Málaga",
        "10, Calle Cantimpla, Málaga",
        "5, Calle Juan Muñoz Herrera, Málaga",
        "30, Calle Arganda, Málaga",
        "1, Calle ALbahaca, Málaga",
        "20, Calle Barroso, Málaga",
        "16, Calle Especería, Málaga",
        "8, Calle Santiago, Málaga",
        "19, Calle Jinetes, Málaga",
        "11, Calle San vicente de Paul, Málaga",
        "47, Calle Rocío, Málaga",	
        "10, Calle Mayoral, Málaga"
    ]

    return random.sample(LOCATIONS, 2)


if __name__ == "__main__":

    iterations = 250
    list_times = []

    # Route to the file where you want to save the results
    result_file = os.path.join(os.path.dirname(__file__), "benchmark_data.csv")

    # If the file already exists, delete it
    if os.path.exists(result_file):
        os.remove(result_file)

    # Write the header to the file
    with open(result_file, "w") as f:
        header = (
            "t00_start,"
            "t01_load_coords,"
            "t02_get_sun_vector,"
            "t03_get_existing_sun_vector,"
            "t04_load_weighted_graph,"
            "t05_calculate_routes,"
            "t06_remove_repeated_routes,"
            "t07_convert_routes_to_coordinates,"
            "t08_calculate_routes_metrics,"
            "t09_load_shadows,"
            "t10_render_map\n"
        )
        f.write(header)

    for i in trange(iterations, desc="Benchmarking"):
        time = random_time()
        origin, destination = random_location()
        date = "2025-06-01"
        info = process_request(origin,
                               destination,
                               date,
                               time,
                               True)
        
        # Process the time subtracting the first time in the info
        base_time = info["t00_start"]
        for key in info.keys():
                info[key] = info[key] - base_time

        # Write the values to the file
        line = ",".join(str(v) for v in info.values())
        with open(result_file, "a") as f:
            f.write(line + "\n")