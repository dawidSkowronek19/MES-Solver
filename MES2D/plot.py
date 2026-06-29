import matplotlib.pyplot as plt

# Struktury do przechowywania zrekonstruowanej siatki
x_coords = []
y_coords = []
elements = []

# Słownik do "sklejania" węzłów o tych samych współrzędnych
unique_nodes = {}
current_node_id = 0

def get_node_id(x, y):
    global current_node_id
    # Zaokrąglenie dla bezpieczeństwa floatów
    coord_key = (round(x, 5), round(y, 5))
    
    if coord_key not in unique_nodes:
        unique_nodes[coord_key] = current_node_id
        x_coords.append(x)
        y_coords.append(y)
        current_node_id += 1
        
    return unique_nodes[coord_key]

# Odczyt danych wygenerowanych przez C++
file_path = "./outdir/trianglesList.dat"

try:
    with open(file_path, "r") as file:
        for line in file:
            # Czyszczenie linii z formy (x,y) na same spacje
            clean_line = line.replace('(', ' ').replace(')', ' ').replace(',', ' ')
            coords = [float(val) for val in clean_line.split()]
            
            if len(coords) == 6:
                id1 = get_node_id(coords[0], coords[1])
                id2 = get_node_id(coords[2], coords[3])
                id3 = get_node_id(coords[4], coords[5])
                
                elements.append([id1, id2, id3])

    print("--- Wczytywanie zakończone ---")
    print(f"Liczba węzłów: {len(x_coords)}")
    print(f"Liczba trójkątów: {len(elements)}")

    # Rysowanie siatki strukturalnej
    plt.figure(figsize=(12, 6))
    
    # Twój poprawiony układ przekątnych narysuje się tutaj idealnie
    plt.triplot(x_coords, y_coords, elements, color='blue', marker='o', markersize=4, linewidth=1)
    
    plt.title("Strukturalna Siatka MES 2D")
    plt.xlabel("Oś X")
    plt.ylabel("Oś Y")
    plt.axis('equal') # Wymuszenie proporcji 1:1, by nie zniekształcać komórek
    plt.grid(True, linestyle='--', alpha=0.6)
    
    plt.show()

except FileNotFoundError:
    print(f"Błąd: Nie znaleziono pliku '{file_path}'. Odpal najpierw kompilację i program C++.")