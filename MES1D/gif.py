import os
import re
import glob
import matplotlib.pyplot as plt
from PIL import Image

def extract_time_from_filename(filename):
    """
    Wyciąga wartość czasu z nazwy pliku, np. z 'outdir/u0.100000.dat' wyciągnie 0.1
    """
    match = re.search(r"u([-+]?\d*\.\d+|\d+)", os.path.basename(filename))
    return float(match.group(1)) if match else -1.0

def create_fem_gif(data_dir="outdir_timeDependent", output_gif="fem_simulation.gif", fps=10):
    # 1. Znajdź wszystkie pliki .dat pasujące do wzorca u*.dat
    search_pattern = os.path.join(data_dir, "u*.dat")
    files = glob.glob(search_pattern)
    
    if not files:
        print(f"[Błąd]: Nie znaleziono plików u*.dat w katalogu '{data_dir}'!")
        return

    # 2. Posortuj pliki chronologicznie według czasu zapisanego w nazwie
    files.sort(key=extract_time_from_filename)
    print(f"Znaleziono {len(files)} plików do animacji.")

    # 3. Pierwsze przejście: znajdź globalne min/max dla osi Y, żeby wykres nie "skakał"
    y_min, y_max = float('inf'), float('-inf')
    x_min, x_max = float('inf'), float('-inf')
    
    for file_path in files:
        with open(file_path, 'r') as f:
            for line in f:
                if line.strip():
                    x_val, y_val = map(float, line.split())
                    if y_val < y_min: y_min = y_val
                    if y_val > y_max: y_max = y_val
                    if x_val < x_min: x_min = x_val
                    if x_val > x_max: x_max = x_val

    # Dodaj mały margines do osi, żeby wykres wyglądał estetycznie
    y_margin = (y_max - y_min) * 0.1 if y_max != y_min else 1.0
    y_lims = (y_min - y_margin, y_max + y_margin)
    x_lims = (x_min, x_max)

    frames = []
    tmp_frame_dir = "tmp_frames"
    os.makedirs(tmp_frame_dir, exist_ok=True)

    print("Generowanie klatek animacji...")
    # 4. Generowanie wykresu dla każdej klatki
    plt.style.use('seaborn-v0_8-whitegrid') # ładny wygląd wykresu
    fig, ax = plt.subplots(figsize=(8, 5))

    for idx, file_path in enumerate(files):
        t_val = extract_time_from_filename(file_path)
        
        # Wczytywanie danych z pliku
        x_data = []
        u_data = []
        with open(file_path, 'r') as f:
            for line in f:
                if line.strip():
                    x_val, u_val = map(float, line.split())
                    x_data.append(x_val)
                    u_data.append(u_val)

        # Rysowanie klatki
        ax.clear()
        ax.plot(x_data, u_data, label=f"u(x, t)", color='royalblue', linewidth=2)
        ax.set_xlim(x_lims)
        ax.set_ylim(y_lims)
        ax.set_title(f"Rozwiązanie niestacjonarne MES (Time Dependent)", fontsize=12, fontweight='bold')
        ax.set_xlabel("Obszar przestrzenny (x)", fontsize=10)
        ax.set_ylabel("Wartość rozwiązania (u)", fontsize=10)
        
        # Formatowanie podpisu czasu do 4 miejsc po przecinku w legendzie
        ax.legend(title=f"Czas t = {t_val:.4f}", loc="upper right")

        # Zapisz tymczasowy obrazek klatki
        frame_path = os.path.join(tmp_frame_dir, f"frame_{idx:04d}.png")
        plt.savefig(frame_path, dpi=150)
        
        # Otwórz obrazek przez PIL i dodaj do listy klatek gifa
        frames.append(Image.open(frame_path))

    plt.close()

    # 5. Budowanie końcowego pliku GIF
    print(f"Zapisywanie animacji do pliku: {output_gif}...")
    duration_ms = int(1000 / fps) # czas trwania jednej klatki w milisekundach
    
    frames[0].save(
        output_gif,
        save_all=True,
        append_images=frames[1:],
        duration=duration_ms,
        loop=0 # 0 oznacza pętlę nieskończoną
    )

    # 6. Czyszczenie plików tymczasowych PNG
    print("Czyszczenie plików tymczasowych...")
    for frame_img in frames:
        frame_img.close() # Zwolnij uchwyty do plików przed usunięciem
        
    for tmp_file in glob.glob(os.path.join(tmp_frame_dir, "*.png")):
        os.remove(tmp_file)
    os.rmdir(tmp_frame_dir)

    print("Sukces! Animacja została utworzona.")

if __name__ == "__main__":
    # Uruchomienie skryptu dla domyślnego folderu "outdir"
    # Możesz zmienić 'fps', aby przyspieszyć lub zwolnić animację
    create_fem_gif(data_dir="outdir_timeDependent", output_gif="fem_simulation.gif", fps=12)