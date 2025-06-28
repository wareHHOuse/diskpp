import os
import vtkmodules.all as vtk
import subprocess

# Obtenez le chemin absolu du fichier en cours
current_file_path = os.path.abspath(__file__)
# Obtenez le répertoire du fichier en cours
current_directory = os.path.dirname(current_file_path)
# Changez le répertoire de travail actuel pour le répertoire du fichier en cours
os.chdir(current_directory)

# Lire le fichier VTK
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName('simplex_l7.vtk')
reader.Update()

# Récupérer les données du maillage
unstructured_grid = reader.GetOutput()

# Lire le fichier VTK
with open('simplex_l7.vtk', "r")  as vtk_file:
    lines = vtk_file.readlines()
    
# Trouver la ligne qui commence par "POINTS"
points_line = None
for line in lines:
    if line.startswith("POINTS"):
        points_line = line
        break

# Si la ligne "POINTS" est trouvée, extraire les coordonnées des points
if points_line:
    num_points = int(points_line.split()[1])
    points_data = [line.split() for line in lines[lines.index(points_line) + 1:lines.index(points_line) + num_points + 1]]

# Trouver la ligne qui commence par "CELLS"
cells_line = None
for line in lines:
    if line.startswith("CELLS"):
        cells_line = line
        break

# Si la ligne "CELLS" est trouvée, extraire les informations des cellules
if cells_line:
    num_cells, _ = map(int, cells_line.split()[1:])  # Ne pas exclure le premier chiffre (nombre de sommets)
    cells_data = [line.split() for line in lines[lines.index(cells_line) + 1:lines.index(cells_line) + num_cells + 1]]
    
    

# Ouvrir un fichier de sortie pour les coordonnées des points en écriture
with open("simplex_l7.txt", "w") as output_file:
    
    # Écrire le nombre de points, le nombre de cellules, et "4" sur la première ligne
    output_file.write(f"{len(points_data)} {num_cells} 4\n")
    
    for point_data in points_data:
        # Les coordonnées sont dans les trois premières colonnes
        x, y, z = map(float, point_data[:3])
        output_file.write(f"{x} {y}\n")

    for cell_data in cells_data:
        adjusted_cell_data = [cell_data[0]] + [str(int(idx) + 1) for idx in cell_data[1:]]
        # Écrire toutes les données des cellules, y compris le nombre de sommets
        output_file.write(" ".join(adjusted_cell_data) + "\n")

    # Extract the points
    points = unstructured_grid.GetPoints()

    # Calculate the dimensions of the domain
    x_coords = [points.GetPoint(i)[0] for i in range(points.GetNumberOfPoints())]
    y_coords = [points.GetPoint(i)[1] for i in range(points.GetNumberOfPoints())]

    min_x, max_x = min(x_coords), max(x_coords)
    min_y, max_y = min(y_coords), max(y_coords)

    print(f"Min X: {min_x}, Max X: {max_x}")
    print(f"Min Y: {min_y}, Max Y: {max_y}")

    # Define lists to store point indices on each side
    indices_left = []
    indices_right = []
    indices_bottom = []
    indices_top = []

    # Define a threshold to determine which points are on the sides
    threshold = 1e-6

    for i in range(points.GetNumberOfPoints()):
        x, y, _ = points.GetPoint(i)
        if abs(x - min_x) < threshold:
            indices_left.append(i+1)
        if abs(x - max_x) < threshold:
            indices_right.append(i+1)
        if abs(y - min_y) < threshold:
            indices_bottom.append(i+1)
        if abs(y - max_y) < threshold:
            indices_top.append(i+1)

    # # Écrire les numéros des points de la frontière dans un fichier
    # with open("bassin.txt", "a") as frontiere_file:
    output_file.write(" ".join(map(str, indices_left)) + "\n")
    output_file.write(" ".join(map(str, indices_right)) + "\n")
    output_file.write(" ".join(map(str, indices_bottom)) + "\n")
    output_file.write(" ".join(map(str, indices_top)) + "\n")

##################################################################################

print("Les coordonnées des points et les informations des cellules ont été écrites dans le fichier 'output.txt'.")
print("Les numéros des points de la frontière ont été enregistrés dans 'frontiere_points.txt'.")
