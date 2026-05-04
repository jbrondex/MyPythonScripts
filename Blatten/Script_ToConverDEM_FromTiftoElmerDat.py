import rasterio
import csv

# Chemin vers votre TIFF découpé
input_tif = "/home/brondexj/Documents/Expertises/Blatten/DEMs/pre-event-fromAdrien-smooth-cropped.tif"   ###"chemin/vers/votre_fichier_decoupe.tif"
# Chemin vers le fichier texte de sortie
output_txt = "/home/brondexj/Documents/Expertises/Blatten/DEMs/DEM_surf_2025_Reconstructed_FromAdrien_MJacquemart_Smoothed.dat"   ###"chemin/vers/sortie.txt"

with rasterio.open(input_tif) as src:
    data = src.read(1)  # Lit la première bande (pour un raster mono-bande)
    transform = src.transform
    height, width = data.shape
    nodata = src.nodata  # Valeur NoData (peut être None)

    # Calculate the Y coordinate for the bottom row (to start from the bottom)
    _, y_bottom = transform * (0, height - 1)  # Y for the last row

    with open(output_txt, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        # Iterate rows from bottom to top (ascending Y)
        for row in range(height - 1, -1, -1):
            for col in range(width):
                x, y = transform * (col, row)  # X, Y coordinates
                z = data[row, col]
                # Write only if not NoData
                if nodata is None or z != nodata:
                    writer.writerow([x, y, z])