from paraview.simple import *
import os
import re

# --- USER SETTINGS ---
input_folder     = "./../data"
output_folder    = "./snap"
field_to_display = "C"
image_resolution = [1920, 1080]
color_map_range  = [-1.0, 1.0]
crop_fraction    = 0.25  # show only 25% of domain height around the interface

# --- Prepare output directory ---
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# --- Load only .vts files with physXXXX.vts pattern ---
files = sorted([
    f for f in os.listdir(input_folder)
    if f.lower().endswith(".vts") and re.match(r"phys\d{4}\.vts", f)
])

for f in files:
    file_path = os.path.join(input_folder, f)
    reader = OpenDataFile(file_path)

    # --- Display & color setup ---
    display = Show(reader)
    display.Representation = 'Surface'
    ColorBy(display, ('POINTS', field_to_display))

    lut = GetColorTransferFunction(field_to_display)
    lut.ApplyPreset("Inferno (matplotlib)", True)
    lut.RescaleTransferFunction(color_map_range)
    HideScalarBarIfNotNeeded(lut)

    # --- Render view setup ---
    view = GetActiveViewOrCreate('RenderView')
    view.Background = [0, 0, 0]

    # compute bounds & center
    b   = reader.GetDataInformation().GetBounds()
    cx  = 0.5 * (b[0] + b[1])
    cy  = 0.5 * (b[2] + b[3])
    cz  = 0.5 * (b[4] + b[5])

    # look down –Y so X→right, Z→up
    dy = b[3] - b[2]
    view.CameraPosition           = [cx, cy + dy, cz]
    view.CameraFocalPoint         = [cx, cy, cz]
    view.CameraViewUp             = [0, 0, 1]
    view.CameraParallelProjection = 1
    view.Update()

    # --- tight zoom around interface ---
    half_x = 0.5 * (b[1] - b[0]) * crop_fraction
    half_z = 0.5 * (b[5] - b[4]) * crop_fraction
    view.CameraParallelScale = max(half_x, half_z)

    # --- Save screenshot with original numeric name ---
    match = re.search(r"phys(\d{4})\.vts", f)
    if match:
        png_name = f"phys{match.group(1)}.png"
    else:
        png_name = f.replace(".vts", ".png")

    SaveScreenshot(
        os.path.join(output_folder, png_name),
        view,
        ImageResolution=image_resolution,
        TransparentBackground=True
    )

    # Clean up reader and display
    Delete(reader)
    del reader
    Hide(display)

print(f"Saved {len(files)} zoomed Inferno‐colormap screenshots to {output_folder}")

