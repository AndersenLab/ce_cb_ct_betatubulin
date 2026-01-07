#!/usr/bin/env python3

import sys
import os

# Output directories
SCRIPT_DIR = "scripts/structure"
IMAGE_DIR = "data/structure/raw_images"


def prereq():
    return """load data/structure/AF3_Predictions/AF3_ABZ/briggsae.cif, briggsae
load data/structure/AF3_Predictions/AF3_ABZ/elegans.cif, elegans
load data/structure/AF3_Predictions/AF3_ABZ/tropicalis.cif, tropicalis
load data/structure/AF3_Predictions/AF3_ABZ/contrortus.cif, contrortus

align briggsae,   elegans
align tropicalis, elegans
align contrortus, elegans

set_view (\
    -0.396658033,    0.840651333,   -0.368715227,\
     0.907275438,    0.297913790,   -0.296806395,\
    -0.139665246,   -0.452259541,   -0.880877316,\
    -0.000083029,    0.000087211, -196.839874268,\
    -6.681622505,    8.408699036,    1.224905014,\
  -10000.000000000, 10000.000000000,  -20.000000000 )

hide all"""


def draw(name, variants, suffix=""):
    mutants = []
    for variant in variants:
        # ignore data for other species
        if variant[0] != name:
            continue

        position = variant[1]
        color = variant[2]
        mutants.append(f"show spheres, {name} and resi {position}\n")
        mutants.append(f"color {color}, {name} and resi {position}\n")

    filename = f"{IMAGE_DIR}/{name}{suffix}" if suffix else f"{IMAGE_DIR}/{name}"

    return f"""
show cartoon, {name}
color 0xFFFFFF, {name} and chain 'A'
show sticks, {name} and chain 'L'
color 0xFFFF00, {name} and chain 'L' and elem C
color 0xFFFF00, {name} and chain 'L' and elem N
color 0xFFFF00, {name} and chain 'L' and elem O
color 0xFFFF00, {name} and chain 'L' and elem S
{"".join(mutants)[:-1]}
png {filename}, width=2000, height=1500, dpi=600, ray=1
hide all"""


def parse(filename):
    variants = []
    with open(filename, "r") as f:
        header = f.readline().strip().split(",")
        col_map = {name: idx for idx, name in enumerate(header)}

        for line in f:
            values = line.strip().split(",")
            species = values[col_map["species"]]
            position = values[col_map["residue"]]

            if not position:
                continue

            color = values[col_map["color_hex_gradient"]].replace("#", "0x")
            label = values[col_map["label"]]
            variants.append((species, position, color, label))
    return variants


def parse_missense(filename):
    variants = []
    with open(filename, "r") as f:
        header = f.readline().strip().split(",")
        col_map = {name: idx for idx, name in enumerate(header)}

        for line in f:
            values = line.strip().split(",")

            # Filter to only Missense variants
            variant_type = values[col_map["ben-1_clean_call"]]
            if variant_type != "Missense":
                continue

            species = values[col_map["species"]]
            position = values[col_map["residue"]]

            if not position:
                continue

            color = values[col_map["color_hex_gradient"]].replace("#", "0x")
            label = values[col_map["label"]]
            variants.append((species, position, color, label))
    return variants


if __name__ == "__main__":
    # Ensure output directories exist
    os.makedirs(IMAGE_DIR, exist_ok=True)

    variants = parse(sys.argv[1])
    print(variants, file=sys.stderr)

    # Build the .pml script content
    pml_content = []
    pml_content.append(prereq())
    pml_content.append(draw("elegans", variants))
    pml_content.append(draw("briggsae", variants))
    pml_content.append(draw("tropicalis", variants))
    pml_content.append(draw("contrortus", variants))

    # Generate missense-focused figures
    missense_variants = parse_missense(sys.argv[1])
    print("\n# Missense variants only", file=sys.stderr)
    print(missense_variants, file=sys.stderr)

    pml_content.append("\n# Missense-focused figures")
    pml_content.append(draw("elegans", missense_variants, "_missense"))
    pml_content.append(draw("briggsae", missense_variants, "_missense"))
    pml_content.append(draw("tropicalis", missense_variants, "_missense"))
    pml_content.append(draw("contrortus", missense_variants, "_missense"))

    # Write the .pml script to scripts/structure/
    pml_output_path = os.path.join(SCRIPT_DIR, "draw_structures.pml")
    with open(pml_output_path, "w") as f:
        f.write("\n".join(pml_content))

    print(f"\nPyMOL script written to: {pml_output_path}", file=sys.stderr)
    print(f"Images will be saved to: {IMAGE_DIR}/", file=sys.stderr)
