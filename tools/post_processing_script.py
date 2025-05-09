import os
import re
import csv

# Please insert directory with your `.out` files, use pwd to discover what it is
DATA_DIR = ['/home/eng/rubinam4/article_simulation_MgOCo/ghost_sim/MgCO/nk_111_sc_1.0000',
            '/home/eng/rubinam4/article_simulation_MgOCo/ghost_sim/MgCO/nk_111_sc_1.0000/results_verified']

# Output CSV path - if you just name it, it will write the file to the current folder you're in
CSV_OUTPUT = 'processed_results_with_lap.csv'

# Regular expressions stuff - to identify relevant results from output files
scf_pattern = re.compile(r'converged SCF energy\s*=\s*(-?\d+\.\d+)')
mp2_pattern = re.compile(r'(E\(KMP2\)|E\(DFKMP2\))\s*=\s*(-?\d+\.\d+)\s+.*E_corr\s*=\s*(-?\d+\.\d+)')
mp2_laplace_pattern = re.compile(r'E\(Laplace-SOS-KMP2\)\s*=\s*(-?\d+\.\d+)')
memory_allocated_pattern = re.compile(r'max_memory\s*(-?\d+)\s*MB')

# Prepare to collect rows in a list
rows = []

# Header row - we can even pass on 'filename' column
header = ['filename', 'structure', 'basis_set', 'SCF_energy', 'MP2_energy', 'E_corr', 'MP2_Laplace', 'Memory_allocated']
rows.append(header)

# Process each .out file
for data_directory in DATA_DIR:
    for filename in os.listdir(data_directory):
        if filename.startswith('MgOCO_') and filename.endswith('.out'):
            filepath = os.path.join(data_directory, filename)
            with open(filepath, 'r') as file:
                scf_energy = ''
                mp2_energy = ''
                mp2_corr = ''
                mp2_laplace = ''
                memory_allocated = ''

                for line in file:
                    if 'converged SCF energy' in line:
                        match = scf_pattern.search(line)
                        if match:
                            scf_energy = match.group(1)
                    elif 'E(KMP2)' in line or 'E(DFKMP2)' in line:
                        match = mp2_pattern.search(line)
                        if match:
                            mp2_energy = match.group(2)
                            mp2_corr = match.group(3)

                    elif 'E(Laplace-SOS-KMP2)' in line:
                        match = mp2_laplace_pattern.search(line)
                        if match:
                            mp2_laplace = match.group(1)
                            # print(mp2_laplace)
                            break

                    elif 'max_memory' in line:
                        match = memory_allocated_pattern.search(line)
                        if match:
                            memory_allocated = match.group(1)
                            # print(mp2_laplace)

                # Extract metadata from filename - which is the structure and the basis set used
                structure_match = re.search(r'(2l1a_\d+(?:_(?:CO|MgO)_ghost)?)', filename)
                basis_match = re.search(r'cc-pv[a-z]{2}', filename)

                structure = structure_match.group(1) if structure_match else ''
                basis_set = basis_match.group(0) if basis_match else ''

                # if we gave up on filename column, put it in comment or something
                rows.append([
                    filename,
                    structure,
                    basis_set,
                    scf_energy,
                    mp2_energy,
                    mp2_corr,
                    mp2_laplace,
                    memory_allocated
                ])

# Write rows to CSV file - open it with windows, not linux; otherwise it's ugly and unreadable
with open(CSV_OUTPUT, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(rows)

# A log that makes sure that everything's alright
print(f"CSV saved to {CSV_OUTPUT}")

"""
import os
import re
import pandas as pd

# Directory with your `.out` files
DATA_DIR = '/home/eng/rubinam4/article_simulation_MgOCo/ghost_sim/MgCO/nk_111_sc_1.0000'

# Regular expressions
scf_pattern = re.compile(r'converged SCF energy\s*=\s*(-?\d+\.\d+)')
mp2_pattern = re.compile(r'(E\(KMP2\)|E\(DFKMP2\))\s*=\s*(-?\d+\.\d+)\s+.*E_corr\s*=\s*(-?\d+\.\d+)')

# To hold extracted data
results = []

# Process each .out file
for filename in os.listdir(DATA_DIR):
    if filename.startswith('MgOCO_verify') and filename.endswith('.out'):
        filepath = os.path.join(DATA_DIR, filename)
        with open(filepath, 'r') as file:
            scf_energy = None
            mp2_energy = None
            mp2_corr = None

            for line in file:
                if 'converged SCF energy' in line:
                    match = scf_pattern.search(line)
                    if match:
                        scf_energy = float(match.group(1))
                elif 'E(KMP2)' in line or 'E(DFKMP2)' in line:
                    match = mp2_pattern.search(line)
                    if match:
                        mp2_energy = float(match.group(2))
                        mp2_corr = float(match.group(3))

            # Extract metadata from filename
            structure = re.search(r'(2l1a_\d+(?:_(?:CO|MgO)_ghost)?)', filename)
            basis_set = re.search(r'cc-pv\w+', filename)

            results.append({
                'filename': filename,
                'structure': structure.group(1) if structure else None,
                'basis_set': basis_set.group(0) if basis_set else None,
                'SCF_energy': scf_energy,
                'MP2_energy': mp2_energy,
                'E_corr': mp2_corr
            })

            # Extract metadata from filename
            structure = re.search(r'(2l1a\w*)', filename)
            basis_set = re.search(r'cc-pv\w+', filename)

            results.append({
                'filename': filename,
                'structure': structure.group(1) if structure else None,
                'basis_set': basis_set.group(0) if basis_set else None,
                'SCF_energy': scf_energy,
                'MP2_energy': mp2_energy,
                'E_corr': mp2_corr
            })

# Convert to DataFrame
df = pd.DataFrame(results)

# Output as CSV or just print
df.to_csv('processed_results.csv', index=False)
print(df)

"""
