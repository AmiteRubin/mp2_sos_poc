def parse_vasp_file(file_path):
    lattice_vectors = []
    cell_structure = []

    try:
        with open(file_path, 'r') as file:
            lines = [line.strip() for line in file if line.strip()]  # Remove empty lines

            # Check minimum number of lines in the file
            if len(lines) < 8:
                raise ValueError("File format is invalid: Too few lines")

            # Extract lattice vectors (lines 2-4)
            try:
                lattice_vectors = [
                    list(map(float, lines[2].split())),
                    list(map(float, lines[3].split())),
                    list(map(float, lines[4].split())),
                ]
            except (ValueError, IndexError):
                raise ValueError("Error parsing lattice vectors. Ensure lines 2-4 contain valid numbers.")

            # Extract atom types and counts (lines 5-6)
            try:
                atom_types = lines[5].split()
                atom_counts = list(map(int, lines[6].split()))
                if len(atom_types) != len(atom_counts):
                    raise ValueError("Mismatch between atom types and counts.")
            except (ValueError, IndexError):
                raise ValueError("Error parsing atom types or counts.")

            # Check coordinate type (line 7)
            coordinate_type = lines[7].strip()
            if coordinate_type.lower() != "cartesian":
                raise ValueError("Only Cartesian coordinates are supported.")

            # Extract atomic positions (starting from line 8)
            expected_atoms = sum(atom_counts)
            atom_positions = lines[8:]
            if len(atom_positions) < expected_atoms:
                raise ValueError(f"Expected {expected_atoms} atom positions, but found {len(atom_positions)}.")
            if len(atom_positions) > expected_atoms:
                raise Warning(f"Extra atom positions found. Only the first {expected_atoms} will be used.")

            ghosts_dict = {atm: number for number, atm in enumerate(atom_types)}
            ghost_count = 0

            # Build the structure
            current_atom_index = 0
            for atom_type, count in zip(atom_types, atom_counts):
                orig_atom_type = atom_type
                for _ in range(count):
                    try:
                        is_ghost = False
                        is_augmented = False
                        line_coordinates = atom_positions[current_atom_index].split()
                        position = tuple(map(float, line_coordinates[0:3]))
                        # Check for special chars if line consists of 4 elements after splitting
                        if len(line_coordinates) == 4:
                            special_chars = line_coordinates[3]
                            if '*' in special_chars:
                                is_augmented = True
                            if 'X' in special_chars:
                                is_ghost = True
                        elif len(position) != 3:
                            raise ValueError(
                                f"Invalid number of coordinates for atom at line {8 + current_atom_index}.")
                        # If it's ghost, or augmented - consider it in unique atom type
                        if is_ghost:
                            atom_type = 'X' + str(ghosts_dict[atom_type])
                            ghost_count += 1
                            # is_ghost = False
                        if is_augmented:
                            atom_type = atom_type + '*'  # atom_type.removesuffix('*') + '*'
                            # is_augmented = False
                        cell_structure.append((atom_type, position))  # x,y,z coordinates of position
                        current_atom_index += 1

                        # Reset atom type to unfiltered type for next iteration
                        atom_type = orig_atom_type
                    except (ValueError, IndexError):
                        raise ValueError(f"Error parsing atom position at line {8 + current_atom_index}.")

    except FileNotFoundError:
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    except Exception as e:
        raise RuntimeError(f"An unexpected error occurred while processing the file: {str(e)}")

    ghosts_dict = {atm: 'X' + str(number) for number, atm in enumerate(atom_types)}

    if ghost_count == 0:
        ghosts_dict = []

    return lattice_vectors, cell_structure, ghosts_dict
