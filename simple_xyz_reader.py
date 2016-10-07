
symbol2atomic_no = {'H':1, 'C':6, 'N':7, 'O':8, 'S':16}

filename = 'dsgdb7ae2.xyz'
all_atomic_no, all_coords, all_energy = [], [], []

with open(filename) as f:
    for line in f:
        words = line.split()
        l = len(words)
        if l == 1:
            atoms_no = int(words[0])
            atoms_count = 0
            atomic_no = []
            coords = []
        elif l == 2:
            all_energy.append(float(words[1]))
        elif l == 7:
            atoms_count += 1
            atomic_no.append(symbol2atomic_no.get(words[0]))
            coords.append([float(words[idx]) for idx in (1,2,3)])
            if atoms_count == atoms_no:
                all_atomic_no.append(list(atomic_no))
                all_coords.append(list(coords))

import soap
train_no = 5
training_atomic_no = all_atomic_no[:train_no]
print(training_atomic_no[0])
training_coords = all_coords[:train_no]
print(training_coords[0])

M = soap.kernelf(training_atomic_no,training_coords,training_atomic_no,training_coords,True)
print(M)

