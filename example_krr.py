
import numpy as np
import random
import soap


symbol2atomic_no = {'H':1, 'C':6, 'N':7, 'O':8, 'S':16}
filename = 'dsgdb7ae2.xyz'
z, r, all_energy = [], [], []

print('read data')
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
                z.append(list(atomic_no))
                r.append(list(coords))

print('pick training and prediction sets')
train_no = 100
pred_no = 30
ind_sample = random.sample(range(len(all_energy)), train_no+pred_no)
ind_train = random.sample(ind_sample, train_no)
ind_pred = list(set(ind_sample).difference(ind_train))
train_energy = np.asarray([all_energy[i] for i in ind_train])
pred_energy = np.asarray([all_energy[i] for i in ind_pred])

print('compute representation')
repr_train = [soap.reprf(zi, ri) for (zi,ri) in zip([z[i] for i in ind_train], [r[i] for i in ind_train])]
repr_pred  = [soap.reprf(zi, ri) for (zi,ri) in zip([z[i] for i in ind_pred ], [r[i] for i in ind_pred ])]

print('compute kernel matrices')
K = np.asarray(soap.kernelf(repr_train))
L = np.asarray(soap.kernelf(repr_train,repr_pred))

print('make predictions')
lmbda = 10**-6
K = K + lmbda * np.eye(len(K))
alpha = np.linalg.solve(K, train_energy)
f = np.dot(L.transpose(),alpha)

print('stats')
print('Relative MAE:', sum(abs((pred_energy-f)/pred_energy))/len(f)*100, '%')
print('MAE:', sum(abs(pred_energy-f))/len(f))
print('RMSE:', (sum((pred_energy-f)**2)/len(f))**0.5)

