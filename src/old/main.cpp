#include <iostream>
#include <algorithm>
#include <math.h>

#include "setup.h"
#include "run.h"
#include "solver.h"
#include "stats.h"
#include "stratify.h"
#include "molecule.h"
#include "descriptor.h"

#include "soap_kernel_function.h"
#include <vector>

#define BELOW_5_NON_H 43

int main()
{
    char filename[] = "dsgdb7ae2.xyz";
    int molecules_no = 7102;
    int train_no = 50;
    int validate_no = 10;

	// setup

	Molecule **mol_arr = read_molecules(filename,molecules_no);

    /* stratify and divide into training and validation arrays */
    std::sort(mol_arr, mol_arr+molecules_no, compare_molecules);
    /* molecules below 5 non-H atoms (lowest energy) all go into training set */
    Molecule *train_mol[train_no];
    Molecule *validate_mol[validate_no];
    if (train_no <= BELOW_5_NON_H)
        stratify(mol_arr, train_mol, train_no, validate_mol, validate_no);
    else
    {
        int large_mol_no = molecules_no - BELOW_5_NON_H;
        Molecule *large_mol_arr[large_mol_no];
        for (int mol_idx=0; mol_idx<molecules_no; mol_idx++)
        {
            if (mol_idx<BELOW_5_NON_H)
                train_mol[mol_idx] = mol_arr[mol_idx];
            else large_mol_arr[mol_idx-BELOW_5_NON_H] = mol_arr[mol_idx];
        }

        int large_train_no = train_no - BELOW_5_NON_H;
        int large_mol_sample_no = large_train_no + validate_no;
        Molecule *take_mol[large_mol_sample_no];
        Molecule **sample_large_mol;
        if (0<large_mol_sample_no && large_mol_sample_no<large_mol_no)
        {
            int leave_no = large_mol_no - large_mol_sample_no;
            Molecule *leave_mol[leave_no];
            stratify(large_mol_arr, take_mol, large_mol_sample_no, leave_mol, leave_no);
            sample_large_mol = take_mol;
        } else if (large_mol_sample_no == large_mol_no)
            sample_large_mol = large_mol_arr;
        else exit(1);

        Molecule *large_train_mol[large_train_no];
        stratify(sample_large_mol, large_train_mol, large_train_no, validate_mol, validate_no);
        for (int mol_idx=0; mol_idx<large_train_no; mol_idx++)
            train_mol[mol_idx+BELOW_5_NON_H] = large_train_mol[mol_idx];
    }


    bnu::vector<double> energy(train_no);
    std::vector<intvector> train_at_no(train_no);
    std::vector<dmatrix> train_coords(train_no);
    for (int mol_idx=0; mol_idx<train_no; mol_idx++)
    {
        energy(mol_idx) = train_mol[mol_idx]->energy;
		train_coords[mol_idx].resize(train_mol[mol_idx]->atoms_no);
        for (int i=0; i<train_mol[mol_idx]->atoms_no; i++)
        {
            train_at_no[mol_idx].push_back(index2atomic_no(train_mol[mol_idx]->atom_types[i]));

			train_coords[mol_idx][i].push_back(train_mol[mol_idx]->ff_coords[i].get<0>());
			train_coords[mol_idx][i].push_back(train_mol[mol_idx]->ff_coords[i].get<1>());
			train_coords[mol_idx][i].push_back(train_mol[mol_idx]->ff_coords[i].get<2>());
        }
    }
    dmatrix Kmat = soap_kernel_function(train_at_no, train_coords, train_at_no, train_coords, true);
    bnu::matrix<double> K(train_no,train_no);
    for (int i=0; i<train_no; i++)
        for (int j=0; j<train_no; j++)
            K(i,j) = Kmat[i][j];



    bnu::vector<double> energy_p(validate_no);
    std::vector<intvector> validate_at_no(validate_no);
    std::vector<dmatrix> validate_coords(validate_no);
    for (int mol_idx=0; mol_idx<validate_no; mol_idx++)
    {
        energy_p(mol_idx) = validate_mol[mol_idx]->energy;
		validate_coords[mol_idx].resize(validate_mol[mol_idx]->atoms_no);
        for (int i=0; i<validate_mol[mol_idx]->atoms_no; i++)
        {
            validate_at_no[mol_idx].push_back(index2atomic_no(validate_mol[mol_idx]->atom_types[i]));

			validate_coords[mol_idx][i].push_back(validate_mol[mol_idx]->ff_coords[i].get<0>());
			validate_coords[mol_idx][i].push_back(validate_mol[mol_idx]->ff_coords[i].get<1>());
			validate_coords[mol_idx][i].push_back(validate_mol[mol_idx]->ff_coords[i].get<2>());
        }
    }
    dmatrix Lmat = soap_kernel_function(train_at_no, train_coords, validate_at_no, validate_coords, false);
    bnu::matrix<double> L(train_no,validate_no);
    for (int i=0; i<train_no; i++)
        for (int j=0; j<validate_no; j++)
            L(i,j) = Lmat[i][j];

    free_mol_array(mol_arr, molecules_no);

	// end setup

    Params params;
    params.lamdba = pow(10,-3);
    params.zeta = 1;
    double diag[] = {1,1,1,1,1};    //H,C,N,O,S
    params.diag = diag;

    /* TRAINING */
#if VERBOSE
    std::cout << "solving linear system" << std::endl;
#endif // VERBOSE
    bnu::identity_matrix<double> eye(train_no);
    K += params.lamdba*eye;    //apply ridge parameter
    bnu::vector<double> alpha (train_no);
    int res = solve_linear_system(K,alpha,energy);
    if (res == 0)
    {
        /* VALIDATION */
#if VERBOSE
        std::cout << "producing predictions" << std::endl;
#endif // VERBOSE
        bnu::vector<double> f = prod(trans(L),alpha);

        /* stats */
#if VERBOSE
        output_plot_data(energy_p,f);
#endif // VERBOSE
        produce_stats(energy_p,f);
    }
    else
    {
#if VERBOSE
        std::cout << "cannot solve the linear system" << std::endl;
#endif // VERBOSE
        Stats s;
        double dmax = std::numeric_limits<double>::max();
        s.re_max = dmax;
        s.mre = dmax;
        s.mae = dmax;
        s.rmse = dmax;
    }

    return 0;
}
