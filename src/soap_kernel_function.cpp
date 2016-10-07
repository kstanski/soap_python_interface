#include <iostream>
#include <limits>
#include <math.h>
#include <boost/numeric/ublas/matrix.hpp>

//#include "solver.h"
#include "descriptor.h"
#include "local_similarity.h"
#include "structural_similarity.h"
//#include "stats.h"
//#include "setup.h"
//#include "run.h"

#include "soap_kernel_function.h"
#include "molecule.h"
#include "power_spectrum.h"

#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


Descriptor *repr2desc(mol_repr_arr repr_arr)
{
	Descriptor *desc_arr = create_descriptor_arr(repr_arr.size());
	for (size_t i=0; i<repr_arr.size(); i++)    //molecules
	{
		for (size_t j=0; j<repr_arr[i].size(); j++)    //atoms
		{
			for (size_t k=0; k<repr_arr[i][j].size(); k++)    //atom types
			{
				Power_spectrum ps = (Power_spectrum) malloc(PS_LEN*sizeof(ps_element_type));
				for (size_t l=0; l<PS_LEN; l++)    //power spectrum elements
				{
					ps[l] = repr_arr[i][j][k][l];
				}
				desc_arr[i][j][k] = ps;
			}
		}
	}
	return desc_arr;
}


mol_repr soap_repr_func(intvector atomic_no, dmatrix coords)
{
    Molecule *mol = vectors2molecule(atomic_no, coords);
    mol_repr repr = molecule2repr(mol);
    free(mol);
	return repr;
}

/*
mol_repr_arr soap_repr_func_parallel(intvector_arr atomic_no, dmatrix_arr coords)
{
    int obj_no = atomic_no.size();
    mol_repr_arr repr_arr(obj_no);
#if VERBOSE
    std::cout << "computing power spectra descriptors" << std::endl;
#endif // VERBOSE
    Molecule **mol = vectors2molecules(atomic_no, coords);
    #pragma omp parallel for schedule(dynamic)
    for (int mol_idx=0; mol_idx<obj_no; mol_idx++)
    {
        repr_arr[mol_idx] = molecule2repr(mol[mol_idx]);
    }
    free_mol_array(mol, obj_no);
	return repr_arr;
}
*/

dmatrix soap_kernel_function(mol_repr_arr train_repr_arr, mol_repr_arr validate_repr_arr = mol_repr_arr(0))
{
	bool training = false;
	if (validate_repr_arr.size() == 0) training = true;
	double zeta = 1;
    double diag[] = {1,1,1,1,1};    //H,C,N,O,S

	int train_no = train_repr_arr.size();
	Descriptor *train_desc = repr2desc(train_repr_arr);

	int validate_no = 0;
	Descriptor *validate_desc;
	if (!training)
	{
		validate_desc = repr2desc(validate_repr_arr);
		validate_no = validate_repr_arr.size();
	}

#if VERBOSE
    std::cout << "computing self similarity" << std::endl;
#endif // VERBOSE
    /* training */
    double **LS = create_local_similarity_array(train_desc,train_no,diag);
    double *SS = create_structural_similarity_array(train_desc,LS,train_no,diag);

    /* validation */
    double **LS_p;
    double *SS_p;
    if (!training)
    {
        LS_p = create_local_similarity_array(validate_desc,validate_no,diag);
        SS_p = create_structural_similarity_array(validate_desc,LS_p,validate_no,diag);
    }

#if VERBOSE
    std::cout << "cross structural similarity:" << std::endl;
#endif // VERBOSE

    if (training)
    {
#if VERBOSE
		std::cout << "training matrix..." << std::endl;
#endif // VERBOSE
		/* cross structural similarity */
		dmatrix K(train_no);
        for (int i=0; i<train_no; i++) K[i].resize(train_no);
		#pragma omp parallel for schedule(dynamic)
		for (int i=0; i<train_no; i++)
		{
		    K[i][i] = 1.0;
		    for (int j=i+1; j<train_no; j++)
		    {
		        double k_ij = structural_similarity(train_desc[i],train_desc[j],LS[i],LS[j],diag);
		        k_ij /= sqrt(SS[i]*SS[j]);
		        k_ij = pow(k_ij,zeta);
		        K[i][j] = k_ij;
		        K[j][i] = k_ij;
		    }
		}
#if VERBOSE
    std::cout << "done" << std::endl;
#endif // VERBOSE
        free_ls_arr(LS,train_no);
        free_ss_arr(SS);
        return K;
    } else
    {
#if VERBOSE
		std::cout << "validation matrix..." << std::endl;
#endif // VERBOSE
		dmatrix L(train_no);
        for (int i=0; i<train_no; i++) L[i].resize(validate_no);
		#pragma omp parallel for schedule(dynamic)
		for (int i=0; i<train_no; i++)
		{
		    for (int j=0; j<validate_no; j++)
		    {
		        double l_ij = structural_similarity(train_desc[i],validate_desc[j],LS[i],LS_p[j],diag);
		        l_ij /= sqrt(SS[i]*SS_p[j]);
		        L[i][j] = pow(l_ij,zeta);
		    }
		}
#if VERBOSE
		std::cout << "done" << std::endl;
#endif // VERBOSE
		free_ls_arr(LS,train_no);
		free_ss_arr(SS);
		free_ls_arr(LS_p,validate_no);
		free_ss_arr(SS_p);
	    return L;
    }
}


namespace py = pybind11;

PYBIND11_PLUGIN(soap) {
    py::module m("soap", "pybind11 soap plugin");

	m.def("kernelf", &soap_kernel_function, "A SOAP kernel function",
		py::arg("train_repr"), py::arg("pred_repr") = mol_repr_arr(0));
	m.def("reprf", &soap_repr_func, "A SOAP representation function",
		py::arg("atomic_no"), py::arg("coords"));

    return m.ptr();
}
