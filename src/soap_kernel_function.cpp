#include "soap_kernel_function.h"

#include "descriptor.h"
#include "local_similarity.h"
#include "molecule.h"
#include "power_spectrum.h"
#include "structural_similarity.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <cmath>
#include <iostream>
#include <limits>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <thread>
#include <vector>

typedef struct kernel_data {
   Descriptor *train_desc;
   Descriptor *validate_desc;
   double **LS;
   double **LS_p;
   double *diag;
   double *SS;
   double *SS_p;
   dmatrix *M;
   double zeta;
   int nthreads;
   int id;
   bool training;
} Kernel_data;

void *compute_kernel_matrix(void *kernel_data_ptr)
{
   Kernel_data *kdp = (Kernel_data *) kernel_data_ptr;
   dmatrix *M_ptr = kdp->M;
   int row_no = (* M_ptr).size();
   int col_no = (* M_ptr)[0].size();
   int entries_no;    //no of unique entries in the matrix
   if (kdp->training) entries_no = (row_no*(row_no-1))/2;
   else entries_no = row_no * col_no;

   int row, col;
   for (int i=kdp->id; i<entries_no; i += kdp->nthreads)
   {
      if (kdp->training)
      {
         row = floor((1+sqrt(1+8*i))/2);
         col = round(i-(row*(row-1))/2);
      } else {
         row = i/col_no;
         col = i%col_no;
      }

	  double val = structural_similarity(kdp->train_desc[row],kdp->validate_desc[col],
	  	                                 kdp->LS[row],kdp->LS_p[col],kdp->diag);
	  val /= sqrt(kdp->SS[row]*kdp->SS_p[col]);
      val = pow(val,kdp->zeta);
      (* M_ptr)[row][col] = val;
      if (kdp->training) (* M_ptr)[col][row] = val;
   }
   pthread_exit(NULL);
}

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

	//threads
	int rc;
	int nthreads = std::thread::hardware_concurrency();
	pthread_t threads[nthreads];
	pthread_attr_t attr;
	void *status;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

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

    if (training)
    {
		/* cross structural similarity */
		dmatrix K(train_no);
        for (int i=0; i<train_no; i++)
        {
        	K[i].resize(train_no);
        	K[i][i] = 1.0;    //diagonal
        }

		Kernel_data kd_arr[nthreads];
		for(int i=0; i<nthreads; i++)
		{
			kd_arr[i].train_desc = train_desc;
			kd_arr[i].validate_desc = train_desc;
			kd_arr[i].LS = LS;
			kd_arr[i].LS_p = LS;
			kd_arr[i].diag = diag;
			kd_arr[i].SS = SS;
			kd_arr[i].SS_p = SS;
			kd_arr[i].M = &K;
			kd_arr[i].zeta = zeta;
			kd_arr[i].nthreads = nthreads;
			kd_arr[i].id = i;
			kd_arr[i].training = training;

			rc = pthread_create(&threads[i], &attr, compute_kernel_matrix, (void *)&kd_arr[i]);
			if (rc){
				std::cerr << "Error:unable to create thread," << rc << std::endl;
				exit(-1);
			}
		}

		//wait for the other threads
		pthread_attr_destroy(&attr);
		for(int i=0; i<nthreads; i++)
		{
			rc = pthread_join(threads[i], &status);
			if (rc){
				std::cerr << "Error:unable to join," << rc << std::endl;
				exit(-1);
			}
		}

/*
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
*/
        free_ls_arr(LS,train_no);
        free_ss_arr(SS);
        return K;
    } else
    {

    	/* cross structural similarity */
		dmatrix L(train_no);
        for (int i=0; i<train_no; i++) L[i].resize(validate_no);

		Kernel_data kd_arr[nthreads];
		for(int i=0; i<nthreads; i++)
		{
			kd_arr[i].train_desc = train_desc;
			kd_arr[i].validate_desc = validate_desc;
			kd_arr[i].LS = LS;
			kd_arr[i].LS_p = LS_p;
			kd_arr[i].diag = diag;
			kd_arr[i].SS = SS;
			kd_arr[i].SS_p = SS_p;
			kd_arr[i].M = &L;
			kd_arr[i].zeta = zeta;
			kd_arr[i].nthreads = nthreads;
			kd_arr[i].id = i;
			kd_arr[i].training = training;

			rc = pthread_create(&threads[i], &attr, compute_kernel_matrix, (void *)&kd_arr[i]);
			if (rc){
				std::cerr << "Error:unable to create thread," << rc << std::endl;
				exit(-1);
			}
		}

		// free attribute and wait for the other threads
		pthread_attr_destroy(&attr);
		for(int i=0; i<nthreads; i++)
		{
			rc = pthread_join(threads[i], &status);
			if (rc){
				std::cerr << "Error:unable to join," << rc << std::endl;
				exit(-1);
			}
		}

/*
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
*/
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
