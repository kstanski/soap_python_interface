#ifndef SOAP_KERNEL_FUNCTION_H_INCLUDED
#define SOAP_KERNEL_FUNCTION_H_INCLUDED

#include <vector>

typedef std::vector<int> intvector;
typedef std::vector<std::vector<double> > dmatrix;
typedef std::vector<intvector> intvector_arr;
typedef std::vector<dmatrix> dmatrix_arr;
typedef std::vector<std::vector<std::vector<float> > > mol_repr;
typedef std::vector<mol_repr> mol_repr_arr;

mol_repr soap_repr_func(intvector atomic_no, dmatrix coords);
dmatrix soap_kernel_function(mol_repr_arr train_repr_arr, mol_repr_arr validate_repr_arr);
void *compute_kernel_matrix(void *kernel_data_ptr);

#endif // SOAP_KERNEL_FUNCTION_H_INCLUDED
