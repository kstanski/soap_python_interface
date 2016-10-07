#ifndef SOAP_KERNEL_FUNCTION_H_INCLUDED
#define SOAP_KERNEL_FUNCTION_H_INCLUDED

//#include "descriptor.h"
#include <vector>

typedef std::vector<int> intvector;
typedef std::vector<std::vector<double> > dmatrix;
typedef std::vector<intvector> intvector_arr;
typedef std::vector<dmatrix> dmatrix_arr;
typedef std::vector<std::vector<std::vector<float> > > mol_repr;
typedef std::vector<mol_repr> mol_repr_arr;

mol_repr_arr soap_repr_func(intvector_arr atomic_no, dmatrix_arr coords);
dmatrix soap_kernel_function(mol_repr_arr train_repr_arr, mol_repr_arr validate_repr_arr);

#endif // SOAP_KERNEL_FUNCTION_H_INCLUDED
