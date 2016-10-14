#include "local_similarity.h"

#include "neighbourhood.h"

#include <math.h>
#include <thread>

#define LOCAL_ZETA 1

typedef struct ls_data {
   Descriptor *desc_arr;
   double **ls_arr;
   double *diag;
   int desc_no;
   int nthreads;
   int id;
} LS_data;

double local_similarity(Power_spectrum *S1, Power_spectrum *S2, double *diag)
{
    int ntypes = ATOM_TYPES;
    double similarity_kernel[ntypes][ntypes];

    /* diagonal matrix */
    for (int i=0; i<ntypes; i++)
    {
        for (int j=0; j<ntypes; j++)
        {
            if (i == j) similarity_kernel[i][i] = diag[i];
            else similarity_kernel[i][j] = 0;
        }
    }

    double sum = 0;
    for (int i=0; i<ntypes; i++)
    {
        for (int j=0; j<ntypes; j++)
        {
            if (similarity_kernel[i][j] != 0)
            {
                double dot = dot_prod(S1[i],S2[j]);
                sum += similarity_kernel[i][j] * pow(dot,LOCAL_ZETA);
            }
        }
    }
    return sum;
}

void *compute_lsa(void *ls_data_ptr)
{
   LS_data *lsdp = (LS_data *) ls_data_ptr;
   for (int i=lsdp->id; i<lsdp->desc_no; i += lsdp->nthreads)
   {
        for (int atom_idx=0; atom_idx<MAX_TOTAL; atom_idx++)
        {
            Power_spectrum *ps_arr = lsdp->desc_arr[i][atom_idx];
            lsdp->ls_arr[i][atom_idx] = local_similarity(ps_arr,ps_arr,lsdp->diag);
        }
   }
   pthread_exit(NULL);
}

double **create_local_similarity_array(Descriptor *desc_arr,int desc_no,double *diag)
{
    //threads
    int rc;
    int nthreads = std::thread::hardware_concurrency();
    pthread_t threads[nthreads];
    pthread_attr_t attr;
    void *status;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    double **ls_arr = (double **) malloc(desc_no*sizeof(double *));
    for (int mol_idx=0; mol_idx<desc_no; mol_idx++)
    {
        ls_arr[mol_idx] = (double *) malloc(MAX_TOTAL*sizeof(double));
    }

    LS_data lsd_arr[nthreads];
    for(int i=0; i<nthreads; i++)
    {
        lsd_arr[i].desc_arr = desc_arr;
        lsd_arr[i].ls_arr = ls_arr;
        lsd_arr[i].diag = diag;
        lsd_arr[i].desc_no = desc_no;
        lsd_arr[i].nthreads = nthreads;
        lsd_arr[i].id = i;

        rc = pthread_create(&threads[i], &attr, compute_lsa, (void *)&lsd_arr[i]);
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
    for (int mol_idx=0; mol_idx<desc_no; mol_idx++)
    {
        for (int atom_idx=0; atom_idx<MAX_TOTAL; atom_idx++)
        {
            Power_spectrum *ps_arr = desc_arr[mol_idx][atom_idx];
            ls_arr[mol_idx][atom_idx] = local_similarity(ps_arr,ps_arr,diag);
        }
    }
*/

    return ls_arr;
}

int free_ls_arr(double **ls_arr, int desc_no)
{
    for (int mol_idx=0; mol_idx<desc_no; mol_idx++)
        free(ls_arr[mol_idx]);

    free(ls_arr);
    return 0;
}

