#include "structural_similarity.h"

#include "descriptor.h"
#include "local_similarity.h"
#include "neighbourhood.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <math.h>
#include <thread>

typedef struct ss_data {
   Descriptor *desc_arr;
   double **ls_arr;
   double *ss_arr;
   double *diag;
   int desc_no;
   int nthreads;
   int id;
} SS_data;

double structural_similarity(Descriptor A, Descriptor B, double *LSA, double *LSB, double *diag)
{
    namespace bnu = boost::numeric::ublas;
    bnu::matrix<double> C(MAX_TOTAL, MAX_TOTAL);
    bnu::matrix<double> Sink(MAX_TOTAL, MAX_TOTAL);

    for (int i=0; i<MAX_TOTAL; i++)
    {
        for (int j=0; j<MAX_TOTAL; j++)
        {
            double ls = local_similarity(A[i],B[j],diag)/sqrt(LSA[i]*LSB[j]);    //normalised
            C(i,j) = ls;
            Sink(i,j) = exp((ls-1)/REG_PARAM);
        }
    }

    /* Sinkhorn algorithm */
    bnu::vector<double> v(MAX_TOTAL), u(MAX_TOTAL), en(MAX_TOTAL), temp(MAX_TOTAL);
    for (int i=0; i<MAX_TOTAL; i++) en(i) = 1.0/MAX_TOTAL;
    v = en;
    for (int counter=0; counter<ITERATIONS_NO; counter++)
    {
        temp = prod(Sink,v);
        for (int i=0; i<MAX_TOTAL; i++) u(i) = en(i)/temp(i);
        temp = prod(trans(Sink),u);
        for (int i=0; i<MAX_TOTAL; i++) v(i) = en(i)/temp(i);
    }

    for (int i=0; i<MAX_TOTAL; i++)
    {
        for (int j=0; j<MAX_TOTAL; j++)
        {
            Sink(i,j) *= u(i)*v(j);
        }
    }

    Sink = prod(trans(Sink),C);
    double trace = 0;
    for (int i=0; i<MAX_TOTAL; i++) trace += Sink(i,i);

    return trace;
}

void *compute_ssa(void *ss_data_ptr)
{
   SS_data *ssdp = (SS_data *) ss_data_ptr;
   for (int i=ssdp->id; i<ssdp->desc_no; i += ssdp->nthreads)
   {
        Descriptor mol_desc = ssdp->desc_arr[i];
        double *LSA = ssdp->ls_arr[i];
        ssdp->ss_arr[i] = structural_similarity(mol_desc,mol_desc,LSA,LSA,ssdp->diag);
   }
   pthread_exit(NULL);
}

double *create_structural_similarity_array(Descriptor *desc_arr, double **ls_arr, int desc_no, double *diag)
{
    //threads
    int rc;
    int nthreads = std::thread::hardware_concurrency();
    pthread_t threads[nthreads];
    pthread_attr_t attr;
    void *status;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    double *ss_arr = (double *) malloc(desc_no*sizeof(double));

    SS_data ssd_arr[nthreads];
    for(int i=0; i<nthreads; i++)
    {
        ssd_arr[i].desc_arr = desc_arr;
        ssd_arr[i].ls_arr = ls_arr;
        ssd_arr[i].ss_arr = ss_arr;
        ssd_arr[i].diag = diag;
        ssd_arr[i].desc_no = desc_no;
        ssd_arr[i].nthreads = nthreads;
        ssd_arr[i].id = i;

        rc = pthread_create(&threads[i], &attr, compute_ssa, (void *)&ssd_arr[i]);
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
        Descriptor mol_desc = desc_arr[mol_idx];
        double *LSA = ls_arr[mol_idx];
        ss_arr[mol_idx] = structural_similarity(mol_desc,mol_desc,LSA,LSA,diag);
    }
*/
    return ss_arr;
}

int free_ss_arr(double *ss_arr)
{
    free(ss_arr);
    return 0;
}
