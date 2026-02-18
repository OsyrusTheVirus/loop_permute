#include "assignment1.h"

// Looping through IJP
void shpc_dgemm_IJP( int m, int n, int k,                                            
                    double *A, int rsA, int csA,                                
                    double *B, int rsB, int csB,                                
                    double *C, int rsC, int csC )
{
    for (int i = 0 ; i < m ; i++)
        for (int j = 0 ; j < n ; j++)
            for (int p = 0 ; p < k ; p++)
                *(C + i*rsC + p*csC) += *(A + i*rsA + j*csA) * *(B + j*rsB + p*csB);
}

// Looping through IPJ
void shpc_dgemm_IPJ( int m, int n, int k,                                            
                    double *A, int rsA, int csA,                                
                    double *B, int rsB, int csB,                                
                    double *C, int rsC, int csC )
{
    for (int i = 0 ; i < m ; i++)
        for (int p = 0 ; p < k ; p++)
            for (int j = 0 ; j < n ; j++)
                *(C + i*rsC + p*csC) += *(A + i*rsA + j*csA) * *(B + j*rsB + p*csB);
}   

// Looping through JIP
void shpc_dgemm_JIP( int m, int n, int k,                                            
                    double *A, int rsA, int csA,                                
                    double *B, int rsB, int csB,                                
                    double *C, int rsC, int csC )
{
    for (int j = 0 ; j < n ; j++)
        for (int i = 0 ; i < m ; i++)
            for (int p = 0 ; p < k ; p++)
                *(C + i*rsC + p*csC) += *(A + i*rsA + j*csA) * *(B + j*rsB + p*csB);
}   

// Looping through JPI
void shpc_dgemm_JPI( int m, int n, int k,                                            
                    double *A, int rsA, int csA,                                
                    double *B, int rsB, int csB,                                
                    double *C, int rsC, int csC )
{
    for (int j = 0 ; j < n ; j++)
        for (int p = 0 ; p < k ; p++)
            for (int i = 0 ; i < m ; i++)
                *(C + i*rsC + p*csC) += *(A + i*rsA + j*csA) * *(B + j*rsB + p*csB);
}   

// Looping through PIJ
void shpc_dgemm_PIJ( int m, int n, int k,                                            
                    double *A, int rsA, int csA,                                
                    double *B, int rsB, int csB,                                
                    double *C, int rsC, int csC )
{
    for (int p = 0 ; p < k ; p++)
        for (int i = 0 ; i < m ; i++)
            for (int j = 0 ; j < n ; j++)
                *(C + i*rsC + p*csC) += *(A + i*rsA + j*csA) * *(B + j*rsB + p*csB);
} 

// Looping through PJI
void shpc_dgemm_PJI( int m, int n, int k,                                            
                    double *A, int rsA, int csA,                                
                    double *B, int rsB, int csB,                                
                    double *C, int rsC, int csC )
{
    for (int p = 0 ; p < k ; p++)
        for (int j = 0 ; j < n ; j++)
            for (int i = 0 ; i < m ; i++)
                *(C + i*rsC + p*csC) += *(A + i*rsA + j*csA) * *(B + j*rsB + p*csB);
}