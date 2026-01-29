
/* Computes C = AB + C */
void shpc_dgemm_IJP( int m, int n, int k, 
                    double *A, int rsA, int csA, 
                    double *B, int rsB, int csB, 
                    double *C, int rsC, int csC );

void shpc_dgemm_IPJ( int m, int n, int k, 
                    double *A, int rsA, int csA, 
                    double *B, int rsB, int csB, 
                    double *C, int rsC, int csC );
                
void shpc_dgemm_JPI( int m, int n, int k, 
                    double *A, int rsA, int csA, 
                    double *B, int rsB, int csB, 
                    double *C, int rsC, int csC );
                
void shpc_dgemm_JIP( int m, int n, int k, 
                    double *A, int rsA, int csA, 
                    double *B, int rsB, int csB, 
                    double *C, int rsC, int csC );
                
void shpc_dgemm_PJI( int m, int n, int k, 
                    double *A, int rsA, int csA, 
                    double *B, int rsB, int csB, 
                    double *C, int rsC, int csC );
                
void shpc_dgemm_PIJ( int m, int n, int k, 
                    double *A, int rsA, int csA, 
                    double *B, int rsB, int csB, 
                    double *C, int rsC, int csC );
                
                
