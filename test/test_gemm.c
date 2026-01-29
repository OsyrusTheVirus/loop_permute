#include "assignment1.h"

int test_gemm( int nrepeats, int first, int last, int inc)
{

    int size, irep;
    int m, n, k;
    int csA, csB, csC;
    int rsA, rsB, rsC;

    double *A, *B, *C, *Cref, *Cold;

    double done = 1.0;


    double t_ref = DBL_MAX;
    double t     = DBL_MAX;
    double t_start; 

    double gflops_ref, gflops_IJP, gflops_IPJ, 
                       gflops_JIP, gflops_JPI, 
                       gflops_PIJ, gflops_PJI ;

    double diff, maxdiff = 0.0;

    printf( "%% --------- DGEMM --------- \n"); 
    printf( "data_dgemm_ref " );
    printf( "( %4lu, 1:5 ) = [ %5lu %5lu %5lu %8.2f %15.4e ];\n",
                         ( unsigned long )( last - first)/inc + 1,
                         ( unsigned long )0,
                         ( unsigned long )0,
                         ( unsigned long )0, 0.0, 0.0 );
    printf( "data_dgemm_IJP " );
    printf( "( %4lu, 1:5 ) = [ %5lu %5lu %5lu %8.2f %15.4e ];\n",
                         ( unsigned long )( last - first)/inc + 1,
                         ( unsigned long )0,
                         ( unsigned long )0,
                         ( unsigned long )0, 0.0, 0.0 );
    printf( "data_dgemm_IPJ " );
    printf( "( %4lu, 1:5 ) = [ %5lu %5lu %5lu %8.2f %15.4e ];\n",
                         ( unsigned long )( last - first)/inc + 1,
                         ( unsigned long )0,
                         ( unsigned long )0,
                         ( unsigned long )0, 0.0, 0.0 );
    printf( "data_dgemm_JPI " );
    printf( "( %4lu, 1:5 ) = [ %5lu %5lu %5lu %8.2f %15.4e ];\n",
                         ( unsigned long )( last - first)/inc + 1,
                         ( unsigned long )0,
                         ( unsigned long )0,
                         ( unsigned long )0, 0.0, 0.0 );
    printf( "data_dgemm_JIP " );
    printf( "( %4lu, 1:5 ) = [ %5lu %5lu %5lu %8.2f %15.4e ];\n",
                         ( unsigned long )( last - first)/inc + 1,
                         ( unsigned long )0,
                         ( unsigned long )0,
                         ( unsigned long )0, 0.0, 0.0 );
    printf( "data_dgemm_PIJ " );
    printf( "( %4lu, 1:5 ) = [ %5lu %5lu %5lu %8.2f %15.4e ];\n",
                         ( unsigned long )( last - first)/inc + 1,
                         ( unsigned long )0,
                         ( unsigned long )0,
                         ( unsigned long )0, 0.0, 0.0 );
    printf( "data_dgemm_PJI " );
    printf( "( %4lu, 1:5 ) = [ %5lu %5lu %5lu %8.2f %15.4e ];\n",
                         ( unsigned long )( last - first)/inc + 1,
                         ( unsigned long )0,
                         ( unsigned long )0,
                         ( unsigned long )0, 0.0, 0.0 );

    // Iterating through the matrix sizes.
    for ( size=last; size>= first; size-=inc )
    {
        /* we will only time cases where all three matrices are square */
        m = n = k = size;
        csA = m; csB = k; csC = m;

        rsA = rsB = rsC = 1;

        A = ( double * ) malloc( csA * k * sizeof( double ) );
        B = ( double * ) malloc( csB * n * sizeof( double ) );
        C = ( double * ) malloc( csC * n * sizeof( double ) );
        Cold = ( double * ) malloc( csC * n * sizeof( double ) );
        Cref = ( double * ) malloc( csC * n * sizeof( double ) );


        bli_drandm( 0, BLIS_DENSE, m, k, A, rsA, csA);
        bli_drandm( 0, BLIS_DENSE, k, n, B, rsB, csB);
        bli_drandm( 0, BLIS_DENSE, m, n, Cold, rsC, csC);


        for ( irep=0; irep<nrepeats; irep++ )
        {
            memcpy( Cref, Cold, csC * n * sizeof( double ) );

            t_start = bli_clock();
        
            bli_dgemm( BLIS_NO_TRANSPOSE, BLIS_NO_TRANSPOSE,  
                        m, n, k, &done, 
                        A, rsA, csA, 
                        B, rsB, csB, 
                        &done, Cref, rsC, csC );    
            t_ref = bli_clock_min_diff( t_ref, t_start );
            
        }

        gflops_ref = 2.0 * m * n * k / ( t_ref * 1.0e9 );
        printf( "data_dgemm_ref");
        printf( "( %4lu, 1:5 ) = [ %5lu %5lu %5lu %8.2f %15.4e ];\n",
                ( unsigned long )(size - first)/inc + 1,
                ( unsigned long )m,
                ( unsigned long )k,
                ( unsigned long )n, gflops_ref, 0.0);

        // Timing IJP 
        for ( irep=0; irep<nrepeats; irep++ )
        {
            memcpy( C, Cold, csC * n * sizeof( double ) );

            t_start = bli_clock();
        
            shpc_dgemm_IJP( m, n, k, 
                        A, rsA, csA, 
                        B, rsB, csB, 
                        C, rsC, csC );  
            
            t = bli_clock_min_diff( t , t_start );
            
        }

        gflops_IJP = 2.0 * m * n * k / ( t * 1.0e9 );
      
        diff    = shpc_maxabsdiff( m, n, C, rsC, csC, Cref, rsC, csC );
        maxdiff = max ( diff, maxdiff );
 
        printf( "data_dgemm_IJP");
        printf( "( %4lu, 1:5 ) = [ %5lu %5lu %5lu %8.2f %15.4e ];\n",
                ( unsigned long )(size - first)/inc + 1,
                ( unsigned long )m,
                ( unsigned long )k,
                ( unsigned long )n, gflops_IJP, diff );

        // Timing IPJ
        t     = DBL_MAX; //reset time
        for ( irep=0; irep<nrepeats; irep++ )
        {
            memcpy( C, Cold, csC * n * sizeof( double ) );

            t_start = bli_clock();
        
            shpc_dgemm_IPJ( m, n, k, 
                        A, rsA, csA, 
                        B, rsB, csB, 
                        C, rsC, csC );  
            
            t = bli_clock_min_diff( t , t_start );
            
        }

        gflops_IPJ = 2.0 * m * n * k / ( t * 1.0e9 );
        
        diff    = shpc_maxabsdiff( m, n, C, rsC, csC, Cref, rsC, csC );

        printf( "data_dgemm_IPJ");
        printf( "( %4lu, 1:5 ) = [ %5lu %5lu %5lu %8.2f %15.4e ];\n",
                ( unsigned long )(size - first)/inc + 1,
                ( unsigned long )m,
                ( unsigned long )k,
                ( unsigned long )n, gflops_IPJ, diff );

        // Timing JPI
        t     = DBL_MAX; //reset time
        for ( irep=0; irep<nrepeats; irep++ )
        {
            memcpy( C, Cold, csC * n * sizeof( double ) );

            t_start = bli_clock();
        
            shpc_dgemm_JPI( m, n, k, 
                        A, rsA, csA, 
                        B, rsB, csB, 
                        C, rsC, csC );  
            
            t = bli_clock_min_diff( t , t_start );
            
        }

        gflops_JPI = 2.0 * m * n * k / ( t * 1.0e9 );
        
        diff    = shpc_maxabsdiff( m, n, C, rsC, csC, Cref, rsC, csC );

        printf( "data_dgemm_JPI");
        printf( "( %4lu, 1:5 ) = [ %5lu %5lu %5lu %8.2f %15.4e ];\n",
                ( unsigned long )(size - first)/inc + 1,
                ( unsigned long )m,
                ( unsigned long )k,
                ( unsigned long )n, gflops_JPI, diff );

        // Timing JIP
        t     = DBL_MAX; //reset time
        for ( irep=0; irep<nrepeats; irep++ )
        {
            memcpy( C, Cold, csC * n * sizeof( double ) );

            t_start = bli_clock();
        
            shpc_dgemm_JIP( m, n, k, 
                        A, rsA, csA, 
                        B, rsB, csB, 
                        C, rsC, csC );  
            
            t = bli_clock_min_diff( t , t_start );
            
        }

        gflops_JIP = 2.0 * m * n * k / ( t * 1.0e9 );
        
        diff    = shpc_maxabsdiff( m, n, C, rsC, csC, Cref, rsC, csC );

        printf( "data_dgemm_JIP");
        printf( "( %4lu, 1:5 ) = [ %5lu %5lu %5lu %8.2f %15.4e ];\n",
                ( unsigned long )(size - first)/inc + 1,
                ( unsigned long )m,
                ( unsigned long )k,
                ( unsigned long )n, gflops_JIP, diff );

        // Timing PIJ
        t     = DBL_MAX; //reset time
        for ( irep=0; irep<nrepeats; irep++ )
        {
            memcpy( C, Cold, csC * n * sizeof( double ) );

            t_start = bli_clock();

            shpc_dgemm_PIJ( m, n, k, 
                A, rsA, csA, 
                B, rsB, csB, 
                C, rsC, csC );  
    
            t = bli_clock_min_diff( t , t_start );
    
        }

        gflops_PIJ = 2.0 * m * n * k / ( t * 1.0e9 );

        diff    = shpc_maxabsdiff( m, n, C, rsC, csC, Cref, rsC, csC );

        printf( "data_dgemm_PIJ");
        printf( "( %4lu, 1:5 ) = [ %5lu %5lu %5lu %8.2f %15.4e ];\n",
                ( unsigned long )(size - first)/inc + 1,
                ( unsigned long )m,
                ( unsigned long )k,
                ( unsigned long )n, gflops_PIJ, diff );


        // Timing PJI 
        t     = DBL_MAX; //reset time
        for ( irep=0; irep<nrepeats; irep++ )
        {
            memcpy( C, Cold, csC * n * sizeof( double ) );

            t_start = bli_clock();
        
            shpc_dgemm_PJI( m, n, k, 
                        A, rsA, csA, 
                        B, rsB, csB, 
                        C, rsC, csC );  
            
            t = bli_clock_min_diff( t , t_start );
            
        }

        gflops_PJI = 2.0 * m * n * k / ( t * 1.0e9 );
        
        diff    = shpc_maxabsdiff( m, n, C, rsC, csC, Cref, rsC, csC );

        printf( "data_dgemm_PJI");
        printf( "( %4lu, 1:5 ) = [ %5lu %5lu %5lu %8.2f %15.4e ];\n",
                ( unsigned long )(size - first)/inc + 1,
                ( unsigned long )m,
                ( unsigned long )k,
                ( unsigned long )n, gflops_PJI, diff );


        free( A );
        free( B );
        free( C );
        free( Cold );
        free( Cref );

    }


}

