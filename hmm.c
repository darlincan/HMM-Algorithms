#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hmm.h"

#define TAILLE_MAX 500


/* ALLOUER DYNAMIQUEMENT UNE MATRICE 2D */

double **allocation_dynamique(int row, int col)
{
    double **ret_val;
    int i;

    ret_val = (double **)malloc(sizeof(double *) * row);
    if (ret_val == NULL) {
        perror("memory allocation failure");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < row; ++i) {
        ret_val[i] = (double *)malloc(sizeof(double) * col);
        if (ret_val[i] == NULL) {
            perror("memory allocation failure");
            exit(EXIT_FAILURE);
        }
    }

    for( int i = 0; i < row; i++ ) {
        for( int j = 0; j < col; j++ ) {
            ret_val[i][j] = 0.0;
        }
    }

    return ret_val;
}


/* ALLOUER DYNAMIQUEMENT UNE MATRICE 3D*/

double*** allocation_dynamique_3D(int X, int Y, int Z) {

   double*** A = (double***)malloc(X * sizeof(double**));
 
    if (A == NULL)
    {
        fprintf(stderr, "Out of memory");
        exit(0);
    }
 
    for (int i = 0; i < X; i++)
    {
        A[i] = (double**)malloc(Y * sizeof(double*));
 
        if (A[i] == NULL)
        {
            fprintf(stderr, "Out of memory");
            exit(0);
        }
 
        for (int j = 0; j < Y; j++)
        {
            A[i][j] = (double*)malloc(Z * sizeof(double));
            if (A[i][j] == NULL)
            {
                fprintf(stderr, "Out of memory");
                exit(0);
            }
        }
    }
 
    // assign values to the allocated memory
    for (int i = 0; i < X; i++)
    {
        for (int j = 0; j < Y; j++)
        {
            for (int k = 0; k < Z; k++) {
                A[i][j][k] = 0.0;
            }
        }
    }

    return A;
}


/* DESALLOUER DYNAMIQUEMENT UNE MATRICE 2D*/

void deallocate_dynamic_matrix(double **matrix, int row)
{
    int i;

    for (i = 0; i < row; ++i)
    {
        free(matrix[i]);
    }
    free(matrix);
}


/* DESALOUER DYNAMIQUEMENT UNE MATRICE 3D */

void deallocate_dynamique_matrix_3D( double ***matrix, int row ) {
    int i,j;

    for( i = 0; i < row; ++i ) {
        for ( j = 0; j < row; ++j )
        {
            free( matrix[i][j] );
        }
        
        free( matrix[i] );
    }

    free( matrix );
}


/* ################################# FONCTION D'AFFICHAGE ################################# */

void display1D( double* array, int N, char* nameArray ) {
    for( int i = 0; i < N; i++ ) {
        printf( "%s[%d] = %f \t", nameArray, i ,array[i] ) ;
    }
    printf( "\n" );
}

void display2D( double** array, int N, int M, char* nameArray ) {
    for( int i = 0; i < N; i++ ) {
        for( int j = 0; j < M; j++ ) {
            printf( "%s[%d][%d] = %f \t", nameArray, i, j, array[i][j] );
        }
        printf( "\n" );
    }
    printf( "\n" );
}

void display3D( double*** array, int X, int Y, int Z, char* nameArray ) {
    for( int i = 0; i < X; i++ )
    {
        for( int j = 0; j < Y; j++ )
        {
            for( int k = 0; k < Z; k++ ) {
                printf( "%s[%d][%d][%d] = %f \t", nameArray, i, j, k, array[i][j][k] );
            }
            printf( "\n" );
        }
        printf( "\n" );
    }
    printf( "\n" );
}



/* ############## FONCTION RETURN INDEX OF SYMBOLE IN SEQUENCE OF OBSERVATION*/
int index_t( int obvCour, int* alphabet, int lenAlpha ) {
    int j = 0;
    int index_t = -1;
    
    while ( index_t == -1 || j < lenAlpha )
    {
        if( obvCour == alphabet[j] ) {
            index_t = j;
        }
        j++;
    }
    return index_t;
}


/* ##################### FORWARD ALGORITHM ##################### */

ReturnTuple forward(double** A, int n, double** B, double* Pi, int T, int* alphabet, int* observation , int lenAlpha) {

    double prob = 0.0;
    /*Initialiser la matrice Forward de taille NxT*/
    double** FORWARD = allocation_dynamique( n, T );

    /* Déclarer la structure de retour */
    ReturnTuple val;

    for( int i = 0; i < n; i++ ) {
        FORWARD[i][0] = Pi[i] * B[i][index_t( observation[0], alphabet, lenAlpha )];
    }

    int t = 1;
    while ( t < T )
    {
        for( int j = 0; j < n; j++ ) {
            for( int i = 0; i < n; i++ ) {
                FORWARD[j][t] += FORWARD[i][t-1] * A[i][j];
            }
        FORWARD[j][t] *= B[j][index_t( observation[t], alphabet, lenAlpha )];
        }
        t=t+1;
    }
    
    for( int i = 0; i < n; i++ ) {
        prob += FORWARD[i][T-1];
    }

    val.matrix = FORWARD;
    val.prob = prob;
    
    // deallocate_dynamic_matrix( FORWARD, n );
    return val;
}


/* ##################### BACKWARD ALGORITHM ##################### */

ReturnTuple backward(double** A, int n, double** B, double* Pi, int T, int* alphabet, int* observation , int lenAlpha) {

    double prob = 0.0;
    /*Initialiser la matrice Forward de taille NxT*/
    double** BACKWARD = allocation_dynamique( n, T );

    /* Déclarer la structure de retour */
    ReturnTuple val;

    for( int i = 0; i < n; i++ ) {
        BACKWARD[i][T-1] = 1.0;
    }

    int t = T-1;

    while ( t >= 1 )
    {
        for( int j = 0; j < n; j++ ) {
            for( int i = 0; i < n; i++ ) {
                BACKWARD[j][t-1] += BACKWARD[i][t] * A[j][i] * B[i][index_t( observation[t], alphabet, lenAlpha )];
            }
        }
        t=t-1;
    }

    for( int i = 0; i < n; i++ ) {
        prob += Pi[i] * B[i][0] * BACKWARD[i][0];
    }

    val.matrix = BACKWARD;
    val.prob = prob;

    //deallocate_dynamic_matrix( BACKWARD, n );
    return val;
}



/* ##################### MATRICE KSI ALGORITHM ##################### */

double*** ksi( int N, int T, double** alpha, double** beta, double** A, double** B, double prob, int* alphabet, int* observation , int lenAlpha) {

    double*** ksi = allocation_dynamique_3D( N, N, T );
    
    for( int t = 0; t < T-1; t++ ) {

        for( int i = 0; i < N; i++ ) {

            for( int j = 0; j < N; j++ ) {

                if( prob == 0.0 ) {
                    ksi[i][j][t] = 0.0;
                }
                else {
                    ksi[i][j][t] = ( alpha[i][t] * A[i][j] * B[j][index_t( observation[t+1], alphabet, lenAlpha )] * beta[j][t+1] ) / prob;
                }
            }
        }
    }

    display3D( ksi, N, N, T, "ksi" );

    return ksi;
} 


/* ##################### GAMMA MATRIX ALGORITHM ##################### */
double** gammaM( int T, int N, double*** ksi ) {

    double** GAMMA = allocation_dynamique( N, T );

    for( int t = 0; t < T-1; t++ ) {

        for( int i = 0; i < N; i++ ) {

            for( int j = 0; j < N; j++ ) {

                GAMMA[i][t] += ksi[i][j][t];
            }
        }          
    }

    display2D( GAMMA, N, T, "GAMMA" );

    return GAMMA;
}



/* ############################ APPRENTISSAGE ########################## */

/* ##################### RE-ESTIMATION VECTEUR ETATS INITIAUX ##################### */

double* rInit( double** GAMMA, int N ) {

    double* pi = (double *)malloc( sizeof(double) * N );

    for( int i = 0; i < N; i++ ) {

        pi[i] = GAMMA[i][0];

    }

    display1D( pi, N, "Pi*" );

    return pi;
}


/* ##################### RE-ESTIMATION MATRICE DE TRANSITION ##################### */

double** rTrans( double** GAMMA, double*** ksi, int N, int T ) {

    double** trans_new = allocation_dynamique( N, N );
    double ksi_sum, gamma_sum;

    for( int i = 0; i < N; i++ ) {

        for( int j = 0; j < N; j++ ) {

            ksi_sum = 0.0;
            gamma_sum = 0.0;

            for( int t = 0; t < T-1; t++ ) {

                ksi_sum +=  ksi[i][j][t];
                gamma_sum += GAMMA[i][t];
            }

            if( gamma_sum == 0.0 ) {
                trans_new[i][j] = 0.0;
            }
            else {
                trans_new[i][j] = ksi_sum / gamma_sum;
                //printf( "trans_new[%d][%d] = %f / %f = %f\n", i, j, ksi_sum, gamma_sum, ksi_sum / gamma_sum );
            }
        }
    }

    display2D( trans_new, N, N, "A*" );

    return trans_new;
}


/* ##################### RE-ESTIMATION MATRICE D'EMISSION ##################### */

double** rEmi( int N, int T, double** GAMMA, int* alphabet, int* observation , int lenAlpha ) {

    double** emi_new = allocation_dynamique( N, T );
    double gamma_sum, gamma_num_sum;

    for( int i = 0; i < N; i++ ) {

        for( int j = 0; j < T; j++ ) {

            gamma_num_sum = 0.0;
            gamma_sum = 0.0;

            for( int t = 0; t < T-1; t++ ) {

                if( index_t( observation[t], alphabet, lenAlpha ) == j ) {

                    gamma_num_sum += GAMMA[i][t];

                }

                gamma_sum += GAMMA[i][t];
            }

            if( gamma_sum == 0.0 ) {
                emi_new[i][j] = 0.0;
            }
            else {
                emi_new[i][j] = gamma_num_sum / gamma_sum;
            }
        }
    }

    display2D( emi_new, N, T, "B*" );
    
    return emi_new;
}



/* ################################# FONCTION DE COPY ############################# */

double** arrayCopy2D( double** arrayInit, int N, int M ) {

    double** arraycopy = allocation_dynamique( N, M );

    for( int i = 0; i < N; i++ ) {
        
        for( int j = 0; j < M; j++ ) {

            arraycopy[i][j] = arrayInit[i][j];
        }
    }
    return arraycopy;
}


double* arrayCopy1D( double* arrayInit, int N ) {

    double* arraycopy = (double*)malloc( N * sizeof(double*) );

    for( int i = 0; i < N; i++ ) {

        arraycopy[i] = arrayInit[i];
    }
    return arraycopy;
}



/* ################### ALGORITHME DE BAUM-WELCH ####################### */

void baumWelch( int N, int T, int lenAlpha, int* alphabet, int* observation, int epochs , double** A, double** B, double* Pi  ) {

    int iter = 0;
    double prob;

    double** trans_new = arrayCopy2D( A, N, N );
    double** emi_new = arrayCopy2D( B, N, T );
    double* init_new = arrayCopy1D( Pi, N );

    double*** KSI;
    double** GAMMA;

    while ( iter < epochs )
    {
        printf( "Epoch = %d\n", iter );
        /* Calcul de Forward et Backward */
        double** FORWARD = forward( trans_new, N, emi_new, init_new, T, alphabet, observation, lenAlpha ).matrix;
        display2D( FORWARD, N, T, "alpha" );

        double** BACKWARD = backward( trans_new, N, emi_new, init_new, T, alphabet, observation, lenAlpha ).matrix;
        display2D( BACKWARD, N, T, "beta" );

        prob = forward( trans_new, N, emi_new, init_new, T, alphabet, observation, lenAlpha ).prob;
        printf( "Prob = %f\n", prob );


        /* Calcul de ksi et gamma */
        
        KSI = ksi( N, T, FORWARD, BACKWARD, trans_new, emi_new, prob, alphabet, observation, lenAlpha );
        GAMMA = gammaM( T, N, KSI );

        /* Ré-estimation des paramètres du modèle */

        init_new = rInit( GAMMA, N );
        trans_new = rTrans( GAMMA, KSI, N, T );
        emi_new = rEmi( N, T, GAMMA, alphabet, observation , lenAlpha );

        iter++;
    }
    
}

hmm read_A_B( char* A_file, char* B_file, char* Pi_file ) {

    hmm model;
    int i, N, T;
    FILE* file_A = NULL;
    FILE* file_B = NULL;
    FILE* file_Pi = NULL;

    char chaine[TAILLE_MAX] = "";
    double** A;
    double** B;
    double* Pi;

    char* ptr;

    /* ################################ LECTURE DE A (matrice de transitions d'états) #################################### */

    file_A = fopen( A_file, "r" );

    if( file_A != NULL ) {

        i = 0;
        while( fgets( chaine, TAILLE_MAX, file_A ) != NULL ) {

            if( i == 0 ) {

                N = atoi( chaine );
                model.number_of_states = N;
                i = 1;
            }

            if( i == 1 ) 
            {
                A = allocation_dynamique( N, N );

                for( int j = 0; j < N; j++ ) {

                    int k = 0;
                    
                    while( k < model.number_of_states ) {

                        if( fgets( chaine, TAILLE_MAX, file_A ) != NULL )
						{
                           A[j][k] = strtod( chaine, NULL );
                           k++;
						}
                    }
                }
            }

        }
        model.A = A;
        fclose( file_A );      
    }

    /* ################################ LECTURE DE B (matrice d'émission de probabilités) #################################### */

    file_B = fopen( B_file, "r" );

    char chaine_B[TAILLE_MAX] = "";

    if( file_B != NULL ) {

        i = 0;
        while( fgets( chaine_B, TAILLE_MAX, file_B ) != NULL ) {

            if( i == 0 ) {

                T = atoi( chaine_B );
                model.number_of_sym = T;
                i = 1;
            }

            if( i == 1 ) 
            {
                B = allocation_dynamique( model.number_of_states, T );

                for( int j = 0; j < N; j++ ) {

                    int k = 0;
                    
                    while( k < model.number_of_sym ) {

                        if( fgets( chaine_B, TAILLE_MAX, file_B ) != NULL )
						{
                           B[j][k] = strtod( chaine_B, &ptr );
                           k++;
						}
                    }
                }
            }

        }
        model.B = B;
        fclose( file_B );   

    }

    /* ################################ LECTURE DE Pi (vecteur de prob d'états initiaux) #################################### */

    file_Pi = fopen( Pi_file, "r" );

    char chaine_P[TAILLE_MAX] = "";

    if( file_Pi != NULL ) {

        i = 0;
        while( fgets( chaine_P, TAILLE_MAX, file_Pi ) != NULL ) {

            Pi =  (double *)malloc( sizeof(double) * model.number_of_states );

            Pi[i] = strtod( chaine_P, &ptr );

            i++;
        }

        model.Pi = Pi;
        fclose( file_Pi );   

    }

    return model;
}
