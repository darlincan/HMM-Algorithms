#include <stdio.h>
#include <stdlib.h>
#include "hmm.h"

int main( int argc, char* argv[] ) {
    
    // double** A = allocation_dynamique( 3, 3 );
    // A[0][0] = 0.45;
    // A[0][1] = 0.35;
    // A[0][2] = 0.20;

    // A[1][0] = 0.10;
    // A[1][1] = 0.50;
    // A[1][2] = 0.40;

    // A[2][0] = 0.15;
    // A[2][1] = 0.25;
    // A[2][2] = 0.60;

    // double** B = allocation_dynamique( 3, 2 );
    // B[0][0] = 1.0;
    // B[0][1] = 0.0;

    // B[1][0] = 0.5;
    // B[1][1] = 0.5;

    // B[2][0] = 0.0;
    // B[2][1] = 1.0;


    // double* Pi =  (double *)malloc(sizeof(double) * 3);
    // Pi[0] = 0.5;
    // Pi[1] = 0.3;
    // Pi[2] = 0.2;

    // int* observation = (int *)malloc(sizeof(int) * 5);
    // observation[0] = 1;
    // observation[1] = 2;
    // observation[2] = 2;
    // observation[3] = 1;
    // observation[4] = 1;

    // int* alphabet = (int *)malloc(sizeof(int) * 2);
    // alphabet[0] = 1; // a
    // alphabet[1] = 2; // b

    // int N, T, lenAlpha, epochs;
    // N = 3;
    // T = 5;
    // lenAlpha = 2;
    // epochs = 16;

    //baumWelch( N, T, lenAlpha, alphabet, observation, epochs , A, B, Pi );
    read_A_B( argv[1], argv[2], argv[3] );

    //display2D( A, 3, 3, "A" );
}