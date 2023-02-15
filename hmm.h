typedef struct ReturnTuple ReturnTuple;
struct ReturnTuple
{
    double prob;
    double** matrix;
};

typedef struct viterti viterbi_t;
struct viterbi
{
    char* chemin;
    double** matrix;
};

typedef struct hmm hmm;
struct hmm
{
    double** A;
    double** B;
    double* Pi;
    unsigned int number_of_states;
    unsigned int number_of_sym;
};


/* ALLOUER DYNAMIQUEMENT UNE MATRICE 2D */

double **allocation_dynamique(int row, int col);


/* ALLOUER DYNAMIQUEMENT UNE MATRICE 3D*/

double*** allocation_dynamique_3D(int X, int Y, int Z);


/* DESALLOUER DYNAMIQUEMENT UNE MATRICE 2D*/

void deallocate_dynamic_matrix(double **matrix, int row);


/* DESALOUER DYNAMIQUEMENT UNE MATRICE 3D */

void deallocate_dynamique_matrix_3D( double ***matrix, int row );


/* ################################# FONCTION D'AFFICHAGE ################################# */

void display1D( double* array, int N, char* nameArray );

void display2D( double** array, int N, int M, char* nameArray );

void display3D( double*** array, int X, int Y, int Z, char* nameArray );


/* ############## FONCTION RETURN INDEX OF SYMBOLE IN SEQUENCE OF OBSERVATION*/

int index_t( int obvCour, int* alphabet, int lenAlpha );


/* ##################### FORWARD ALGORITHM ##################### */

ReturnTuple forward(double** A, int n, double** B, double* Pi, int T, int* alphabet, int* observation , int lenAlpha);


/* ##################### BACKWARD ALGORITHM ##################### */

ReturnTuple backward(double** A, int n, double** B, double* Pi, int T, int* alphabet, int* observation , int lenAlpha);


/* ##################### MATRICE KSI ALGORITHM ##################### */

double*** ksi( int N, int T, double** alpha, double** beta, double** A, double** B, double prob, int* alphabet, int* observation , int lenAlpha);


/* ##################### GAMMA MATRIX ALGORITHM ##################### */

double** gammaM( int T, int N, double*** ksi );


/* ############################ APPRENTISSAGE ########################## */

/* ##################### RE-ESTIMATION VECTEUR ETATS INITIAUX ##################### */

double* rInit( double** GAMMA, int N );


/* ##################### RE-ESTIMATION MATRICE DE TRANSITION ##################### */

double** rTrans( double** GAMMA, double*** ksi, int N, int T );


/* ##################### RE-ESTIMATION MATRICE D'EMISSION ##################### */

double** rEmi( int N, int T, double** GAMMA, int* alphabet, int* observation , int lenAlpha );


/* ################################# FONCTION DE COPY ############################# */

double** arrayCopy2D( double** arrayInit, int N, int M );

double* arrayCopy1D( double* arrayInit, int N );


/* ################### ALGORITHME DE BAUM-WELCH ####################### */

void baumWelch( int N, int T, int lenAlpha, int* alphabet, int* observation, int epochs , double** A, double** B, double* Pi  );


/* ############################# LECTURE DES DONNEES DU MODELE ############################# */

hmm read_A_B( char* A_file, char* B_file, char* Pi_file );
