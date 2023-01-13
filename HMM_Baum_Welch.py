from helper import *
from HMM_Forward import forward
from HMM_Backward import backward
def baumWelch(A, B, Pi, O, Q) :
    """
        :param A: Matrice des probabilités de transition d'états
        :param B: Matrice des probabilités d'émission
        :param Pi: Matrice des probabilités des états initiaux
        :param O: Séquence d'observations de taille T
        :param Q: Ensemble des états de taille N

    """
    n = len(Q)
    T = len(O)
    Xi = matrixVariable(n, Q, T)
    GAMMA = matrix(n, T, Q)

    FORWARD = forward(A, B, Pi, O, Q)
    BACKWARD = backward(A, B, Pi, O, Q)

    # E-Step
    for t in range(T-1):
        for i in Q:
            for j in Q:
                Xi[i][t][j] = varaiblesXi(t, i, j, FORWARD, BACKWARD, A, B, Q)
            GAMMA[i][t] = gamma(t, i, FORWARD, BACKWARD, A, B, Q)

    # M-STEP ( Ré-estimation des paramètres )
    # Ré-estimation de Pi et A
    for s in Q:
        Pi[s] = GAMMA[s][0]
        for j in Q:
            num = 0
            den = 0
            for i in range(T-1):
                num += Xi[s][i][j]
                den += GAMMA[s][i]
            A[s][j] += num / den

    # Ré-estimation de B
    for k in Q:
        for o in O:
            numer = denom = 0.0
            for p in range(T):
                if O(p) == o:
                    numer += GAMMA[k][p]
                else:
                    numer += 0
                denom += GAMMA[k][t]
        B[k][O.index(o)] = (numer) / denom