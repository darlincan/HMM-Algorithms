from test import matrix

def backward(A, B, Pi, O, Q) :
    """
    :param A: Matrice des probabilités de transition d'états
    :param B: Matrice des probabilités d'émission
    :param Pi: Matrice des probabilités des états initiaux
    :param O: Séquence d'observations de taille T
    :param Q: Ensemble des états de taille N
    :return: P(O) qui est la probabilité d'observation de la séquence O sachant le modèle

    """
    n = len(Q)
    T = len(O)
    prob = 0
    # Initialiser la martice Backward de taille NxT
    Backward = matrix(n, T, Q)
    print('Matrice Backward initialisée : {} \n'.format(Backward))

    for s in Q:
        Backward[s][T-1] = 1

    for t in range(T-2, 0):
        print('Pour t = {}'.format(t))
        for s in Q:
            for e in Q:
                Backward[s][t] += (Backward[e][t+1] * A[e][s]) * B[s][t+1]
                print('Bacward')

    for s in Q:
        print('Val :'.format(Backward[s][0]))
        prob += Pi[s] * B[s][0] * Backward[s][0]

    return prob