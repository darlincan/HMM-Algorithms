from helper import matrix

def forward(A, B, Pi, O, Q) :
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
    # Initialiser la martice forward de taille NxT
    Forward = matrix(n, T, Q)
    print('Matrice Forward initialisée : {} \n'.format(Forward))

    print('INITIALISATION \n')
    for s in Q:
        Forward[s][0] = Pi[s] * B[s][0]
        print('Pour {} :\n Forward[{},{}] = {} * {} = {}'.format(s, s, 0, Pi[s], B[s][0], Forward[s][0]))

    print('\n')

    print('RECURSION \n')
    for t in range(1,T):
        print('Pour t = {} :\n'.format(t))
        for s in Q:
            print('Pour {}:\n'.format(s))
            for e in Q:
                Forward[s][t] += (Forward[e][t-1] * A[e][s]) * B[s][t]

    for s in Q:
        prob += Forward[s][T-1]


    return prob