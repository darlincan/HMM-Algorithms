from helper import matrix

def viterbi(A, B, Pi, O, Q) :
    """

    :param A: Matrice des probabilités de transition d'états
    :param B: Matrice des probabilités d'émission
    :param Pi: Matrice des probabilités des états initiaux
    :param O: Séquence d'observations de taille T
    :param Q: Ensemble des états de taille N
    :return: bestPathProb qui est la probabilité suivant le chemin optimal
             bestPath qui est le chemin optimal
    """
    n = len(Q)
    T = len(O)
    bestPath = []

    # Create a path probability matrix
    Viterbi = matrix(n, T, Q)
    # Create a bestpath matrix
    backpointer = matrix(n, T, Q)

    for s in Q:
        Viterbi[s][0] = Pi[s] * B[s][0]
        backpointer[s][0] = 0

    for t in range(1, T) :
        for s in Q:
            temp = []
            temp_dico = {}
            for e in Q:
                temp.append(Viterbi[e][t-1] * A[e][s] * B[s][t])
                temp_dico[e] = Viterbi[e][t-1] * A[e][s] * B[s][t]
            Viterbi[s][t] = max(temp)
            backpointer[s][t] = max(temp_dico, key=temp_dico.get)

        bestPath.append(backpointer[s][t])

    t = []
    dic = {}
    for s in Q:
        t.append(Viterbi[s][T-1])
        dic[s] = Viterbi[s][T-1]

    bestPathProb = max(t)
    bestPathPointer = max(dic, key=dic.get)
    bestPath.insert(0, bestPathPointer)
    bestPath = '-'.join(bestPath)

    return (bestPath, bestPathProb)


