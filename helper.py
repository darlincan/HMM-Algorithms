def matrix(row, col, Q) :
    M = {}

    for i in range(row):
        temp = []
        for j in range(col):
            temp.append(0)
        M[Q[i]] = temp
    return M

def varaiblesXi(t, i, j, alpha, beta, A, B, Q) :
    """
    la variable Xi rattachée aux index i et j à l'instant t est la probabilité d'avoir emprunté un chemin
    par lequel de l'instant t à l'instant t+1 il y'a une transition d'état ei -> ej (P(qt = ei, qt+1 = ej | O, modèle)
    :param t: l'instant t de la variableXi
    :param i: l'index i (l'état qi)
    :param j: l'index j (l'état qj)
    :param alpa: la matrice des variables forward
    :param beta: la matrice des variables backward
    :param A: la matrice des probabilités de transition d'états
    :param B: la matrice des probabilités d'émission
    :return: variableXi(i, j) à 'linstant t
    """
    n = len(Q)

    numerateur = alpha[i][t-1] * A[i][j] * B[j][t] * beta[j][t]
    denominateur = 0

    for k in Q:
        for l in Q:
            denominateur += alpha[k][l] * A[k][l] * B[l][t] * beta[l][t]

    return numerateur / denominateur


def gamma(t, i, alpha, beta, A, B, Q):
    Gamma = 0
    for j in Q:
         Gamma += varaiblesXi(t, i, j, alpha, beta, A, B, Q)

    return Gamma






M = matrix(2, 3, ['h','c'])
dico = {'h': [5, 0, 0], 'c': [0, 0, 0]}

for i in range(1, -1, -1):
    print(i)

#print(dico['h'][0])