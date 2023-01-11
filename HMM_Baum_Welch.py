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
