from HMM_Forward import forward
from HMM_Viterbi import viterbi
from HMM_Backward import backward

A = {
        'h' : {
                'h': 0.6,
                'c': 0.4
              },

        'c' : {
                'h': 0.5,
                'c': 0.5
                }
    }

B = {
        'c' : [0.5, 0.4, 0.1],
        'h' : [0.2, 0.4, 0.4]
    }

Pi = {
        'h' : 0.8,
        'c' : 0.2
     }

O = [1, 2, 3]

Q = ['h', 'c']

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    #prob = forward(A, B, Pi, O, Q)
    #path, prob = viterbi(A, B, Pi, O, Q)
    prob = backward(A, B, Pi, O, Q)

    #print('Chemin optimal : bestPath = {}'.format(path))
    #print('Probabilité du chemin optimal : P(O) = {}'.format(prob))

    print('P(O)1 = {}'.format(prob))


