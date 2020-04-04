import matplotlib.pyplot as plt
import matplotlib
import numpy as np


def show_mat(mat):
    vmax = np.max(mat[:])
    im = plt.imshow(mat, cmap='Blues', interpolation='nearest', vmax=vmax * 1.2)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            text = plt.text(j, i, '%.2f' % mat[i, j],
                            ha="center", va="center", color="black", fontsize=15)
    im.axes.get_xaxis().set_visible(False)
    im.axes.get_yaxis().set_visible(False)
    plt.tight_layout()
    plt.show()


def mult_seq(seq, pwm, height):
    onehot = make_onehot(seq, height)
    mult = np.multiply(pwm, onehot)
    response = np.sum(np.max(mult, 0))
    return response


def make_onehot(seq, height):
    onehot = np.zeros((height, len(seq)))
    onehot[np.array(seq), np.arange(len(seq))] = 1
    return onehot


HEIGHT = 4
WIDTH = 8
SNR = 1.5

long_seq = np.random.randint(0, HEIGHT, (100,), np.int)
target_seq_onehot = make_onehot(long_seq[20: 20 + WIDTH], HEIGHT)

pwm = np.random.rand(HEIGHT, WIDTH)
pwm += SNR * target_seq_onehot
pwm = pwm / np.sum(pwm, 0)
pwm = np.log(pwm / 0.25)

##########################
# RESPONSE
##########################
responses = np.array([mult_seq(long_seq[i:i + WIDTH], pwm, HEIGHT) for i in range(len(long_seq) - WIDTH)])
plt.fill_between(np.arange(len(responses)), responses, min(responses) - 1, color=[0.5, 0, 0])
plt.show()


target_mult = np.multiply(pwm, target_seq_onehot)
show_mat(pwm)
show_mat(target_seq_onehot)
show_mat(target_mult)
print(responses[20])