import matplotlib.pyplot as plt
import matplotlib
import numpy as np


def show_mat(mat):
    vmax = np.max(mat[:])
    im = plt.imshow(mat, cmap='Blues', interpolation='nearest', vmax=vmax*1.2)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            text = plt.text(j, i, '%.2f' % mat[i, j],
                           ha="center", va="center", color="black", fontsize=55)
    im.axes.get_xaxis().set_visible(False)
    im.axes.get_yaxis().set_visible(False)
    plt.tight_layout()
    plt.show()


def mult_seq(seq, pwm, height):
    onehot = make_onehot(seq, height)
    mult = np.multiply(pwm, onehot)
    response = np.log(np.prod(np.max(mult, 0)))
    return response


def make_onehot(seq, height):
    onehot = np.zeros((height, len(seq)))
    onehot[np.array(seq), np.arange(len(seq))] = 1
    return onehot


height = 4
width = 8
long_seq = np.random.randint(0, height, (100,), np.int)
target_seq_onehot = make_onehot(long_seq[20: 20 + width], height)

pwm = np.random.rand(height, width)
pwm += 3 * target_seq_onehot
pwm = pwm / np.sum(pwm, 0)

responses = np.array([mult_seq(long_seq[i:i+width], pwm, height) for i in range(len(long_seq) - width)])
print(responses[20])
plt.fill_between(np.arange(len(responses)), responses, min(responses)-1, color=[0.5, 0, 0])
plt.show()


target_mult = np.multiply(pwm, target_seq_onehot)
show_mat(target_mult)
show_mat(target_seq_onehot)
show_mat(pwm)