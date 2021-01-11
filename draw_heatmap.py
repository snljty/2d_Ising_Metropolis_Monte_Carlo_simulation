#! /usr/bin/env python3
# -*- Coding: UTF-8 -*-

r"""
Draw heat map from 'Metropolis_Monte_Carlo_result.txt'.
"""

import numpy as np
import matplotlib.pyplot as plt

in_name  = 'Metropolis_Monte_Carlo_result.txt'
out_name = 'Metropolis_Monte_Carlo_result.png'

try:
    a = np.loadtxt('Metropolis_Monte_Carlo_result.txt')
except OSError:
    while True:
        print('Cannot found "{in_name:s}". Input file name now:'.format(in_name = in_name))
        try:
            a = np.loadtxt(input())
            break
        except OSError:
            ...
fig, ax = plt.subplots()
ax.imshow(a, plt.cm.hot)
ax.set_xticks(list())
ax.set_yticks(list())
fig.savefig(out_name)
print('Figure saved to "{out_name:s}"'.format(out_name = out_name))
plt.show()

