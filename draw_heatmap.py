#! /usr/bin/env python3
# -*- Coding: UTF-8 -*-

r"""
Draw heat map from 'Metropolis_Monte_Carlo_result.txt'.
"""

import numpy as np
import matplotlib.pyplot as plt

a = np.loadtxt('Metropolis_Monte_Carlo_result.txt')
fig, ax = plt.subplots()
ax.imshow(a, plt.cm.hot)
ax.set_xticks(list())
ax.set_yticks(list())
fig.savefig('Metropolis_Monte_Carlo_result.png')
plt.show()

