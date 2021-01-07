#! /usr/bin/env python3
# -*- Coding: UTF-8 -*-

r"""
Metropolis Monte Carlo simulation for 2-dimensional Ising Model.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse
import tqdm

spin = {'up': 1, 'down': -1}

length_default = 20
step_default = int(1E6)
temperature_default = 2.0
filename_default = ''
description = 'This program performs a Metropolis Monte Carlo simulation of 2-dimensional Ising model.'
epilog = '''
The default argumets are "-l 20 -s 1000000 -t 2.0".
And the default initial status is about the left half
is up-spin, and about the right half is down-spin.

Please pay attention that the Curie temperature 
of the analytic solution is 2/ln(1+sqrt(2)) = 2.269 J/k_B.
'''

def main():
    r'''
Execute the process.
'''
    global args
    figs= [None, None]
    axs = [None, None]
    Read_data()
    figs[0], axs[0] = Print_matrix('$T = %.2f\\, J/k_\\mathrm{B}$, Step = 0' % args.temperature)
    figs[0].savefig('Metropolis_Monte_Carlo_initial.png')
    print('Initial status: total magnetic moment is %5d m, total energy is %5d J.' % (Calc_magnet(), 
                                                                                      Calc_energy()))
    print('The initial status has been saved to "Metropolis_Monte_Carlo_initial.png".')
    for n in tqdm.trange(args.step):
        Operate_matrix_one_time()
    figs[1], axs[1] = Print_matrix('$T = %.2f\\, J/k_\\mathrm{B}$, Step = %d' % (args.temperature, args.step))
    figs[1].savefig('Metropolis_Monte_Carlo_result.png')
    print('Final   status: total magnetic moment is %5d m, total energy is %5d J.' % (Calc_magnet(), 
                                                                                      Calc_energy()))
    print('The final status has been saved to "Metropolis_Monte_Carlo_result.png".')
    print(u'Minimum energy: total magnetic moment is \u00b1%4d m, total energy is %5d J.' % 
        (board.shape[0] * board.shape[1], -2 * board.shape[0] * board.shape[1]))
    print('The final status has also been saved to "Metropolis_Monte_Carlo_result.txt".')
    Print_matrix_simple('Metropolis_Monte_Carlo_result.txt')
    plt.show()
    return

def Read_data():
    r'''
parse command line arguments, and initial the board.
'''
    parser = argparse.ArgumentParser(prog = '2d_Ising_Metropolis_Monte_Carlo_simulation', 
        description = description, epilog = epilog)
    parser.add_argument('-l', '--length', type = int, default = length_default, 
        help = 'number of atoms on each side.')
    parser.add_argument('-s', '--step', type = int, default = step_default, 
        help = 'number of total steps.')
    parser.add_argument('-t', '--temperature', type = float, default = temperature_default, 
        help = 'temperature (in unit J/k_Boltzmann)')
    parser.add_argument('-f', '--filename', type = str, default = filename_default, 
        help = 'name of file contains initial status, "-" means to read from stdin.')
    global args
    args = parser.parse_args()
    if args.length <= 1:
        raise ValueError('Error! Length should be at least 2, but it is %d.' % args.length)
    if args.step <= 0:
        raise ValueError('Error! Step should be positive, but it is %d.' % args.step)
    if args.temperature <= 0:
        raise ValueError('Error! Temperature should be positive, but it is %.2f.' % args.temperature)
    print('Length = %d, Step = %d, Temperature = %.2f, filename = \"%s\"' % (
        args.length, args.step, args.temperature, args.filename))
    global board
    board = np.empty((args.length, args.length), dtype = int)
    if args.filename == '':
        half = int(np.ceil(board.shape[1] / 2))
        if (board.shape[1] % 2):
            for i in range(half):
                for j in range(half):
                    board[i, j] = spin['up']
                for j in range(half, board.shape[1]):
                    board[i, j] = spin['down']
            for i in range(half, board.shape[0]):
                for j in range(half - 1):
                    board[i, j] = spin['up']
                for j in range(half - 1, board.shape[1]):
                    board[i, j] = spin['down']
        else:
            for i in range(board.shape[0]):
                for j in range(half):
                    board[i, j] = spin['up']
                for j in range(half, board.shape[1]):
                    board[i, j] = spin['down']
    else:
        if args.filename == '-':
            print('Reading initial status from stdin, please input:')
            board_tmp = np.loadtxt(sys.stdin)
        else:
            board_tmp = np.loadtxt(args.filename)
        if board_tmp.shape != board.shape:
            raise ValueError('Error! Shape of input (%d, %d) does not match ' + 
                'that specified by "length": (%d, %d).' % board_tmp.shape + board.shape)
        for i in range(board.shape[0]):
            for j in range(board.shape[1]):
                board[i, j] = spin['up'] if board_tmp[i, j] > 0 else spin['down']
    return

def Calc_one_energy(posi, posj):
    r'''
calculates energy at a specific position.
'''
    global board
    ene = 0
    if (posi >= board.shape[0]) or (posj >= board.shape[1]):
        raise ValueError("Error! Position (%d, %d) out of range (0, 0) to (%d, %d)." % (
            posi, posj, board.shape[0] - 1, board.shape[1] - 1))
    ene -= board[posi, posj] * board[posi - 1 if posi > 0 else board.shape[0] - 1, posj]
    ene -= board[posi, posj] * board[posi + 1 if posi < board.shape[0] - 1 else 0, posj]
    ene -= board[posi, posj] * board[posi, posj - 1 if posj > 0 else board.shape[1] - 1]
    ene -= board[posi, posj] * board[posi, posj + 1 if posj < board.shape[1] - 1 else 0]
    return ene

def Calc_energy():
    r'''
calculates energy of the entire board.
'''
    global board
    ene = 0
    for i in range(board.shape[0]):
        for j in range(board.shape[1]):
            ene += Calc_one_energy(i, j)
    ene //= 2 # Each element has been calculated twice.
    return ene

def Calc_magnet():
    r'''
calculates the moment of magnet of the entire board.
'''
    global board
    return board.sum()

def Operate_matrix_one_time():
    r'''
do one operation on the board
'''
    global board, args
    pos = np.random.randint(board.shape[0] * board.shape[1])
    posi, posj = pos // board.shape[0], pos % board.shape[1]
    delta_ene = Calc_one_energy(posi, posj) * (-2)
    if (delta_ene <= 0) or (np.exp(- delta_ene / args.temperature) >= np.random.rand()):
        board[posi, posj] = - board[posi, posj]
    return

def Print_matrix(title):
    r'''
draw the current status
'''
    global board
    fig, ax = plt.subplots()
    ax.imshow(board, cmap = plt.cm.hot, vmin = -1, vmax = 1)
    ax.set_xticks(list())
    ax.set_yticks(list())
    ax.set_title(title)
    return fig, ax

def Print_matrix_simple(outputname):
    r'''
out put the result to 'outputname' for continue running or other use.
'''
    global board
    with open(outputname, 'wt') as f:
        for i in range(board.shape[0]): print(* ['%2d' % board[i][j] for j in range(board.shape[1])], file = f)

if __name__ == '__main__':
    main()

