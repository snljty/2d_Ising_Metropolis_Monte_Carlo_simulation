Metropolis Monte Carlo simulation of a 2-dimensional Ising model.

The interactive energy between two particles are:
-1, if they have the same  spin.
 1, if they have different spin.

During each operation, a random particle is picked, and it is flipped at the probability:
                    1, if the total energy reduced.
 exp(-beta * delta_E), if the total energy raised.

At a temperature lower than Curie temperature, there will be a spontaneous magnetization phenomenon.
To observe this phenomenon, the simulation steps should be about at least 1E6.

The analytic solution of 2-dimensional Ising model is T_c = 2/ln(1+sqrt(2)) = 2.269 J/k_B.

Use command arguments to control this program. Use "-h" argument to see help.

Do not just double-click this program, run it within a command console instead.

before running, make sure the encoding of the console is set to UTF-8, 
which you can achieve by input "chcp 65001" and press <Enter> under Windows.

Under Windows OS, do not run under a cmd or powershell environment, as it does not support
the escape code by default. A cygwin console, a cmder console, or a "Windows terminal" is suggested.

Under Linux OS, compile with "-m" option while using gcc.

Looking forward to feedback by email to snljty@sina.com .

