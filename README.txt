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

Before running, make sure the encoding of the console is set to UTF-8, 
which you can achieve by input "chcp 65001" and press <Enter> under Windows.

Under Windows OS, the previous version cannot be run on a cmd or powershell environment, as they 
do not support the escape code by default. A cygwin console, a cmder console, or a "Windows terminal" is recommended. In the current version, the Windows handle technics is used to print color, hence it
should has a normal performance in all situations.

Under Linux OS, compile with "-m" option while using gcc to link with functions in "math.h".

There was a Python version years before, but the drawing method is not elegant enough.
M ay be pushed later.

"step" may sometimes be greater than 2^32-1, which is usually the largest "unsigned int", but change 
it to "unsigned long long int" may cause a lot of troubles. Hence if you want more steps, just run a
shorter time and then continue running with "-f Metropolis_Monte_Carlo_result.txt" using the old 
result until reaching your expecting steps.

Looking forward to feedback by email to snljty@sina.com .

