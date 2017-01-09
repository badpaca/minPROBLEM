# minPROBLEM
Solver for a minimization problem and smoother using iterative methods.

Solves the minimization problem f(x)=1/2x^TAx-b^Tx across the real
domain. A is SPD and large. I use various methods to solve this.

Then, I use some of the solver algorithms to smooth some rough data
iteratively. See comments in code for details.

C 2016, badpaca
