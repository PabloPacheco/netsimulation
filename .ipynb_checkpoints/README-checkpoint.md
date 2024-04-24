# Net Simulation

Python program based on the code from the book:
- Brebbia, C. A., & Ferrante, A. J. (2013). Computational hydraulics. Butterworth-Heinemann.

<img src = "BrebbiaBook.jpg" width = "333">

## Introduction

Net Simulation is a code that solves systems connected by one-dimensional elements, as shown in the figure:

<img src="netExample.svg" width="400">

where the one-dimensional element is represented as follows:

<img src="pipeElement.svg" width="300">

Each element of the system is mathematically related as shown below, using the local matrix of the element:

$$\begin{bmatrix}b_{k}^{i} \\ b_{j}^{i}\end{bmatrix}=k^{i}\begin{bmatrix}1 & -1 \\ -1 & 1\end{bmatrix}\begin{bmatrix}x_{k}^{i} \\ x_{j}^{i}\end{bmatrix}$$

where $x_{k}^{i}$ and $x_{j}^{i}$ are the unknowns, $b_{k}^{i}$, $b_{j}^{i}$ and $k^i$ are constants.
Each system of one-dimensional elements contains a certain number of nodes and elements. Considering that each element has an associated local matrix, a global matrix of the entire system can be obtained. This matrix is formed according to the connectivity of the elements. This global matrix is represented as follows:

$$ \mathbf{A} \mathbf{x} = \mathbf{b} $$

There are problems where the matrix $\mathbf{A}$ depends on the values of the unknowns stored in the vector $\mathbf{x}$. These cases correspond to nonlinear problems. Net Simulation allows solving this type of nonlinear problems, which are represented as follows:

$$ \mathbf{A}\left(\mathbf{x}\right) \mathbf{x} = \mathbf{b} $$ 

where the global matrix $\mathbf{A}\left(\mathbf{x}\right)$ is a function of the vector $\mathbf{x}$


## Solution Method

### Linear Systems

The solution of the linear system is carried out using the partition method. Depending on the type of physics of the problem to be solved and the type of boundary conditions of the system, the unknowns may be in the vector $\mathbf{x}$ or in the vector $\mathbf{b}$. For this reason, the global system of equations is partitioned as follows:

$$ \mathbf{A} \mathbf{x} = \mathbf{b} \quad \rightarrow \quad \left[\begin{matrix} \mathbf{A}_E & \mathbf{A}_{EF} \\  \mathbf{A}_{EF}^T & \mathbf{A}_F\end{matrix} \right]
\left[\begin{matrix} \mathbf{x}_E \\ \mathbf{x}_F \end{matrix} \right] = \left[\begin{matrix} \mathbf{b}_E \\ \mathbf{b}_F\end{matrix}\right] $$

where:
- $\mathbf{x}_F$: nodes where $x$ is unknown
- $\mathbf{x}_E$: nodes where $x$ is known
- $\mathbf{b}_F$: nodes where $b$ is known
- $\mathbf{b}_E$: nodes where $b$ is unknown

From the partitioned system, the following relationship is obtained:

$$ \mathbf{A}_{EF}^T \mathbf{x}_E + \mathbf{A}_F \mathbf{x}_F = \mathbf{b}_F$$

Organizing the equation conveniently for its solution:

$$\mathbf{A}_F \mathbf{x}_F = \mathbf{b}_F -  \mathbf{A}_{EF}^T \mathbf{x}_E $$

In this way, the vector $\mathbf{x}_F$ can be obtained by solving using some solver. Once $\mathbf{x}_F$ is calculated, the vector $\mathbf{b}_E$ can be computed:

$$ \mathbf{b}_E = \mathbf{A}_E \mathbf{x}_E + \mathbf{A}_{EF} \mathbf{x}_F $$ 

The solvers used in this code correspond to iterative solvers of [scipy sparse](https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html#module-scipy.sparse.linalg), which are shown below:

<p>Available iterative solvers:</p>
<table class="autosummary longtable table autosummary">
<tbody>
<tr class="row-odd"><td><p><a class="reference internal" href="generated/scipy.sparse.linalg.bicg.html#scipy.sparse.linalg.bicg" title="scipy.sparse.linalg.bicg"><code class="xref py py-obj docutils literal notranslate"><span class="pre">bicg</span></code></a>(A,&nbsp;b[,&nbsp;x0,&nbsp;tol,&nbsp;maxiter,&nbsp;M,&nbsp;callback,&nbsp;...])</p></td>
<td><p>Use BIConjugate Gradient iteration to solve <code class="docutils literal notranslate"><span class="pre">Ax</span> <span class="pre">=</span> <span class="pre">b</span></code>.</p></td>
</tr>
<tr class="row-even"><td><p><a class="reference internal" href="generated/scipy.sparse.linalg.bicgstab.html#scipy.sparse.linalg.bicgstab" title="scipy.sparse.linalg.bicgstab"><code class="xref py py-obj docutils literal notranslate"><span class="pre">bicgstab</span></code></a>(A,&nbsp;b,&nbsp;*[,&nbsp;x0,&nbsp;tol,&nbsp;maxiter,&nbsp;M,&nbsp;...])</p></td>
<td><p>Use BIConjugate Gradient STABilized iteration to solve <code class="docutils literal notranslate"><span class="pre">Ax</span> <span class="pre">=</span> <span class="pre">b</span></code>.</p></td>
</tr>
<tr class="row-odd"><td><p><a class="reference internal" href="generated/scipy.sparse.linalg.cg.html#scipy.sparse.linalg.cg" title="scipy.sparse.linalg.cg"><code class="xref py py-obj docutils literal notranslate"><span class="pre">cg</span></code></a>(A,&nbsp;b[,&nbsp;x0,&nbsp;tol,&nbsp;maxiter,&nbsp;M,&nbsp;callback,&nbsp;...])</p></td>
<td><p>Use Conjugate Gradient iteration to solve <code class="docutils literal notranslate"><span class="pre">Ax</span> <span class="pre">=</span> <span class="pre">b</span></code>.</p></td>
</tr>
<tr class="row-even"><td><p><a class="reference internal" href="generated/scipy.sparse.linalg.cgs.html#scipy.sparse.linalg.cgs" title="scipy.sparse.linalg.cgs"><code class="xref py py-obj docutils literal notranslate"><span class="pre">cgs</span></code></a>(A,&nbsp;b[,&nbsp;x0,&nbsp;tol,&nbsp;maxiter,&nbsp;M,&nbsp;callback,&nbsp;...])</p></td>
<td><p>Use Conjugate Gradient Squared iteration to solve <code class="docutils literal notranslate"><span class="pre">Ax</span> <span class="pre">=</span> <span class="pre">b</span></code>.</p></td>
</tr>
<tr class="row-odd"><td><p><a class="reference internal" href="generated/scipy.sparse.linalg.gmres.html#scipy.sparse.linalg.gmres" title="scipy.sparse.linalg.gmres"><code class="xref py py-obj docutils literal notranslate"><span class="pre">gmres</span></code></a>(A,&nbsp;b[,&nbsp;x0,&nbsp;tol,&nbsp;restart,&nbsp;maxiter,&nbsp;M,&nbsp;...])</p></td>
<td><p>Use Generalized Minimal RESidual iteration to solve <code class="docutils literal notranslate"><span class="pre">Ax</span> <span class="pre">=</span> <span class="pre">b</span></code>.</p></td>
</tr>
<tr class="row-even"><td><p><a class="reference internal" href="generated/scipy.sparse.linalg.lgmres.html#scipy.sparse.linalg.lgmres" title="scipy.sparse.linalg.lgmres"><code class="xref py py-obj docutils literal notranslate"><span class="pre">lgmres</span></code></a>(A,&nbsp;b[,&nbsp;x0,&nbsp;tol,&nbsp;maxiter,&nbsp;M,&nbsp;...])</p></td>
<td><p>Solve a matrix equation using the LGMRES algorithm.</p></td>
</tr>
<tr class="row-odd"><td><p><a class="reference internal" href="generated/scipy.sparse.linalg.minres.html#scipy.sparse.linalg.minres" title="scipy.sparse.linalg.minres"><code class="xref py py-obj docutils literal notranslate"><span class="pre">minres</span></code></a>(A,&nbsp;b[,&nbsp;x0,&nbsp;shift,&nbsp;tol,&nbsp;maxiter,&nbsp;M,&nbsp;...])</p></td>
<td><p>Use MINimum RESidual iteration to solve Ax=b</p></td>
</tr>
<tr class="row-even"><td><p><a class="reference internal" href="generated/scipy.sparse.linalg.qmr.html#scipy.sparse.linalg.qmr" title="scipy.sparse.linalg.qmr"><code class="xref py py-obj docutils literal notranslate"><span class="pre">qmr</span></code></a>(A,&nbsp;b[,&nbsp;x0,&nbsp;tol,&nbsp;maxiter,&nbsp;M1,&nbsp;M2,&nbsp;...])</p></td>
<td><p>Use Quasi-Minimal Residual iteration to solve <code class="docutils literal notranslate"><span class="pre">Ax</span> <span class="pre">=</span> <span class="pre">b</span></code>.</p></td>
</tr>
<tr class="row-odd"><td><p><a class="reference internal" href="generated/scipy.sparse.linalg.gcrotmk.html#scipy.sparse.linalg.gcrotmk" title="scipy.sparse.linalg.gcrotmk"><code class="xref py py-obj docutils literal notranslate"><span class="pre">gcrotmk</span></code></a>(A,&nbsp;b[,&nbsp;x0,&nbsp;tol,&nbsp;maxiter,&nbsp;M,&nbsp;...])</p></td>
<td><p>Solve a matrix equation using flexible GCROT(m,k) algorithm.</p></td>
</tr>
<tr class="row-even"><td><p><a class="reference internal" href="generated/scipy.sparse.linalg.tfqmr.html#scipy.sparse.linalg.tfqmr" title="scipy.sparse.linalg.tfqmr"><code class="xref py py-obj docutils literal notranslate"><span class="pre">tfqmr</span></code></a>(A,&nbsp;b[,&nbsp;x0,&nbsp;tol,&nbsp;maxiter,&nbsp;M,&nbsp;callback,&nbsp;...])</p></td>
<td><p>Use Transpose-Free Quasi-Minimal Residual iteration to solve <code class="docutils literal notranslate"><span class="pre">Ax</span> <span class="pre">=</span> <span class="pre">b</span></code>.</p></td>


### Nonlinear Systems

Nonlinear systems are solved using optimization algorithms from [scipy optimize](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html). As mentioned before, in nonlinear systems the global matrix depends on the values of the vector $\mathbf{x}$.

The steps to solve this problem are summarized below:

- 1) Set an initial value of $k^i$.
- 2) Compute an initial solution to the problem.
- 3) Use the initial values of $\mathbf{x}$ to compute a new estimate of $\mathbf{A}\left(\mathbf{x}\right)$. For this, the user defines a function where the input data are the vector of unknowns $\mathbf{x}$ and the respective object of the NetSimulation class.
- 4) Solve.
    - If convergence is reached, the program terminates.
    - If convergence is not reached, return to step 3.

The available optimizers are shown below:


- ‘Nelder-Mead’
- ‘Powell’
- ‘CG’ 
- ‘BFGS’ 
- ‘Newton-CG’ 
- ‘L-BFGS-B’ 
- ‘TNC’ 
- ‘COBYLA’
- ‘SLSQP’ 
- ‘trust-constr’
- ‘dogleg’
- ‘trust-ncg’ 
- ‘trust-exact’ 
- ‘trust-krylov’ 

## Requirements

The code works for me with the following:

- python 3.7.0
- numpy 1.21.5
- scipy 1.4.1
- numba 0.55.1
- matplotlib 3.5.3
