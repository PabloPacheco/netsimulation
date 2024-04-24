import numpy as np
import scipy as sp


def nonlinear_problem(xF, mynet):
    dx = np.zeros(mynet.n_elements)
    inddt = np.zeros(mynet.n_nodes, dtype = int)
    inddt[mynet.nodes_x_known[:]] = 1
    x = np.zeros(mynet.n_nodes) 
    x[mynet.nodes_x_known[:]] = mynet.values_x_known[:]
    x[inddt == 0] = xF[:]
    xE = mynet.values_x_known

    for i in range(mynet.n_elements):
        c0 = mynet.connectivity[i,0]
        c1 = mynet.connectivity[i,1]
        dx[i] = x[c0] - x[c1]    
        
    mynet.k_element = mynet.k_function(x, mynet) 
    mynet.matrix_assembly()  
    vn_nodes = np.linspace(0, mynet.n_nodes - 1, mynet.n_nodes, dtype = int) 

    inddt[mynet.nodes_x_known[:]] = 1  
    indAFx, indAFy = np.meshgrid(vn_nodes[inddt == 0], vn_nodes[inddt == 0])
    indAEFx, indAEFy = np.meshgrid(vn_nodes[inddt == 1], vn_nodes[inddt == 0])
    indAEx, indAEy = np.meshgrid(vn_nodes[inddt == 1], vn_nodes[inddt == 1])
    AF = mynet.A[indAFx, indAFy]
    bF = mynet.values_b_known
    AEF = mynet.A[indAEFx, indAEFy]

    res = AEF.dot(xE) + AF.dot(xF) - bF
    return np.linalg.norm(res)
    
    
def k_calc(x, mynet):
    dx = np.zeros(mynet.n_elements)
    for i in range(mynet.n_elements):
        c0 = mynet.connectivity[i,0]
        c1 = mynet.connectivity[i,1]
        dx[i] = x[c0] - x[c1]   
    chw = 110
    diam = np.array([0.1524, 
                     0.1524, 
                     0.1524, 
                     0.1270, 
                     0.1016, 
                     0.1270, 
                     0.1016, 
                     0.1524, 
                     0.1270, 
                     0.1270, 
                     0.1016, 
                     0.1524])   
    k = 0.2784 * chw  * diam[:]**2.63 / (mynet.length[:]**0.54) / (np.abs(dx)**0.46)
    return k
        

class NetSimulation:
    def __init__(self):
        #input data 
        self.n_elements = 5
        self.n_nodes = 4
        self.connectivity = np.array([[0, 1], 
                                     [0, 2],
                                     [1, 2],
                                     [1, 3],
                                     [2, 3]])
        self.cxy = np.zeros((self.n_elements, 2))
        
        self.nodes_x_known = np.array([0, 3])
        self.values_x_known = np.array([20, 10])
        self.nodes_b_known = np.array([1, 2])
        self.values_b_known = np.array([10, 10])
        self.k_element = np.ones(self.n_elements)
        self.length = np.array([1000, 1000, 2000, 2000, 2000])
        self.solver = "cg"
        self.solver_tol = 1e-6
        #spsolve, bicg, bicgstab, cg, cgs, gmres, lgmres, minres, qmr, gcrotmk, tfqmr
        self.solver_nonlinear = "CG"
        self.solver_tol_nonlinear = 1e-5
        # Nelder-Mead, Powell, CG, BFGS, Newton-CG, L-BFGS-B, TNC, 
        # COBYLA, SLSQP, trust-constr, dogleg, trust-ncg, trust-exact, trust-krylov
        self.k_function = k_calc
        self.nonlinear_optimizer_function = None
        
        
    def calculate_length_coordinates(self):
        self.length = np.ones(self.n_elements)
        for i in range(self.n_elements):
            c1 = self.connectivity[i, 0]
            c2 = self.connectivity[i, 1]
            p1 = self.cxy[c1]
            p2 = self.cxy[c2]
            diff = p2 - p1
            self.length[i] = np.linalg.norm(diff) 
            
    def matrix_assembly(self):
        Ae1d = np.zeros((self.n_elements, 2, 2))
        Le1d = np.zeros((self.n_elements, 2, self.n_nodes))
        A = sp.sparse.lil_matrix((self.n_nodes, self.n_nodes))
        
        for i in range(self.n_elements):
            Ae1d[i, 0, 0] = self.k_element[i]
            Ae1d[i, 0, 1] = -self.k_element[i]
            Ae1d[i, 1, 0] = -self.k_element[i]
            Ae1d[i, 1, 1] = self.k_element[i]
            
        for i in range(self.n_elements):
            c0 = self.connectivity[i,0]
            c1 = self.connectivity[i,1]
            Le1d[i, 0, c0] = 1.
            Le1d[i, 1, c1] = 1.
            
        for i in range(self.n_elements):
            LL = Le1d[i]
            kk = Ae1d[i]
            m1 = np.dot(np.transpose(LL), kk)
            KK = np.dot(m1, LL)
            A[:,:] += KK
        self.A = A
        
        
    def solve(self):
        vn_nodes = np.linspace(0, self.n_nodes - 1, self.n_nodes, dtype = int)
        b = np.zeros(self.n_nodes)
        x = np.zeros(self.n_nodes)        
        Q = np.zeros(self.n_elements)        
        dx = np.zeros(self.n_elements)       
        inddt = np.zeros(self.n_nodes, dtype = int)
        inddt[self.nodes_x_known[:]] = 1
        xF0 = np.ones_like(np.where(inddt == 0)[0])
        bE0 = np.ones_like(np.where(inddt == 1)[0])
        for i in range(len(self.nodes_b_known)):
            ind = self.nodes_b_known[i]
            b[ind] = self.values_b_known[i]
        xE = self.values_x_known
        indAFx, indAFy = np.meshgrid(vn_nodes[inddt == 0], vn_nodes[inddt == 0])
        indAEFx, indAEFy = np.meshgrid(vn_nodes[inddt == 1], vn_nodes[inddt == 0])
        indAEx, indAEy = np.meshgrid(vn_nodes[inddt == 1], vn_nodes[inddt == 1])
        ##########################
        # partition method
        ##########################
        AF = self.A[indAFx, indAFy]
        bF = b[inddt == 0]
        AEF = self.A[indAEFx, indAEFy]
        AE = self.A[indAEx, indAEy]
        if self.solver == "cg":
            xF, solver_info = sp.sparse.linalg.cg(AF.tocsr(), bF.reshape((len(bF), 1)) - AEF.dot(xE.reshape((len(xE), 1))), atol = self.solver_tol)
        elif self.solver == "bicg":
            xF, solver_info = sp.sparse.linalg.bicg(AF.tocsr(), bF.reshape((len(bF), 1)) - AEF.dot(xE.reshape((len(xE), 1))), atol = self.solver_tol)
        elif self.solver == "bicgstab":
            xF, solver_info = sp.sparse.linalg.bicgstab(AF.tocsr(), bF.reshape((len(bF), 1)) - AEF.dot(xE.reshape((len(xE), 1))), atol = self.solver_tol)
        elif self.solver == "cgs":
            xF, solver_info = sp.sparse.linalg.cgs(AF.tocsr(), bF.reshape((len(bF), 1)) - AEF.dot(xE.reshape((len(xE), 1))), atol = self.solver_tol)           
        elif self.solver == "gmres":
            xF, solver_info = sp.sparse.linalg.gmres(AF.tocsr(), bF.reshape((len(bF), 1)) - AEF.dot(xE.reshape((len(xE), 1))), atol = self.solver_tol) 
        elif self.solver == "minres":
            xF, solver_info = sp.sparse.linalg.minres(AF.tocsr(), bF.reshape((len(bF), 1)) - AEF.dot(xE.reshape((len(xE), 1))), atol = self.solver_tol) 
        elif self.solver == "lgmres":
            xF, solver_info = sp.sparse.linalg.lgmres(AF.tocsr(), bF.reshape((len(bF), 1)) - AEF.dot(xE.reshape((len(xE), 1))), atol = self.solver_tol) 
        else:
            xF, solver_info = sp.sparse.linalg.cg(AF.tocsr(), bF.reshape((len(bF), 1)) - AEF.dot(xE.reshape((len(xE), 1))), atol = self.solver_tol)
            
        xF = xF.reshape(len(xF), 1)

        bE = AE.dot(xE.reshape((len(xE), 1))) + AEF.transpose().dot(xF.reshape((len(xF), 1)))

        x[self.nodes_x_known[:]] = self.values_x_known[:]
        x[inddt == 0] = xF[:,0]
        b[inddt == 1] = bE[:,0]
        
        for i in range(self.n_elements):
            c0 = self.connectivity[i,0]
            c1 = self.connectivity[i,1]
            dx[i] = x[c0] - x[c1]
            Q[i] = self.k_element[i] * dx[i]   
            
        self.b = b
        self.x = x
        self.dx = dx
        self.Q = Q        
        
    def solve_nonlinear(self):
        # calculations for initial estimation of xF
        self.matrix_assembly()
        self.solve()
        inddt = np.zeros(self.n_nodes, dtype = int)
        inddt[self.nodes_x_known[:]] = 1
        vn_nodes = np.linspace(0, self.n_nodes - 1, self.n_nodes, dtype = int)
        xF = self.x[inddt == 0]
        # calculation of xF
        if self.nonlinear_optimizer_function == None:
            result = sp.optimize.minimize(nonlinear_problem, 
                                          xF, 
                                          args = (self), 
                                          options = {'disp': False, 'maxiter': 1000}, 
                                          tol = self.solver_tol_nonlinear,
                                          method = self.solver_nonlinear)
            xF = result.x
        else:
            result = self.nonlinear_optimizer_function(self, xF)
            xF = result.x
        
        self.x[inddt == 0] = xF
        self.k_element = self.k_function(self.x, self)
        self.matrix_assembly()
        # dx, Q calculation
        for i in range(self.n_elements):
            c0 = self.connectivity[i,0]
            c1 = self.connectivity[i,1]
            self.dx[i] = self.x[c0] - self.x[c1]
            self.Q[i] = self.k_element[i] * self.dx[i]  
        # k_element
        xE = self.values_x_known
        indAEFx, indAEFy = np.meshgrid(vn_nodes[inddt == 1], vn_nodes[inddt == 0])
        indAEx, indAEy = np.meshgrid(vn_nodes[inddt == 1], vn_nodes[inddt == 1])
        AEF = self.A[indAEFx, indAEFy]
        AE = self.A[indAEx, indAEy]
        bE = AE.dot(xE.reshape((len(xE), 1))) + AEF.transpose().dot(xF.reshape((len(xF), 1)))
        self.b[inddt == 1] = bE[:,0]
        return

        

        