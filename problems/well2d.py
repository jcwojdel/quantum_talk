from general import solvers, wavefunctions

POTENTIAL = -10.0

class Potential(object):
    def __init__(self):
        self._v = POTENTIAL

    def operate(self, wf):
        return self._v * wf._psi

class Solver(solvers.Solver2D):
    def __init__(self):
        wf = wavefunctions.GridWavefunction2D()
        potential = Potential()
        super(Solver, self).__init__(wf, potential)
