from general import solvers, wavefunctions

POTENTIAL = -5.0

class Potential(object):
    def __init__(self):
        self._v = POTENTIAL

    def operate(self, wf):
        return self._v * wf._psi

class Solver(solvers.Solver):
    def __init__(self):
        wf = wavefunctions.GridWavefunction1D()
        potential = Potential()
        super(Solver, self).__init__(wf, potential)
