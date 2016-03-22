import numpy
from general import solvers, wavefunctions


class Potential(object):
    def __init__(self, wf):
        self._v = numpy.zeros([wf.BASE_SIZE])
        for i in range(wf.BASE_SIZE):
            self._v[i] = -1.*i*(wf.BASE_SIZE-i-1)/(wf.BASE_SIZE*wf.BASE_SIZE) - 5.

    def operate(self, wf):
        return self._v * wf._psi

class Solver(solvers.Solver1D):
    def __init__(self):
        wf = wavefunctions.GridWavefunction1D()
        potential = Potential(wf)
        super(Solver, self).__init__(wf, potential)
