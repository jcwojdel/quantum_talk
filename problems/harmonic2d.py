from general import solvers, wavefunctions
import numpy

POTENTIAL = -10.0

class Potential(object):
    def __init__(self, wf):
        self._v = numpy.zeros([wf.BASE_SIZE, wf.BASE_SIZE])
        for i in range(wf.BASE_SIZE):
            for j in range(wf.BASE_SIZE):
                self._v[i, j] = -1.*i*(wf.BASE_SIZE-i-1)/(wf.BASE_SIZE*wf.BASE_SIZE)-1.*j*(wf.BASE_SIZE-j-1)/(wf.BASE_SIZE*wf.BASE_SIZE) - 5.

    def operate(self, wf):
        return self._v * wf._psi

class Solver(solvers.Solver2D):
    def __init__(self):
        wf = wavefunctions.GridWavefunction2D()
        potential = Potential(wf)
        super(Solver, self).__init__(wf, potential)
