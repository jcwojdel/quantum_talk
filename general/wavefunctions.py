import numpy


class GridWavefunction(object):
    def density(self, r):
        return self.psi(r) * numpy.conjugate(self.psi(r))

    def normalize(self):
        norm = numpy.linalg.norm(self._psi)
        self._psi = self._psi / norm
        return norm

    def hessian(self, potential):
        laplacian = self.laplacian()

        new_wf = self.__class__(-laplacian + potential.operate(self))
        return new_wf
    
    def set_psi(self, psi):
        self._psi[:] = psi[:]
        self.normalize()
    
    def projection(self, other_wf):
        return numpy.sum(self._psi*other_wf._psi)

    def difference(self, other_wf):
        return numpy.linalg.norm(self._psi - other_wf._psi)

    def orthogonalize(self, other_wf):
        self._psi -= self.projection(other_wf) * other_wf._psi


class GridWavefunction1D(GridWavefunction):
    BASE_SIZE = 301
    MIN_X = -1.
    MAX_X = 1.

    def __init__(self, psi=None):
        self._psi = numpy.zeros([self.BASE_SIZE])
        if psi is not None:
            self._psi[:] = psi[:]
        else:
            self.set_initial()

    def set_initial(self):
        for i in range(self.BASE_SIZE):
            self._psi[i] = (i *  (self.BASE_SIZE - i))
        self.normalize()

    def psi(self, x):
        if x < self.MIN_X or x > self.MAX_X:
            return 0.0
        relative = (self.BASE_SIZE - 1)*(x - self.MIN_X)/(self.MAX_X - self.MIN_X)
        offset = numpy.floor(relative)
        partial = relative - offset
        result = self._psi[offset]*(1-partial)
        if partial > 1e-9:
            result += self._psi[offset+1]*partial
        return result 
        
    def normalize(self):
        self._psi[0] = 0.0
        self._psi[-1] = 0.0
        return GridWavefunction.normalize(self)

    def laplacian(self):
        result = numpy.zeros([self.BASE_SIZE])
        result[0] = 0. + self._psi[1] - 2 * self._psi[0]
        result[-1] = self._psi[-2] + 0.0 - 2 * self._psi[-1]
        result[1:-1] = self._psi[0:-2] + self._psi[2:] - 2 * self._psi[1:-1]
        return result
    
    def plot_psi(self, plt):
        xvals = numpy.linspace(self.MIN_X, self.MAX_X, 200)
        plt.plot(xvals, [self.psi(x) for x in xvals])

    def plot_density(self, plt):
        xvals = numpy.linspace(self.MIN_X, self.MAX_X, 200)
        plt.plot(xvals, [self.density(x) for x in xvals])


class GridWavefunction2D(GridWavefunction):
    BASE_SIZE = 51
    MIN_X = -1.
    MAX_X = 1.
    MIN_Y = -1.
    MAX_Y = 1.

    def __init__(self, psi=None):
        self._psi = numpy.zeros([self.BASE_SIZE, self.BASE_SIZE])
        if psi is not None:
            self._psi[:,:] = psi[:,:]
        else:
            self.set_initial()

    def set_initial(self):
        for x in range(self.BASE_SIZE):
            for y in range(self.BASE_SIZE):
                self._psi[x, y] = (x *  (self.BASE_SIZE - x))*(y *  (self.BASE_SIZE - y))
        self.normalize()

    def psi(self, r):
        x, y = r
        if x < self.MIN_X or x > self.MAX_X:
            return 0.0
        if y < self.MIN_Y or y > self.MAX_Y:
            return 0.0

        relativex = (self.BASE_SIZE - 1)*(x - self.MIN_X)/(self.MAX_X - self.MIN_X)
        relativey = (self.BASE_SIZE - 1)*(y - self.MIN_Y)/(self.MAX_Y - self.MIN_Y)
        offsetx = numpy.floor(relativex)
        offsety = numpy.floor(relativey)
        partialx = relativex - offsetx
        partialx = relativex - offsetx
        result = self._psi[offsetx, offsety]
        #*(1-partial)
        #if partial > 1e-9:
        #    result += self._psi[offset+1]*partial
        return result 
    
    def normalize(self):
        self._psi[0, :] = 0.0
        self._psi[-1, :] = 0.0
        self._psi[:, 0] = 0.0
        self._psi[:, -1] = 0.0
        return GridWavefunction.normalize(self)

    def laplacian(self):
        result = numpy.zeros([self.BASE_SIZE, self.BASE_SIZE])
        result[1:, :] += self._psi[0:-1, :]
        result[:-1,:] += self._psi[1:,   :]
        result[:, 1:] += self._psi[:, 0:-1]
        result[:,:-1] += self._psi[:,   1:]
        result[1:-1, 1:-1] -= 4 * self._psi[1:-1, 1:-1]
        return result

    def plot_psi(self, plt):
        vmax = numpy.max(numpy.abs(self._psi[:]))

        plt.imshow(self._psi, vmax=vmax, vmin=-vmax)

    def plot_density(self, plt):
        plt.imshow(self._psi * self._psi)
