import numpy
import matplotlib.pyplot as plt
from general import wavefunctions

class Solver(object):
    def __init__(self, wf, potential):
        self.wavefunction = wf
        self.potential = potential
        self.last_energy = 0
        self.wf_list = []

    def lanczos(self):
        # Loop long enough, 
        for i in range(40000):
            # Just apply Hessian operator for a given potential
            new_wf = self.wavefunction.hessian(self.potential)
            
            # Orthogonalise with respect to previously found solutions
            for old_wf in self.wf_list:
                new_wf.orthogonalize(old_wf)
            
            # Figure out the proportionality constant (energy)
            self.last_energy = self.wavefunction.projection(new_wf)
            
            # Normalize the resulting wavefunction
            new_wf.normalize()
            
            # Check if we already converged
            if self.wavefunction.difference(new_wf) < 0.1 and i > 100:
                break
            
            # Store the new wavefunction and loop over
            self.wavefunction = new_wf

        print("E:", self.last_energy)
        self.wf_list.append(self.wavefunction)


class Solver1D(Solver):
    def plot(self):
        plt.subplot(211)
        plt.title('Particle density')
        for wf in self.wf_list:
            wf.plot_density(plt)
    
        plt.subplot(212)
        plt.title('Particle wavefunction')
        for wf in self.wf_list:
            wf.plot_psi(plt)
        
        plt.show()

class Solver2D(Solver):
    def plot(self):
        numpl = len(self.wf_list)
        width = int(numpy.sqrt(numpl))
        height = numpl / width

        if width*height < numpl:
            height += 1
            
        for i, wf in enumerate(self.wf_list):
            subplot_string = '{}{}{}'.format(width, height, i+1)
            plt.subplot(subplot_string)
            wf.plot_psi(plt)
        
        plt.show()