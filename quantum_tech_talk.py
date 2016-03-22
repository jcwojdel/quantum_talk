#!/usr/bin/env python3
from argparse import ArgumentParser
from problems import well, harmonic, well2d, harmonic2d

def well_problem():
    solver = well.Solver()
    return solver

def well2d_problem():
    solver = well2d.Solver()
    return solver

def harmonic_problem():
    solver = harmonic.Solver()
    return solver
    
def harmonic2d_problem():
    solver = harmonic2d.Solver()
    return solver

def main():
    parser = ArgumentParser(description='Quantum TechTalk thingy')
    parser.add_argument('-p', '--problem', default='well')
    parser.add_argument('-n', '--number', type=int, default=5)
    args = parser.parse_args()

    if args.problem == 'well':
        solver = well_problem()
    elif args.problem == 'harmonic':
        solver = harmonic_problem()
    elif args.problem == 'well2d':
        solver = well2d_problem()
    elif args.problem == 'harmonic2d':
        solver = harmonic2d_problem()

    for _ in range(args.number):
        solver.lanczos()

    solver.plot()

if __name__ == "__main__":
    main()