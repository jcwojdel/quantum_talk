#!/usr/bin/env python3
from argparse import ArgumentParser
from problems import well, harmonic

def well_problem():
    solver = well.Solver()
    return solver

def harmonic_problem():
    solver = harmonic.Solver()
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

    for _ in range(args.number):
        solver.lanczos()

    solver.plot()

if __name__ == "__main__":
    main()