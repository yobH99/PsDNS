"""Single-mode Rayleigh-Taylor
Navier-Stokes equations with Boussinesq buoyancy terms set up for a
single mode Rayleigh-Taylor simulation.  The initial disturbance is a
velocity disturbance that roughly corresponds to the linear
instability eigenfunction.
"""
import numpy as np
from psdns import *
from psdns.equations.navier_stokes import Boussinesq

#define parameter for the output 
tdump_data = 0.1
tdump = 0.01
Lx = 4.0 *numpy.pi
Ly = 4.0 *numpy.pi
Lz = 4.0 *numpy.pi

data_folder = '/home/yobh/Desktop/LANL/PsDNS/caso_buoyant/data/'
spec_folder = '/home/yobh/Desktop/LANL/PsDNS/caso_buoyant/spectra_2D/'
grid = SpectralGrid(
    sdims=[127, 127, 127],
    pdims=[191, 191, 191],
    box_size=[Lx, Ly, Lz]
    )
grid.checkpoint(data_folder + "data.grid")
equations = Boussinesq(Re=100, g=-1, At=0.1, L=grid.box_size[0])

x = grid.x[:2,:,:,0]
solver = RungeKutta(
    dt=0.01,
    tfinal=10,
    equations=equations,
    ic=equations.perturbed_interface(
        grid,
        equations.band(grid, 1,12,1e-4),
        1,
        0.5,
        ),
    diagnostics=[
        FieldDump(tdump=tdump_data, grid=grid, filename= data_folder + "data{:04g}"),
        StandardDiagnostics(
            tdump=tdump, grid=grid,
            fields=['divU','tke','rho2','cavg'],
            outfile= data_folder + "std.dat"
            ),
            spectral_length_2D(tdump=tdump, grid=grid, outfile=spec_folder + 'spectra2'),
            spectral_length_3D(tdump=100, grid=grid, outfile=data_folder+'spectra.dat')
        ],
    )
solver.run()
solver.print_statistics()
