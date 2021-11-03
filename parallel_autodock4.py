#!/usr/bin/env python
# -*- coding: utf-8 -*-

#    This script is of the Autodock deep-learning toxicology GUI developed by McCarthy
#
#    Copyright (C) 2021  Michael J McCarthy
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>

"""

Run with (making sure MPI and mpi4py are installed):

$ mpirun -n X python parallel_autodock4.py

where X is the number of processes you want to run this on.

This is a modified script from the original author Lion Krischer who had made it freely available

"""
from mpi4py import MPI
import os
import glob
from sys import platform

COMM = MPI.COMM_WORLD
def split(container, count):
    return [container[_i::count] for _i in range(count)]


if COMM.rank == 0:
   if os.path.isfile("list-left-to-dock.txt"):
        jobs = open("list-left-to-dock.txt").readlines()
        jobs = split(jobs, COMM.size)
   else:
        jobs = open("list-to-dock.txt").readlines()
        jobs = split(jobs, COMM.size)
else:
   jobs = None

# Scatter jobs across cores.
jobs = COMM.scatter(jobs, root=0)

for job in jobs:
    # determine the OS this is running on and execute
    if platform == 'darwin':
        os.system("autodock4 -p " + job )
    elif platform == "win32":
        os.system(".\\autodock4.exe -p " + job )
    elif platform == "linux" or "linux2":
        os.system("autodock4 -p " + job )
