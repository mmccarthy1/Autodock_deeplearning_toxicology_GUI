

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

from mpi4py import MPI
import fileinput
import re
import os
import sys
from os import path
import subprocess
from sys import platform
import json
import ast
import errno

COMM = MPI.COMM_WORLD


def split(container, count):
    return [container[_i::count] for _i in range(count)]
if COMM.rank == 0:
    if os.path.isfile("list_for_binana.json"):
        with open("list_for_binana.json", 'rb') as jsonfile:
            data = jsonfile.read()
        jobs = json.loads(data)
        jobs = split(jobs, COMM.size)

else:
    jobs = None
#jobs = COMM.bcast(jobs, root=0)
jobs = COMM.scatter(jobs, root=0)

for job in jobs:
    for first_part, second_part in job.items():
        if platform == 'darwin':
            print("found mac")
            firsta, firstb, firstc, firstd, firste = first_part.split()
            subprocess.run(["python",firsta,firstb,firstc,firstd,firste], stdout=open(second_part,"w"))
        elif platform == "win32":
            print("found windows")
            firsta, firstb, firstc, firstd, firste = first_part.split()
            subprocess.run(["python",firsta,firstb,firstc,firstd,firste], stdout=open(second_part,"w"))
        elif platform == "linux" or "linux2":
            print("found linux")
            firsta, firstb, firstc, firstd, firste = first_part.split()
            subprocess.run(["python",firsta,firstb,firstc,firstd,firste], stdout=open(second_part,"w"))
            #subprocess.run(["python",first_part],stdout=open(second_part,"w"))
