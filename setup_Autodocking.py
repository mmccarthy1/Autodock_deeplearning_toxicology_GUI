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

import fileinput
import re
import os
import sys
import rdkit
import time
import datetime
import shutil
import glob
from sys import platform
import subprocess
from openbabel import openbabel
from openbabel import pybel
from openbabel import OBChargeModel


class setup_docking():

    def __init__(self, *args, **kwargs):
            #super.__init__(self)
            #self.working_project = working_project
            self.args = args
            self.kwargs = kwargs

    def molecule_title(self, file_name, mol_name):
        with open(file_name, 'r+') as curr_file:
            old_content = curr_file.read()
            curr_file.seek(0,0)
            curr_file.write("REMARK  Name = " + mol_name + "\n" + old_content)
            #print(curr_file)
            curr_file.close()
        return

    def convert_pdb(self, split_sdf_file):
        #ob_charge_model = openbabel.OBChargeModel.FindType("gasteiger")
        for find_name in pybel.readfile("sdf", split_sdf_file):
            name = find_name.title
            print("the name that was found: ", name)
        if platform == "darwin":
            obConversion = openbabel.OBConversion()
            obConversion.SetInAndOutFormats("sdf", "pdbqt")
            mol_dock = openbabel.OBMol()
            ob_charge_model = openbabel.OBChargeModel.FindType("gasteiger")
            obConversion.ReadFile(mol_dock, str(split_sdf_file))
            mol_dock.AddHydrogens()
            mol_dock.CorrectForPH()
            ob_charge_model.ComputeCharges(mol_dock)
            ob_charge_model.GetPartialCharges()
            obConversion.WriteFile(mol_dock, name+'.pdbqt')
            pdbqt_file_name = name+'.pdbqt'
            if os.path.isfile(pdbqt_file_name):
                self.molecule_title(pdbqt_file_name, mol_name)
            else:
                pass
        else:
            print("not MAC")
            obConversion = openbabel.OBConversion()
            obConversion.SetInAndOutFormats("sdf", "pdb")
            mol_dock = openbabel.OBMol()
            obConversion.ReadFile(mol_dock, str(split_sdf_file))
            mol_dock.AddHydrogens()
            obConversion.WriteFile(mol_dock, name+'.pdb')
        return
    def convert_pdbqt(self, pdb_file, pythonsh_loc, prepare_ligand):
        with open(pdb_file, "r") as ifile:
            for line in ifile:
                if line.startswith("COMPND"):
                   mol_name = re.search(r'COMPND    (\w+\-+.+)$', line)
                   if mol_name == None:
                       mol_name = pdb_file.replace(".pdb","")
                       subprocess.run([pythonsh_loc,prepare_ligand,"-l",pdb_file, "-A","Hydrogens"])
                       pdbqt_file_name = pdb_file.replace(".pdb",".pdbqt")
                       if os.path.isfile(pdbqt_file_name):

                          self.molecule_title(pdbqt_file_name, mol_name)
                       else:
                          continue
                   else:
                       subprocess.run([pythonsh_loc,prepare_ligand,"-l",pdb_file, "-A","Hydrogens"])
                       pdbqt_file_name = pdb_file.replace(".pdb",".pdbqt")
                       if os.path.isfile(pdbqt_file_name):
                          self.molecule_title(pdbqt_file_name, mol_name.group(1))
                       else:
                          continue
    #ifile.close()

    def setup_mol_file_for_docking(self, molecule_files, install_dir):
    ## The for-loop breaks up the molecule file
       os.chdir(molecule_files)
       for sdf_file in os.listdir(molecule_files):
           if sdf_file.endswith(".sdf"):
              for mol in pybel.readfile("sdf", sdf_file):
                  name = mol.title
                  if name == '':
                      continue
                  else:
                      mol.write("sdf", "%s.sdf" % name)
                  ##the following was the old method to do above and can probably be removed
                  #os.system("obabel " + sdf_file + " -O " + sdf_file + " -m ")
              # change the next two line so they don't use system call use mkdir and move
           os.mkdir("ligand_files_"+datetime.datetime.now().strftime("%m-%d-%y"))
           shutil.move(sdf_file, "./ligand_files_"+datetime.datetime.now().strftime("%m-%d-%y"))
       return
    def generate_pdb(self, molecule_files, install_dir):
       for split_sdf_file in os.listdir(molecule_files):
            if split_sdf_file.endswith(".sdf"):
               self.convert_pdb(split_sdf_file)
       #os.chdir(molecule_files)
       return
    def generate_pdbqt(self, molecule_files, install_dir):
        if platform == "win32":
            print("windows was found")
            prepare_ligand = os.path.join(install_dir,"mgltools_win32_1.5.6\\AutoDockTools\\Utilities24\\prepare_ligand4.py")
            pythonsh_loc = os.path.join(install_dir,"python.exe")
        elif platform == "linux" or "linux2":
            print("linux was found")
            prepare_ligand = os.path.join(install_dir,"mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py")
            pythonsh_loc = os.path.join(install_dir,"mgltools_x86_64Linux2_1.5.6/bin/pythonsh")
        for pdb_file in os.listdir(molecule_files):
            if pdb_file.endswith(".pdb"):
               self.convert_pdbqt(pdb_file, pythonsh_loc, prepare_ligand)
            else:
               pass
        return
    def setup_dpf_files(self, molecule_files, protein_files):
        print("running setup_dpf_files")
        for file in os.listdir(molecule_files):
            filename = os.fsdecode(file)
            if filename.endswith(".pdbqt"):
               full_name = os.path.join(molecule_files, filename)
               #print(filename)
               mol_info = next(pybel.readfile("pdbqt", filename))
               rot_bonds = mol_info.OBMol.NumRotors()
               ##This will adjust the number of energy evaluations based on the number of rotatable bonds
               if 0 <= rot_bonds <= 1:
                  set_num = 1000000
               if 2 <= rot_bonds < 5:
                  set_num = 2500000
               if 5 <= rot_bonds <= 10:
                  set_num = 4500000
               if 11 <= rot_bonds:
                  set_num = 10000000

               for file_pdbqt, subdirs, files in os.walk(protein_files):
                  for name in files:
                     output = os.path.join(file_pdbqt, name)
                     pdbqt_file = os.fsdecode(output)
                     if pdbqt_file.endswith(".pdbqt"):
                           protein_file = os.path.splitext(pdbqt_file)[0]
                           new_file = open(os.path.splitext(pdbqt_file)[0]+"_"+os.path.splitext(filename)[0]+".dpf", "w")
                           new_file.write("autodock_parameter_version 4.2       # used by autodock to validate parameter set \n"
                           "outlev 1                             # diagnostic output level\n"
                           "intelec                              # calculate internal electrostatics\n"
                           "seed pid time                        # seeds for random generator\n"
                           "ligand_types A Br C Cl HD F I N NA OA P SA S Zn  # atoms types in ligand\n"
                           "fld  "+protein_file+".maps.fld         # grid_data_file\n"
                           "map  "+protein_file+".A.map            # atom-specific affinity map\n"
                           "map  "+protein_file+".Br.map           # atom-specific affinity map\n"
                           "map  "+protein_file+".C.map            # atom-specific affinity map\n"
                           "map  "+protein_file+".Cl.map           # atom-specific affinity map\n"
                           "map  "+protein_file+".HD.map           # atom-specific affinity map\n"
                           "map  "+protein_file+".F.map            # atom-specific affinity map\n"
                           "map  "+protein_file+".I.map            # atom-specific affinity map\n"
                           "map  "+protein_file+".N.map            # atom-specific affinity map\n"
                           "map  "+protein_file+".NA.map           # atom-specific affinity map\n"
                           "map  "+protein_file+".OA.map           # atom-specific affinity map\n"
                           "map  "+protein_file+".P.map            # atom-specific affinity map\n"
                           "map  "+protein_file+".SA.map           # atom-specific affinity map\n"
                           "map  "+protein_file+".S.map            # atom-specific affinity map\n"
                           "map  "+protein_file+".Zn.map            # atom-specific affinity map\n"
                           "elecmap "+protein_file+".e.map         # electrostatics map\n"
                           "desolvmap "+protein_file+".d.map       # desolvation map\n"
                           "move "+full_name+"  # small molecule\n"
                           "tran0 random                         # initial coordinates/A or random\n"
                           "quaternion0 random                   # initial orientation\n"
                           "dihe0 random                         # initial dihedrals (relative) or random\n"
                           "torsdof 5                            # torsional degrees of freedom\n"
                           "rmstol 2.0                           # cluster_tolerance/A\n"
                           "extnrg 1000.0                        # external grid energy\n"
                           "e0max 0.0 10000                      # max initial energy; max number of retries\n"
                           "ga_pop_size 150                      # number of individuals in population\n"
                           "ga_num_evals "+str(set_num)+"        # maximum number of energy evaluations\n"
                           "ga_num_generations 27000             # maximum number of generations\n"
                           "ga_elitism 1                         # number of top individuals to survive to next generation\n"
                           "ga_mutation_rate 0.02                # rate of gene mutation\n"
                           "ga_crossover_rate 0.8                # rate of crossover\n"
                           "ga_window_size 10                    # \n"
                           "ga_cauchy_alpha 0.0                  # Alpha parameter of Cauchy distribution\n"
                           "ga_cauchy_beta 1.0                   # Beta parameter Cauchy distribution\n"
                           "set_ga                               # set the above parameters for GA or LGA\n"
                           "sw_max_its 300                       # iterations of Solis & Wets local search\n"
                           "sw_max_succ 4                        # consecutive successes before changing rho\n"
                           "sw_max_fail 4                        # consecutive failures before changing rho\n"
                           "sw_rho 1.0                           # size of local search space to sample\n"
                           "sw_lb_rho 0.01                       # lower bound on rho\n"
                           "ls_search_freq 0.06                  # probability of performing local search on individual\n"
                           "set_psw1                             # set the above pseudo-Solis & Wets parameters\n"
                           "unbound_model bound                  # state of unbound ligand\n"
                           "ga_run 10                            # do this many hybrid GA-LS runs\n"
                           "analysis                             # perform a ranked cluster analysis\n ")
                           new_file.close()
        return

    def move_files_for_docking(self, protein_files, molecule_files, dpf_files):
        for dpf_file, subdirs, files in os.walk(protein_files):
            for filename_dpf in files:
                if filename_dpf.endswith(".dpf"):
                    #print(filename_dpf)
                    shutil.move(dpf_file+"/"+filename_dpf, dpf_files)

        for mol_file, subdirs, files in os.walk(molecule_files):
             for filename_mol in files:
                #print("This is the filename_mol:", filename_mol)
                if filename_mol.endswith(".pdbqt"):
                    #print(filename_mol)
                    shutil.copy(mol_file+"/"+filename_mol, dpf_files)

        return

    def make_list_to_dock(self, dpf_files, working_project):
        os.chdir(dpf_files)
        with open("list-to-dock.txt", "w") as new_file:
            for work_file, subdirs, files in os.walk(dpf_files):
                for filename_work in files:
                   if filename_work.endswith(".dpf"):
                       #print(filename_work)
                       dpfile = os.path.join(dpf_files, filename_work)
                       if platform == "win32":
                            new_file.write(str(dpfile)+'\n')
                       else:
                            new_file.write(str(dpfile) + os.linesep)
        new_file.close()
        listofdocked = open("list_of_docked.txt", "w")
        listofdocked.write("")
        listofdocked.close()
        for list in os.listdir(dpf_files):
            if list.startswith('list_of_docked'):
                #print("location of working_project: ", working_project)
                shutil.move(list, working_project)
            if list.startswith('list-to-dock'):
                shutil.move(list, working_project)
        return
