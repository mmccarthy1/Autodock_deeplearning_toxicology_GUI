# -*- coding: utf-8 -*-


#    This will start and run all functions of the Autodock deep-learning toxicology GUI developed by McCarthy
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

import time
import sys
import os
import tarfile
import datetime
import shutil
import csv
import re
#import h5py
import json
import multiprocessing
from os.path import join
from openbabel import openbabel as OB
#from openbabel import pybel as PB
from sys import platform
import pubchempy as pcp
import numpy as np
import pandas as pd
import plotly as py
import plotly.graph_objs as go
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import PandasTools
from rdkit.Chem import Draw
from tensorflow.keras.preprocessing.image import ImageDataGenerator
from keras.utils import np_utils
import tensorflow as tf
from tensorflow.keras.models import load_model
import setup_Autodocking
import parse_dlg
import report_of_results_pylatex
import count_contacts_bina
import generate_3Dligands_10_28_2020
import errno
import subprocess
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtPrintSupport import *
from PyQt5 import QtWebEngineWidgets
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5 import QtWidgets#, uic

QtWidgets.QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
QtWidgets.QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)

class WorkerSignals(QObject):
    '''
    Defines the signals available from a running worker thread.

    Supported signals are:

    finished
        No data

    error
        `tuple` (exctype, value, traceback.format_exc() )

    result
        `object` data returned from processing, anything

    progress
        `int` indicating % progress

    '''
    started = pyqtSignal()
    finished = pyqtSignal()
    error = pyqtSignal(tuple)
    result = pyqtSignal(object)
    progress = pyqtSignal(int)

class check_progres(QThread):
    ##This class is just going to check for progress of docking while docking.
    prog_status = pyqtSignal(int)
    message_status = pyqtSignal(str)
    def run(self):
        curr_path = os.getcwd()
        dlg_files = os.path.join(curr_path, "dpf_files/")
        dock_to_do = len(open("list-to-dock.txt").readlines())
        time.sleep(5)
        tot_done = 0
        while tot_done <= dock_to_do:
            good = 0
            bad = 0
            ## it looks through the dlg_files directory for dlg files and reads through them for the final output from docking which is a determination on successful or unsuccessful
            for dlg_file_name in os.listdir(dlg_files):
                if dlg_file_name.endswith(".dlg"):
                    try:
                        with open(os.path.join(dlg_files, dlg_file_name)) as dlg_file:
                            info = dlg_file.read()
                            outcome_good = re.search("Successful", info)
                            outcome_bad = re.search("Unsuccessful", info)
                            if outcome_good:
                                # dlg_found = os.path.join(dlg_files, dlg_file_name)
                                good = good + 1
                            elif outcome_bad:
                                #dlg_found = os.path.join(dlg_files, dlg_file_name)
                                bad = bad + 1
                            else:
                                pass
                    except IOError as exc:
                        if exc.errno != errno.EISDIR: # Do not fail if a directory is found, just ignore it.
                            raise # Propagate other kinds of IOError.
            self.message = "Docking - Successful: "+str(good)+" Unsuccessful: "+str(bad)
            tot_done = good + bad
            prog = float(tot_done) / dock_to_do
            self.prog_status.emit(int(prog*100))
            self.message_status.emit(self.message)
            if (prog*100) == 100:
               break

class process_setup_progress(QThread):
    ##This class handles checking for progress of setting up a docking project
    prog_status = pyqtSignal(int)
    message_status = pyqtSignal(str)
    def __init__(self, **kwargs):
        QObject.__init__(self)
        self.kwargs = kwargs
    def run(self):
        if self.kwargs.get('param') == 'separate_mols':
            self.message = "step 1. spearating master molecule file into individual files"
            self.message_status.emit(self.message)
        elif self.kwargs.get('param') == 'gen_pdb':
            self.message = "step 2. generating PDB files"
            self.message_status.emit(self.message)
        elif self.kwargs.get('param') == 'gen_pdbqt':
            self.message = "step 3. generating PDBQT files with AutodockTools"
            self.message_status.emit(self.message)
        elif self.kwargs.get('param') == 'dpf_files':
            self.message = "step 4. Writing docking parameter files (dpf)"
            self.message_status.emit(self.message)
        elif self.kwargs.get('param') == 'moving_files':
            self.message = "step 5. Moving dpf docking files"
            self.message_status.emit(self.message)
        elif self.kwargs.get('param') == 'writing_docking_list':
            self.message = "step 6. Writing docking lists"
            self.message_status.emit(self.message)
        return

class process_results_progress(QThread):
    ##This class handles checking for progress of processing of results
    prog_status = pyqtSignal(int)
    message_status = pyqtSignal(str)
    def __init__(self, **kwargs):
        QObject.__init__(self)
        self.kwargs = kwargs
    def run(self):
            print("starting process results")
            parlog_path = os.path.join(working_project, "dlg_files/")
            tot_files = 0
            tot_dlg_files = 0
            tot_dpf_files = 0
            tot_full_log_file = 0
            file_dir = os.path.join(working_project, "dlg_files/")
            dpf_dir = os.path.join(working_project, "dpf_files/")
            result_dir = os.path.join(working_project, "results_files/")
            molecule_dir = os.path.join(working_project, "molecule_files/")
            for i in os.listdir(parlog_path):
                if i.endswith("parsed_log.txt"):
                   tot_files = tot_files + 1
                if i.endswith(".dlg"):
                   try:
                        with open(os.path.join(parlog_path, i)) as dlg_file:
                            info = dlg_file.read()
                            outcome_good = re.search("Successful", info)
                            #outcome_bad = re.search("Unsuccessful", info)
                            dlg_file.close()
                            if outcome_good:
                                tot_dlg_files = tot_dlg_files + 1
                   except IOError as exc:
                            if exc.errno != errno.EISDIR: # Do not fail if a directory is found, just ignore it.
                                raise # Propagate other kinds of IOError.
                                #break
                            else:
                                pass

                if i.endswith("full_log.txt"):
                   tot_full_log_file += 1
            for i in os.listdir(dpf_dir):
                if i.endswith(".dpf"):
                    tot_dpf_files = tot_dpf_files +1
            if self.kwargs.get('param') == 'rename':
                print("Starting to process results")
                self.message = "step 1. move docking logs to dlg_files directory"
                file_dlg_moved = []
                while len(file_dlg_moved) < tot_dlg_files:
                      self.message_status.emit(self.message)
                      for i in os.listdir(file_dir):
                          if i.endswith(".dlg")and i not in file_dlg_moved:
                              file_dlg_moved.append(i)
                              prog = len(file_dlg_moved) / float(tot_dlg_files)
                              self.prog_status.emit(int(prog*100))
            if self.kwargs.get('param') == 'write_le_pose':
                self.message = "step 2. writing lowest energy pose"
                files_pdbqt_done = []
                while len(files_pdbqt_done) < tot_dlg_files:
                      self.message_status.emit(self.message)
                      for i in os.listdir(file_dir):
                          if i.endswith(".pdbqt") and i not in files_pdbqt_done:
                             files_pdbqt_done.append(i)
                             prog = len(files_pdbqt_done) / float(tot_dlg_files)
                             self.prog_status.emit(int(prog*100))
            if self.kwargs.get('param') == 'get_score':
                self.message = "step 3. Making log of docking scores"
                count = 0
                lines = []
                while len(lines) < tot_dlg_files:
                      self.message_status.emit(self.message)
                      for i in os.listdir(file_dir):
                          try:
                              if i.startswith("docking_summary_log"):
                                  with open(i) as info_file:
                                       file_info = info_file.read()
                                       if file_info.split("\n") not in lines:
                                          lines.append(file_info.split("\n"))
                                          prog = len(lines) / float(tot_dlg_files)
                                          self.prog_status.emit(int(prog*100))
                                       elif int(prog*100) >= 95:
                                          break
                          except:
                              break
            if self.kwargs.get('param') == 'get_score_redo':
                self.message = "redo step 3. Adding info to log of docking scores"
                count = 0
                lines = []
                add_lines = []
                for i in os.listdir(working_project):
                    if i.startswith("list-mols_left_to_get-info-on"):
                       file_to_open = os.path.join(working_project, i)
                       with open(file_to_open) as info_file:
                            file_info = info_file.readlines()
                            for line in file_info:
                                if line not in add_lines:
                                   add_lines.append(file_info)
                       info_file.close()
                while len(lines) < len(add_lines)+1:
                      self.message_status.emit(self.message)
                      for i in os.listdir(file_dir):
                          if i.startswith("docking_summary_log"):
                              print("file i before opening: ", i)
                              print(os.getcwd())
                              with open(i) as info_file:
                                   file_info = info_file.read()
                                   if file_info.split("\n") not in lines:
                                      lines.append(file_info.split("\n"))
                                      prog = len(lines) / float(len(add_lines)+1)
                                      self.prog_status.emit(int(prog*100))
            if self.kwargs.get('param') == 'protein_key':
                self.message = "step 3. Adding protein information to docking scores log"
                self.message_status.emit(self.message)
            if self.kwargs.get('param') == 'run_binana':
                print("Running binana to find molecule and protein interactions")
                self.message = "step 4. looking for docking interactions using Binana starting"
                self.message_status.emit(self.message)
                files_full_log_done = []
                while len(files_full_log_done) < tot_dlg_files:
                      self.message = "step 4. looking for docking interactions using Binana running"
                      self.message_status.emit(self.message)
                      for i in os.listdir(file_dir):
                          if i.endswith("_full_log.txt") and i not in files_full_log_done:
                             files_full_log_done.append(i)
                             prog = len(files_full_log_done) / float(tot_dlg_files)
                             self.prog_status.emit(int(prog*100))
            if self.kwargs.get('param') == 'collect_binana':
                self.message = "step 5. parsing logs from Binana and generating other results summary tables"
                self.message_status.emit(self.message)
                files_parsed_log_done= []
                while len(files_parsed_log_done) < tot_full_log_file:
                      for i in os.listdir(file_dir):
                          if i.endswith("parsed_log.txt") and i not in files_parsed_log_done:
                             files_parsed_log_done.append(i)
                             prog = len(files_parsed_log_done) / float(tot_full_log_file)
                             self.prog_status.emit(int(prog*100))
            if self.kwargs.get('param') == 'sankey-html':
                self.message = "step 6. generating interaction diagrams, sankey-html"
                self.message_status.emit(self.message)
                files_html_done = []
                while len(files_html_done) < tot_files:
                     for i in os.listdir(file_dir):
                         if i.endswith("sankey.html") and i not in files_html_done:
                            files_html_done.append(i)
                            prog = len(files_html_done) / float(tot_files)
                            self.prog_status.emit(int(prog*100))
                            #bar.next()
            if self.kwargs.get('param') == 'sankey-png':
                self.message = "step 6. generating interaction diagrams, sankey-png"
                self.message_status.emit(self.message)
                files_png_done= []
                while len(files_png_done) < tot_files:
                      for i in os.listdir(file_dir):
                          if i.endswith("sankey.png") and i not in files_png_done:
                             files_png_done.append(i)
                             prog = len(files_png_done) / float(tot_files)
                             self.prog_status.emit(int(prog*100))
            if self.kwargs.get('param') == 'heatmaps':
                print("generating interaction diagrams, heatmaps")
                self.message = "step 6. generating interaction diagrams, heatmaps"
                files_done = []
                while len(files_done) < tot_files:
                      self.message_status.emit(self.message)
                      for i in os.listdir(file_dir):
                          if i.endswith("heatmap.png") and i not in files_done:
                             files_done.append(i)
                             prog = len(files_done) / float(tot_files)
                             self.prog_status.emit(int(prog*100))

            if self.kwargs.get('param') == 'find_results':
                self.message = "step 7. generating model_input_image_dataframe"
                count = 0
                results_file = []
                while len(results_file) < 22:
                      self.message_status.emit(self.message)
                      for i in os.listdir(result_dir):
                          if i not in results_file:
                              results_file.append(i)
                          prog = len(results_file) / float(22)
                          self.prog_status.emit(int(prog*100))
            if self.kwargs.get('param') == 'model_run':
                self.message = "step 8. running deep-learing models"
                count = 0
                lines = []
                tot_mols = []
                tot_models = []
                for i in os.listdir(molecule_dir):
                    if i.endswith(".pdbqt") and i not in tot_mols:
                        tot_mols.append(i)
                for j in os.listdir(result_dir):
                    if j.startswith("parsed_") and j not in tot_models:
                        tot_models.append(j)
                ##I looks like the way I am writing the table_of_prediciton an new line character is put in with each set writien after the models... so I am only going to count those
                tot_pred = len(tot_models)
                while len(lines) < tot_pred:
                      self.message_status.emit(self.message)
                      for i in os.listdir(result_dir):
                          if i.startswith("table_of_prediction"):
                              with open(i) as info_file:
                                   file_info = info_file.read()
                                   if file_info.split("\n") not in lines:
                                      lines.append(file_info.split("\n"))
                                      prog = len(lines) / float(tot_pred)
                                      self.prog_status.emit(int(prog*100))
            if self.kwargs.get('param') == 'final_work':
                self.message = "step 9. Updating results, genreating heatmaps and other plots"
                count = 0
                more_results = []
                while len(more_results) < 5:
                      self.message_status.emit(self.message)
                      for i in os.listdir(result_dir):
                          if i.endswith("heatmap.png") or i.endswith("interaction_info.csv") and i not in more_results:
                              more_results.append(i)
                              prog = len(more_results) / float(5)
                              self.prog_status.emit(int(prog*100))
            if self.kwargs.get('param') == 'check_restart':
                self.message = "Checking the progress of this project"
                os.chdir(working_project)
                lines = []
                while len(lines) < 2:
                      self.message_status.emit(self.message)
                      for i in os.listdir(working_project):
                          if i.startswith("log_proccessing_results_progress"):
                             with open(i) as info_file:
                                   file_info = info_file.read()
                                   if file_info.split("\n") not in lines:
                                      lines.append(file_info.split("\n"))
                                      prog = len(lines) / 2
                                      self.prog_status.emit(int(prog*100))

            return

class restart_progress_check(QObject):
    def __init__(self, *args, **kwargs):
        QObject.__init__(self)
        self.args = args
        self.kwargs = kwargs

    def check_status(self, working_project, install_dir):
        os.chdir(working_project)
        dlg_in_dpf = os.path.join(working_project, "dpf_files/")
        results_in_dlg = os.path.join(working_project, "dlg_files/")
        list_done = open("list_of_docked.txt", "w")
        list_success = open("successfully_docked.txt", "w")
        list_ndone = open("list-to-dock.txt", "w")

        ##looking through DPF directoru for information
        ##At the end of docking Autodock will indicate if docking was completed successfully. This section will look for that information
        ##This is done in the DPF directory first to catch any files that weren't moved.  Or to determine where to start if docking hadn't finished
        for dlg_file_name in os.listdir(dlg_in_dpf):
            if dlg_file_name.endswith(".dpf"):
                 dlg_found = os.path.join(dlg_in_dpf, dlg_file_name)

                 list_ndone.write(str(dlg_found)+os.linesep)
            if dlg_file_name.endswith(".dlg"):
                try:
                    with open(os.path.join(dlg_in_dpf, dlg_file_name)) as dlg_file:
                        info = dlg_file.read()
                        outcome_good = re.search("Successful", info)
                        outcome_bad = re.search("Unsuccessful", info)
                        if outcome_good:
                            dpf_found = os.path.join(dlg_in_dpf, dlg_file_name)
                            list_done.write(str(dpf_found)+os.linesep)
                            list_success.write(str(dpf_found)+os.linesep)
                        elif outcome_bad:
                           dpf_found = os.path.join(dlg_in_dpf, dlg_file_name)
                           list_done.write(str(dpf_found)+os.linesep)
                except IOError as exc:
                    if exc.errno != errno.EISDIR: # Do not fail if a directory is found, just ignore it.
                        raise # Propagate other kinds of IOError.

        ##Then we look though DLG directory for information
        ##Same as above but now in the DLG directory to make sure all docking was completed
        for dlg_file_name in os.listdir(results_in_dlg):
            if dlg_file_name.endswith(".dlg"):
                try:
                    with open(os.path.join(results_in_dlg, dlg_file_name)) as dlg_file:
                        info = dlg_file.read()
                        outcome_good = re.search("Successful", info)
                        outcome_bad = re.search("Unsuccessful", info)
                        if outcome_good:
                            dpf_found = os.path.join(dlg_in_dpf, dlg_file_name)
                            list_done.write(str(dpf_found)+os.linesep)
                            list_success.write(str(dpf_found)+os.linesep)

                        elif outcome_bad:
                           dpf_found = os.path.join(dlg_in_dpf, dlg_file_name)
                           list_done.write(str(dpf_found)+os.linesep)

                except IOError as exc:
                    if exc.errno != errno.EISDIR: # Do not fail if a directory is found, just ignore it.
                        raise # Propagate other kinds of IOError.
                        #break
        ##I find it is best to explicitly close files
        list_done.close()
        list_ndone.close()
        list_success.close()
        new_list_done = open("list_of_docked.txt", "r").readlines()
        list_compare = open("list-to-dock.txt","r").readlines()
        success_docked = open("successfully_docked.txt", "r").readlines()
        list_left = open("list-left-to-dock.txt","w")
        newcomplete = []
        for i in sorted(new_list_done):#, sorted(list_compare)):
            newcomplete.append(i.replace(".dlg",".dpf"))
        try:
            new_filtered_list = []
            new_filtered_list.append((set(newcomplete)-set(list_compare)) | (set(list_compare)-set(newcomplete)))
        except IOError as errno:
            print("I/O error({0}): {1}".format(errno))
        except IOError as strerror:
            print("I/O error({0}): {1}".format(strerror))
        except ValueError:
            print("Could not convert data to an integer.")
        except:
            print("Unexpected error:", sys.exc_info()[0])
            raise
        for i in new_filtered_list:
            for j in i:
                list_left.write(j)

        list_left.close()
        if not os.path.isfile("log_proccessing_results_progress.txt"):
            prog_log = open("log_proccessing_results_progress.txt", "w")
        else:
            prog_log = open("log_proccessing_results_progress.txt", "a")

        list_left = open("list-left-to-dock.txt","r").readlines()
        prog_log.write("DOCKING PROGRESS: \n")

        if  len(list_left) == 0:
            prog_log.write("INFO: DOCKING RESULTS: "+str(len(list_compare))+" Molecules for docking were attempted and "+str(len(success_docked))+" molecules were successfully docked. Any remaining were unsuccessful"+"\n" )
        else:
            prog_log.write("INFO: DOCKING RESULTS: "+str(len(list_compare))+" Molecules for docking were attempted and "+str(len(success_docked))+" molecules were successfully docked. The remaining "+str(len(list_left))+" can still be docked"+"\n" )
        prog_log.close()
        #lists of files found
        dlg_found_dlg_dir = []
        pdbqt_found_dlg_dir = []
        parsed_log_dlg_dir = []
        parsed_heatmap_dlg_dir = []
        full_log_dlg_dir = []
        sankey_html_dlg_dir = []
        sankey_png_dlg_dir = []
        dlg_unsuccessfull = []
        #look for files
        for file in os.listdir(results_in_dlg):
            if file.endswith(".dlg"):
                try:
                    with open(os.path.join(results_in_dlg, file)) as dlg_file: # No need to specify 'r': this is the default.
                        info = dlg_file.read()
                        outcome_good = re.search("Successful", info)
                        outcome_bad = re.search("Unsuccessful", info)
                        dlg_file.close()
                        if outcome_good:
                           dlg_found_dlg_dir.append(file.replace(".dlg",""))
                        elif outcome_bad:
                           dlg_unsuccessfull.append(file.replace(".dlg",""))
                except IOError as exc:
                    if exc.errno != errno.EISDIR: # Do not fail if a directory is found, just ignore it.
                        raise # Propagate other kinds of IOError.
                        #break
                    else:
                        pass

        for file in os.listdir(results_in_dlg):
            if file.endswith(".pdbqt"):
                pdbqt_found_dlg_dir.append(file.replace(".pdbqt",""))
            elif file.endswith("_full_log.txt"):
                full_log_dlg_dir.append(file.replace("_full_log.txt",""))
            elif file.endswith("_parsed_log.txt"):
                parsed_log_dlg_dir.append(file.replace("_parsed_log.txt",""))
            elif file.endswith("_parsed_heatmap.png"):
                parsed_heatmap_dlg_dir.append(file.replace("_parsed_heatmap.png",""))
            elif file.endswith("_sankey.html"):
                sankey_html_dlg_dir.append(file.replace("_sankey.html",""))
            elif file.endswith("_sankey.png"):
                sankey_png_dlg_dir.append(file.replace("_sankey.png",""))

        for list in ['pdbqt_found_dlg_dir','parsed_log_dlg_dir','parsed_heatmap_dlg_dir','full_log_dlg_dir', 'sankey_html_dlg_dir','sankey_png_dlg_dir']:
            if not pdbqt_found_dlg_dir:
                pdbqt_found_dlg_dir.append("no data pdbqt")
            elif not full_log_dlg_dir:
                full_log_dlg_dir.append("no data full log")
            elif not parsed_log_dlg_dir:
                parsed_log_dlg_dir.append("no data parsed log")
            elif not parsed_heatmap_dlg_dir:
                parsed_heatmap_dlg_dir.append("no data heatmap")
            elif not sankey_html_dlg_dir:
                sankey_html_dlg_dir.append("no data sankey html")
            elif not sankey_png_dlg_dir:
                sankey_png_dlg_dir.append("no data sankey png")
            else:
                break


        dockinglog = os.path.join(results_in_dlg, "docking_summary_log.txt")
        docking_info = pd.read_csv(dockinglog, delimiter=',', header=0, engine='python')
        molecule_name = []
        receptor_name = []
        dlgf_from_docking_log = []
        dlg_protmol_key = []
        ##The next part looks through the docking log to see if is complete
        for dlg_file_name in docking_info['lowestEnergy_dlgfn']:
            if dlg_file_name not in dlgf_from_docking_log:
                dlgf_from_docking_log.append(dlg_file_name)
        del docking_info
        protmol_log = os.path.join(results_in_dlg, "docking_summary_protein_key.txt")
        protmol_info = pd.read_csv(protmol_log, delimiter=',', header=0, engine='python')
        for protkey_info_dlg, mol_name, prot_name in zip(protmol_info['lowestEnergy_dlgfn'], protmol_info['mol_name'], protmol_info['protein_id']):
            if protkey_info_dlg not in dlg_protmol_key:
                dlg_protmol_key.append(protkey_info_dlg)
                #sometime names get chopped down but once they make it into a unique dpf file they will be maintained by adding numbers to the end, however that won't work for this method of appending the list so i am going to associate the molecule name with the dlg name
                molecule_name.append(mol_name)
            if prot_name not in receptor_name:
                receptor_name.append(prot_name)
        del protmol_info
        protein_dir = os.path.join(install_dir, "protein_files/")
        proteins_found = []
        for prot_file, subdirs, files in os.walk(protein_dir):
            for prot_file in files:
                if prot_file.endswith(".pdbqt"):
                    proteins_found.append(prot_file)
        ##I need to change this to a dictionary
        list_from_dlg_dir = [dlg_found_dlg_dir, pdbqt_found_dlg_dir, dlgf_from_docking_log, molecule_name, receptor_name, full_log_dlg_dir, parsed_log_dlg_dir, parsed_heatmap_dlg_dir, sankey_html_dlg_dir, sankey_png_dlg_dir]
        list_id_name = ["dlg files", "pdbqt files (extracted from docking logs)", "docking summary log from dlg files", "molecule entries in next log", "protein entries in next log", "full interaction logs from binana", "parsed interaction logs", "parsed interaction heatmaps", "sankey interaction diagrams(html)", "sankey interaction diagrams(png)"]

        def gen_list_left(dlg_found_dlg_dir, list, **kwargs):
            if kwargs.get('param') == 'pdbqt files (extracted from docking logs)':
                new_filtered_list = []
                if 'no data pdbqt' in list:
                    list.remove('no data pdbqt')
                else:
                    pass
                new_filtered_list.append((set(list)-set(dlg_found_dlg_dir)) | (set(dlg_found_dlg_dir)-set(list)))# (list(list(set(parsed_heatmap_dlg_dir)-set(dlgf_from_docking_log)) + list(set(dlgf_from_docking_log)-set(parsed_heatmap_dlg_dir))))
                list_left = open("list-mols_left_to_extract.txt","w")
                for i in new_filtered_list:
                    for j in i:
                        list_left.write(j+".dlg"+os.linesep)
            elif kwargs.get('param') == 'docking summary log from dlg files':
                new_filtered_list = []
                if 'no data' in list:
                    list.remove('no data')
                else:
                    pass
                new_filtered_list.append((set(list)-set(dlg_found_dlg_dir)) | (set(dlg_found_dlg_dir)-set(list)))# (list(list(set(parsed_heatmap_dlg_dir)-set(dlgf_from_docking_log)) + list(set(dlgf_from_docking_log)-set(parsed_heatmap_dlg_dir))))
                list_left = open("list-mols_left_to_get-info-on.txt","w")
                for i in new_filtered_list:
                    for j in i:
                        list_left.write(j+".dlg"+os.linesep)
                list_left.close()
            elif kwargs.get('param') == 'full interaction logs from binana':
                new_filtered_list = []
                if 'no data full log' in list:
                    list.remove('no data full log')
                else:
                    pass
                new_filtered_list.append((set(list)-set(dlg_found_dlg_dir)) | (set(dlg_found_dlg_dir)-set(list)))
                list_left = open("list-mols_left_for_binana.txt","w")
                for i in new_filtered_list:
                    for j in i:
                        #print(j)
                        file_plus_path = os.path.join(results_in_dlg,j)
                        list_left.write(file_plus_path+".dlg"+os.linesep)
            elif kwargs.get('param') == "parsed interaction logs":
                new_filtered_list = []
                if 'no data parsed log' in list:
                    list.remove('no data parsed log')
                else:
                    pass
                new_filtered_list.append((set(list)-set(dlg_found_dlg_dir)) | (set(dlg_found_dlg_dir)-set(list)))
                list_left = open("list-mols_left_for_counting.txt","w")
                for i in new_filtered_list:
                    for j in i:
                        #print(j)
                        file_plus_path = os.path.join(results_in_dlg,j)
                        list_left.write(file_plus_path+"_full_log.txt"+os.linesep)
            elif kwargs.get('param') == 'parsed interaction heatmaps':
                new_filtered_list = []
                if 'no data heatmap' in list:
                    list.remove('no data heatmap')
                else:
                    pass
                new_filtered_list.append((set(list)-set(dlg_found_dlg_dir)) | (set(dlg_found_dlg_dir)-set(list)))# (list(list(set(parsed_heatmap_dlg_dir)-set(dlgf_from_docking_log)) + list(set(dlgf_from_docking_log)-set(parsed_heatmap_dlg_dir))))
                list_left = open("list-left-to-process_heatmaps.txt","w")
                for i in new_filtered_list:
                    for j in i:
                        file_plus_path = os.path.join(results_in_dlg,j)
                        list_left.write(file_plus_path+"_parsed_log.txt"+os.linesep)
            elif  kwargs.get('param') == 'sankey interaction diagrams(png)':
                new_filtered_list = []
                if 'no data sankey png' in list:
                    list.remove('no data sankey png')
                else:
                    pass
                new_filtered_list.append((set(list)-set(dlg_found_dlg_dir)) | (set(dlg_found_dlg_dir)-set(list)))
                list_left = open("list-left-to-process_sankey-png.txt","w")
                for i in new_filtered_list:
                    for j in i:
                        #print(j)
                        file_plus_path = os.path.join(results_in_dlg,j)
                        list_left.write(file_plus_path+"_parsed_log.txt"+os.linesep)
            elif kwargs.get('param') == 'sankey interaction diagrams(html)':
                new_filtered_list = []
                if 'no data sankey html' in list:
                    list.remove('no data sankey html')
                else:
                    pass
                new_filtered_list.append((set(list)-set(dlg_found_dlg_dir)) | (set(dlg_found_dlg_dir)-set(list)))
                list_left = open("list-left-to-process_sankey-html.txt","w")
                for i in new_filtered_list:
                    for j in i:
                        #print(j)
                        file_plus_path = os.path.join(results_in_dlg,j)
                        list_left.write(file_plus_path+"_parsed_log.txt"+os.linesep)
            return

        if not os.path.isfile("log_proccessing_results_progress.txt"):
            prog_log = open("log_proccessing_results_progress.txt", "w")
        else:
            prog_log = open("log_proccessing_results_progress.txt", "a")
        i = 0
        prog_log.write("\nPROCESSING DOCKING RESULTS PROGRESS: \n")
        for list in list_from_dlg_dir:

            if len(list) == 0:
                # the first variable was list_id_name[i] and changed to str(list_id_name[list_from_dlg_dir.index(list)]) ### I think i can delete this comment the change worked and it works as expected as it now!
                prog_log.write("INFO: No "+str(list_id_name[list_from_dlg_dir.index(list)])+" are found, there should be "+str(len(success_docked))+"\n")
                i = i + 1
            #These are compared to the number successful because they will probably be less than total attmepted
            elif len(list) == len(success_docked) and not list_from_dlg_dir.index(list) == 4:
                prog_log.write("INFO: generating "+str(list_id_name[list_from_dlg_dir.index(list)])+" has finished, there are: "+str(len(list))+ " and should be: "+str(len(success_docked))+"\n")
                # what is going to happen if the following ins't there
                i = i + 1
            elif len(list) < len(success_docked) and not list_from_dlg_dir.index(list) == 4:
                gen_list_left(sorted(dlg_found_dlg_dir), sorted(list), param=str(list_id_name[list_from_dlg_dir.index(list)]))
                prog_log.write("INFO: generating "+str(list_id_name[list_from_dlg_dir.index(list)])+" hasn't finished, so far only: "+str(len(list))+ " there should be "+str(len(success_docked))+"\n")
                i = i + 1
               # what is this statement going to do if to happen if "4" ins't there or no proteins are found
            elif list_from_dlg_dir.index(list) == 4 and len(list) == len(proteins_found):
                output = str(list)
                prog_log.write("INFO: generating new line "+str(list_id_name[list_from_dlg_dir.index(list)])+" has finished, there are: "+str(len(list))+ " and should be: "+str(len(proteins_found))+"\n")
                prog_log.write("INFO: Receptors docked on: " +str(receptor_name).strip("['").strip("\']")+"\n")
                i = i + 1
                continue
        prog_log.close()
        return

class process_input(QObject):
    prog_status = pyqtSignal(int)
    message_status = pyqtSignal(str)
    def __init__(self, *args, **kwargs):
        QObject.__init__(self)
        #self.working_project = working_project
        self.args = args
        self.kwargs = kwargs

    def search_pubchem(self, compounds_dict):
            cid_array = []
            #cmpd_list = []
            #num = 0
            for name, cid in compounds_dict.items():
                if cid != [0] and len(cid) == 1:
                    cid_array.append(str(cid).replace("[","").replace("]",""))
                    try:
                       print("looking for cid ", cid)
                       self.message ="looking through PubChem for cid: "+str(cid)
                       self.message_status.emit(self.message)
                       lookup = pcp.download('SDF', name +"_pubchem-version"+'.sdf', cid, record_type='3d')
                    except:
                        try:
                            print("Request not found error, likely no 3D structure. \n Trying for any (probably 2D) structure")
                            self.message = "Request not found error, likely no 3D structure. \n Trying for any (probably 2D) structure"
                            self.message_status.emit(self.message)
                            pcp.download('SDF', name +"_pubchem-version_2D"+'.sdf', cid)
                        except:
                            pass
                    else:
                       pass
                elif cid != [0] and len(cid) > 1:
                    for i in cid:
                         print("There more than one CID number this is one: ", name)
                         cid_array.append(str(i).replace("[","").replace("]",""))
                    try:
                       print("looking for cid ", cid)
                       self.message ="looking through PubChem for cid: "+str(cid)
                       self.message_status.emit(self.message)
                       lookup = pcp.download('SDF', name+"_pubchem-version"+'.sdf', cid, record_type='3d')
                    except:
                        try:
                            print("Request not found error, likely no 3D structure. \n Trying for any (probably 2D) structure")
                            self.message = "Request not found error, likely no 3D structure. \n Trying for any (probably 2D) structure"
                            self.message_status.emit(self.message)
                            pcp.download('SDF', name+"_pubchem-version_2D"+'.sdf', cid)
                        except:
                            pass
                    else:
                       pass
            ###despite the fact that it looks like an array of CID numbers is require it can actaully be a list
            pcp.download('CSV', 'compounds_found.csv', cid_array, operation='property/MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,InChI,InChIKey,IUPACName,XLogP,ExactMass,TPSA')

    def check_molecules(self, molecule_files, compounds_dict):
        for sdf_output_file in os.listdir(molecule_files):
            if sdf_output_file.endswith(".sdf"):
                suppl = Chem.SDMolSupplier(sdf_output_file)
                temp_file = sdf_output_file.replace(".sdf","_temp.sdf").replace(" ","")
                writer = Chem.SDWriter(temp_file)
                for mol in suppl:
                    print("this is prop: ", mol.GetProp('_Name'))
                    nam = mol.GetProp('_Name').replace("-","_")
                    prop = Chem.PropertyMol.PropertyMol(mol)
                    #prop.SetProp('image', str(nam)+".png")
                    mwt = Descriptors.MolWt(mol)
                    tpsa = Descriptors.TPSA(mol)
                    logp = Descriptors.MolLogP(mol)
                    rotbond = Descriptors.NumRotatableBonds(mol)
                    hbdaccept = Descriptors.NumHAcceptors(mol)
                    hbddon = Descriptors.NumHDonors(mol)
                    fractCSP3 = Descriptors.FractionCSP3(mol)
                    prop.SetProp('MolWt', str(mwt))
                    prop.SetProp('tpsa', str(tpsa))
                    prop.SetProp('logp', str(logp))
                    prop.SetProp('rotatable_bonds', str(rotbond))
                    prop.SetProp('hbond_accept', str(hbdaccept))
                    prop.SetProp('hbond_donor', str(hbddon))
                    prop.SetProp('fractionCSP3', str(fractCSP3))
                    if mol.GetConformer().Is3D() == False:
                        print("this molecule is 2D: ", sdf_output_file)
                        m2 = Chem.AddHs(mol)
                        ids = AllChem.EmbedMultipleConfs(m2, numConfs=1, params=AllChem.ETKDG())
                        mwt = Descriptors.MolWt(m2)
                        tpsa = Descriptors.TPSA(m2)
                        logp = Descriptors.MolLogP(m2)
                        rotbond = Descriptors.NumRotatableBonds(m2)
                        hbdaccept = Descriptors.NumHAcceptors(m2)
                        hbddon = Descriptors.NumHDonors(m2)
                        fractCSP3 = Descriptors.FractionCSP3(m2)
                        m2.SetProp('MolWt', str(mwt))
                        m2.SetProp('tpsa', str(tpsa))
                        m2.SetProp('logp', str(logp))
                        m2.SetProp('rotatable_bonds', str(rotbond))
                        m2.SetProp('hbond_accept', str(hbdaccept))
                        m2.SetProp('hbond_donor', str(hbddon))
                        m2.SetProp('fractionCSP3', str(fractCSP3))
                        m2.SetProp('_Name', sdf_output_file.replace("pubchem-version_2D.sdf","")+str(nam))
                        m2.SetProp('image', sdf_output_file.replace("pubchem-version_2D.sdf","")+str(nam)+".png")
                        temp_file_2D = sdf_output_file.replace("_2D.sdf","_temp.sdf")
                        writer2D = Chem.SDWriter(temp_file_2D)
                        for id in ids:
                            writer2D.write(m2, confId=id)
                        tmp=AllChem.Compute2DCoords(mol)
                        pic_file_name = sdf_output_file.replace("pubchem-version_2D.sdf","")+str(nam)+".png"
                        Draw.MolToFile(mol, pic_file_name)
                    elif "2D.sdf" not in sdf_output_file:
                        prop.SetProp('_Name', sdf_output_file.replace("pubchem-version.sdf","")+str(nam))
                        prop.SetProp('image', sdf_output_file.replace("pubchem-version.sdf","")+str(nam)+".png")
                        writer.write(prop)
                        tmp=AllChem.Compute2DCoords(mol)
                        pic_file_name = sdf_output_file.replace("pubchem-version.sdf","")+str(nam)+".png"
                        Draw.MolToFile(mol, pic_file_name)
        #writer.close()
        #writer2D.close()
        return

    def concatinate_files(self, molecule_files):
        concat_file_list = []
        for file in os.listdir(molecule_files):
            if file.endswith("temp.sdf") and file not in concat_file_list:
                concat_file_list.append(file)

        with open("all_molecules_generated.sdf", "w") as outfile:
            for filename in concat_file_list:
                with open(filename, "r") as infile:
                    line_data = infile.read()
                    outfile.write(line_data)
        #outfile.close()
        #infile.close()
        return
    def generate_output_table(self, final_output_file):
        frame = PandasTools.LoadSDF(final_output_file,
                                         smilesName='SMILES',
                                         includeFingerprints=False)
        frame.set_index(['ID'], inplace=True)
        frame.to_csv('excel_molecules.csv', sep=':')
        return
    def start_the_search(self, input_file):
            chemicals = pd.read_csv(input_file, sep=',', header=0, index_col=False, engine='python')
            compounds_dict = {}
            for i, j in zip(chemicals['Compound'], chemicals['Name']) :
                if Chem.MolFromSmiles(i):
                    format = "smiles"
                    #self.search_pubchem(i, format, j)
                    cid = pcp.get_cids(i, format, list_return='flat')
                    if cid and j not in compounds_dict:
                        compounds_dict[j] = cid
                elif Chem.inchi.MolFromInchi(i):
                    format = "inchi"
                    #self.search_pubchem(i, format, j)
                    cid = pcp.get_cids(i, format, list_return='flat')
                    if cid and j not in compounds_dict:
                        compounds_dict[j] = cid
                elif not Chem.inchi.MolFromInchi(i) and not Chem.MolFromSmiles(i):
                    format = "name"
                    #self.search_pubchem(i, format, name)
                    cid = pcp.get_cids(i, format, list_return='flat')
                    if cid and j not in compounds_dict:
                        compounds_dict[j] = cid
                    print("guessing format is name")
            self.search_pubchem(compounds_dict)
            self.check_molecules(molecule_files, compounds_dict)
            self.concatinate_files(molecule_files)
            final_output_file = "all_molecules_generated.sdf"
            self.generate_output_table(final_output_file)
            os.mkdir("startup_files_from_table_"+datetime.datetime.now().strftime("%m-%d-%y"))
            for file in os.listdir(molecule_files):
                if file.endswith(".png") or file.endswith(".csv"):
                    shutil.move(file, working_project)
                elif not file.startswith("all_molecules_generated"):
                    shutil.move(file, "./startup_files_from_table_"+datetime.datetime.now().strftime("%m-%d-%y"))
                else:
                    pass
            self.message = "Gathering molecules has finished. \n Setup docking in tools menu is the next step"
            self.message_status.emit(self.message)

class process_results(QObject):
    prog_status = pyqtSignal(int)
    message_status = pyqtSignal(str)
    def __init__(self, *args, **kwargs):
        QObject.__init__(self)
        #self.working_project = working_project
        self.args = args
        self.kwargs = kwargs

    def sumarize_results_and_generate_tables(self, dlg_files, protein_files):
            list_of_logs = []
            list_of_proteins = []
            list_of_sankey = []
            for mol_log_location in os.listdir(dlg_files):
                if mol_log_location.endswith("parsed_log.txt"):
                    sankey_full_files = os.path.join(dlg_files, mol_log_location)
                    list_of_logs.append(mol_log_location)
            for dir, subdirs, files in os.walk(protein_files):
                for name in files:
                    protein_name = os.fsdecode(name)
                    if protein_name.endswith(".pdbqt"):
                       list_of_proteins.append(protein_name.replace(".pdbqt",""))

            for sankey_pics_location in os.listdir(dlg_files):
                if sankey_pics_location.endswith("sankey.png"):
                    sankey_full_files = os.path.join(dlg_files, sankey_pics_location)
                    list_of_sankey.append(sankey_pics_location)
            k=-1
            hbonds = []
            hydophob = []
            pipi = []
            tstack = []
            picat = []
            sbridge = []
            ori_log = []
            curr_protein = []
            curr_mol = []
            report_mol_name = []
            merge_col = []
            for j in list_of_proteins:
               k = k +1
               for i, elem in enumerate(list_of_logs):
                if list_of_proteins[k] in elem:

                    file_with_path = os.path.join(dlg_files, elem)
                    log_contents = pd.read_csv(file_with_path, index_col=0, sep=',', engine='python')
                    file = open(file_with_path, 'r')
                    contents = file.read()
                    protein_name = list_of_proteins[k]

                    step_one = str(elem).replace("_parsed_log.txt", "")

                    mol_name_ori = str(step_one).replace(list_of_proteins[k]+"_","")
                    #this is a little silly but with "_" in the molecule name it is too long and connot to text wrap. the same is probably true for "-" but I will have to fis that later
                    if str("_") in mol_name_ori:
                        mol_name = mol_name_ori.replace("_", " ")
                    else:
                        mol_name = mol_name_ori
                    merge_mol_prot = protein_name+"_"+mol_name_ori
                    ori_log.append(elem)
                    curr_protein.append(protein_name)
                    curr_mol.append(mol_name)
                    report_mol_name.append(mol_name_ori)
                    merge_col.append(merge_mol_prot)
                    with open(file_with_path) as search:
                        for line in search:
                           line = line.rstrip()
                           if re.match("Hydrogen bonds", line):
                               hbonds.append(line.replace("Hydrogen bonds,","").replace(","," "))
                           if re.match("Hydrophobic contacts", line):
                               hydophob.append(line.replace("Hydrophobic contacts,","").replace(","," "))
                           if re.match("pi-pi stacking interactions", line):
                               pipi.append(line.replace("pi-pi stacking interactions,","").replace(","," "))
                           if re.match("T-stacking", line):
                               tstack.append(line.replace("T-stacking,","").replace(","," "))
                           if re.match("cation-pi", line):
                               picat.append(line.replace("cation-pi,","").replace(","," "))
                           if re.match("Salt Bridges", line):
                               sbridge.append(line.replace("Salt Bridges,","").replace(","," "))

            report_data = pd.DataFrame({
                        "Protein Name" : curr_protein,
                        "Report Molecule Name" : report_mol_name,
                        "Molecule Name" : curr_mol,
                        "Hydrogen bonds" : hbonds,
                        "pi-pi stacking interactions" : pipi,
                        "T-stacking" : tstack,
                        "cation-pi" : picat,
                        "Salt Bridges" : sbridge,
                        "merge_column": merge_col})
            report_data.set_index(['Protein Name'], inplace=True)
            if not os.path.isfile('table_of_interaction_info.csv'):
                report_data.to_csv('table_of_interaction_info.csv', mode='w')
            else:
                report_data.to_csv('table_of_interaction_info.csv', mode='a', header=False)
            return

    def plot_merg_scores_and_prediction(self, scores_info_file, prediction_info_file):
        docking_scores = pd.read_csv(scores_info_file, delimiter=',', header=0, engine='python')
        prediction_info = pd.read_csv(prediction_info_file, delimiter=',', header=0, engine='python')

        full_results_log = pd.merge(left=docking_scores, right=prediction_info, how='outer', left_on='lowestEnergy_dlgfn', right_on='Filename')
        full_results_log.fillna(0, inplace=True)

        perdiction_per_score = []
        #here I am "zipping" through the columns selected and iterating them at the same rate to find info to append to perdiction_per_score list then adding that to colum "Score_per_Prediction"
        for i, j, k, l in zip(full_results_log['Score_per_Mass'], full_results_log['Predictions'], full_results_log['probability'], full_results_log['Predciton_heatmap_format']):
                if i <= -0.025 and j <= 2 and k >= 70:
                     perdiction_per_score.append(l)
                elif j == 0 and k == 0:
                     perdiction_per_score.append(5)
                else:
                     perdiction_per_score.append(3)

        full_results_log['Score_per_Prediction'] = perdiction_per_score
        full_results_log.to_csv("docking_summary_and_prediction_key.txt")

        final_results = pd.DataFrame(full_results_log[full_results_log['Score_per_Prediction']<5])
        final_results = final_results.sort_values(by=['Protein_name'])
        #final_results[:14].apply(str).str.replace("-","_") I am not sure I need this but it was how I accessed the y values or labels
        fig = go.Figure(go.Heatmap(z=final_results['Predciton_heatmap_format'].values,x=final_results['Protein_name'],y=final_results['Molecule'].apply(str).str.replace("-","_"),zmax=5,zmin=0,colorscale=[[0.0, 'rgb(205, 38, 38)'], [0.02, 'rgb(255, 255, 204)'], [0.22, 'rgb(205, 38, 38)'], [0.22, 'rgb(255, 255, 204)'],
            [0.41, 'rgb(205, 38, 38)'], [0.41, 'rgb(255, 255, 204)'], [0.6, 'rgb(205, 38, 38)'], [0.6, 'rgb(51, 161, 201)'], [0.8, 'rgb(51, 161, 201)'], [1.0, 'rgb(51, 161, 201)']]))#["lavender", "pink", "purple", "lightblue", "lightblue", "lightblue"]))
        py.offline.plot(fig, auto_open=False, output_type='file', filename="predictions_heatmap.html")
        py.io.write_image(fig, "predictions_heatmap.png", format='png', scale=None, width=1000, height=800)

        fig1 = go.Figure(go.Heatmap(z=final_results['Score_per_Prediction'],x=final_results['Protein_name'],y=final_results['Molecule'].apply(str).str.replace("-","_"),zmax=5,zmin=0,colorscale=[[0.0, 'rgb(205, 38, 38)'], [0.02, 'rgb(255, 255, 204)'], [0.22, 'rgb(205, 38, 38)'], [0.22, 'rgb(255, 255, 204)'],
            [0.41, 'rgb(205, 38, 38)'], [0.41, 'rgb(255, 255, 204)'], [0.6, 'rgb(205, 38, 38)'], [0.6, 'rgb(51, 161, 201)'], [0.8, 'rgb(51, 161, 201)'], [1.0, 'rgb(51, 161, 201)']]))#["lavender", "pink", "purple", "lightblue", "lightblue", "lightblue"]))
        py.offline.plot(fig1, auto_open=False, output_type='file', filename="score_predictions_heatmap.html")
        py.io.write_image(fig1, "score_predictions_heatmap.png", format='png', scale=None, width=1000, height=800)

        fig2 = go.Figure(go.Heatmap(z=full_results_log['LE'].values,x=full_results_log['protein_id'],y=full_results_log['mol_name'].apply(str).str.replace("-","_"),zmax=0,zmin=full_results_log['LE'].values.min(),colorscale='RdYlBu'))
        py.offline.plot(fig2, auto_open=False, output_type='file', filename="docking_summary_protein_key_heatmap.html")
        py.io.write_image(fig2, "docking_summary_protein_key_heatmap.png", format='png', scale=None, width=1000, height=800)

        fig3 = go.Figure(go.Heatmap(z=full_results_log['Score_per_Mass'].values,x=full_results_log['protein_id'],y=full_results_log['mol_name'].apply(str).str.replace("-","_"),zmax=0,zmin=full_results_log['Score_per_Mass'].values.min(),colorscale='RdYlBu'))
        py.offline.plot(fig3, auto_open=False, output_type='file', filename="docking_summary_protein_key_score-per-mass_heatmap.html")
        py.io.write_image(fig3, "docking_summary_protein_key_score-per-mass_heatmap.png", format='png', scale=None, width=1000, height=800)

    def make_figures(self, table_of_prediction):
        new_data_frame = pd.read_csv(table_of_prediction, delimiter=',')

        fig = go.Figure(go.Heatmap(z=new_data_frame.iloc[:, 5].values,x=new_data_frame.iloc[:, 2].values,y=new_data_frame.iloc[:,1].apply(str).str.replace("-","_"),zmax=5,zmin=0,colorscale=[[0.0, 'rgb(205, 38, 38)'], [0.02, 'rgb(255, 255, 204)'], [0.22, 'rgb(205, 38, 38)'], [0.22, 'rgb(255, 255, 204)'],
            [0.41, 'rgb(205, 38, 38)'], [0.41, 'rgb(255, 255, 204)'], [0.6, 'rgb(205, 38, 38)'], [0.6, 'rgb(51, 161, 201)'], [0.8, 'rgb(51, 161, 201)'], [1.0, 'rgb(51, 161, 201)']]))#["lavender", "pink", "purple", "lightblue", "lightblue", "lightblue"]))
        py.offline.plot(fig, auto_open=False, output_type='file', filename="predictions_heatmap.html")
        py.io.write_image(fig, "predictions_heatmap.png", format='png', scale=None, width=1200, height=1000)

    def run_models(self, dframe, PATH, model_to_load, model_type_match, results_files):
        IMG_HEIGHT = 210
        IMG_WIDTH = 210
        prediction_image_generator = ImageDataGenerator(rescale=1./255)
        pred_data_gen = prediction_image_generator.flow_from_dataframe(dataframe=dframe,
                                                                       directory=PATH,
                                                                       x_col="full_file_name",
                                                                       y_col=None,
                                                                       batch_size=1,
                                                                       shuffle=False,
                                                                       target_size=(IMG_HEIGHT, IMG_WIDTH),
                                                                       class_mode=None)
        #######################################################################################################
        ### Tensorflow model Predictions mode
        #######################################################################################################
        tf.keras.backend.clear_session()
        model_X = load_model(str(model_to_load))
        predictionstenflo = model_X.predict(pred_data_gen, batch_size=1, verbose=0)
        score = tf.nn.softmax(predictionstenflo,axis=1)
        proba = score.numpy()
        predict_class_indices=np.argmax(predictionstenflo,axis=1).tolist()
        filenames=pred_data_gen.filenames
        mol_name = []
        short_filenames = []
        for elem in filenames:
            step_one = str(elem).replace("_parsed_heatmap.png", "")
            step_two = str(step_one).replace(model_type_match+"_","")
            mol_name.append(step_two)
            short_filenames.append(str(elem).replace("_parsed_heatmap.png", ""))
        results = pd.DataFrame({"Filename" : short_filenames,
                           "Molecule" : mol_name,
                           "Protein_name" : model_type_match,
                           "Predictions" : predict_class_indices,
                           "probability" : 100*np.max(proba,axis=1),
                           "Predciton_heatmap_format" : predict_class_indices+np.max(proba,axis=1)-0.001 })
        results.set_index(['Filename'], inplace=True)
        if not os.path.isfile('table_of_prediction.csv'):
            results.to_csv('table_of_prediction.csv', mode='w')
        else:
            results.to_csv('table_of_prediction.csv', mode='a', header=False)
        tf.keras.backend.clear_session()
        return

    def start_predictions(self, install_dir, results_files, dlg_files):
       tensor_model_dir = os.path.join(install_dir, "image_tensorflow_models/")
       if platform == "win32":
          print("windows was found")
          tenflow_dict = os.path.join(install_dir,"tensorflow_dictionary_hard-copy.txt")
       elif platform == "darwin":
          print("MACOS was found")
          tenflow_dict = os.path.join(install_dir,"tensorflow_dictionary_hard-copy.txt")
       elif platform == "linux" or "linux2":
          print("linux was found")
          tenflow_dict = os.path.join(install_dir,"tensorflow_dictionary_hard-copy.txt")
          print("this is pythonsh_loc from setup_Autodocking script: ", tenflow_dict)
      ## Right NOW 2-4-2021 this has to be in results_files and I want it in install_dir
       with open(str(tenflow_dict), 'rb') as handle:
           data = handle.read()
       model_dict = json.loads(data)
       for files in os.listdir():
          if files.startswith("parsed_"):
             dframe = pd.read_csv(files, header=0, low_memory=False)
             model_type_match = files.lstrip("parsed_").rstrip(".csv")
             PATH = dlg_files+"/"
             ## I need to change this path to install_dir + image_tensorflow_models
             ## I made the change on 2-4-2021 but I am not sure it  is final
             model_to_load = os.path.join(tensor_model_dir, model_dict.get(files.lstrip("parsed_").rstrip(".csv")))
             self.run_models(dframe, PATH, model_to_load, model_type_match, results_files)
          else:
             continue
       return

    def generate_deep_learning_input(self, dlg_files, results_files):
           first = True
           for find_parsed_heatmap_file in os.listdir(dlg_files):
                if find_parsed_heatmap_file.endswith("parsed_heatmap.png") and first==True:
                   first = False
                   with open("model_input_image_dataframe.csv", "w") as molecule_prot_info:
                              colnames = ['full_path', 'full_file_name']
                              writer = csv.DictWriter(molecule_prot_info, fieldnames=colnames)
                              writer.writeheader()
                              writer.writerow({'full_path' : os.getcwd(), 'full_file_name' : find_parsed_heatmap_file})
                elif find_parsed_heatmap_file.endswith("parsed_heatmap.png"):
                   with open("model_input_image_dataframe.csv", "a") as molecule_prot_info:
                              colnames = ['full_path', 'full_file_name']
                              writer = csv.DictWriter(molecule_prot_info, fieldnames=colnames)
                              writer.writerow({'full_path' : os.getcwd(),'full_file_name' : find_parsed_heatmap_file})
           ### added 3-17-2021
           molecule_prot_info.close()
           for find_image_dataframe in os.listdir(dlg_files):
               if find_image_dataframe.endswith("dataframe.csv"):
                  shutil.copy(find_image_dataframe, results_files)
               elif find_image_dataframe.endswith("key.txt"):
                  shutil.copy(find_image_dataframe, results_files)
               elif find_image_dataframe.endswith("summary_log.txt"):
                  shutil.copy(find_image_dataframe, results_files)
               else:
                  continue
           os.chdir(results_files)
           if platform == "win32":
               print("windows was found")
               tenflow_dict = os.path.join(install_dir,"tensorflow_dictionary_hard-copy.txt")
           elif platform == "darwin":
               print("MACOS was found")
               tenflow_dict = os.path.join(install_dir,"tensorflow_dictionary_hard-copy.txt")
           elif platform == "linux" or "linux2":
               print("linux was found")
               tenflow_dict = os.path.join(install_dir,"tensorflow_dictionary_hard-copy.txt")
           with open(str(tenflow_dict), 'rb') as handle:
                data = handle.read()
           # reconstructing the data as a dictionary
           model_dict = json.loads(data)
           # make a list of keys
           checkfor = []
           for key, values in model_dict.items():
               checkfor.append(key)

           for file in os.listdir(results_files):
                 if file.startswith("model_input_image_dataframe"):
                    masterfile = pd.read_csv(file, header=0, low_memory=False)
                 else:
                     continue

           for line in checkfor:
               find_lines = masterfile[masterfile.full_file_name.apply(str).str.contains(line)]
               find_lines.to_csv('parsed_'+line+'.csv', index=False)

    def run_summarize_interactions(self, install_dir, prot_def, lig):
        binana_loc = os.path.join(install_dir, "binana.py")
        subprocess.run(["python",str(binana_loc),"-receptor",str(prot_def),"-ligand",str(lig)],stdout=open(lig.replace('.pdbqt','_full_log.txt'),"w"))
        #os.system("mpirun --use-hwthread-cpus -n 20 python "+str(binana_loc)+" -receptor "+str(prot_def)+" -ligand "+str(lig) +" -output_file "+ str(lig.replace('.pdbqt','_full_log.txt')))
        return

    def rename_file(self, renamefile):
       path, filename = os.path.split(renamefile)
       os.chdir(path)

       try:
          with open(filename) as line:
            # this is reading in the contents of the dlg file with the .read() function and that is stored in the "data" variable
           data = line.read()
           # the regular expression search function is going to look for the ligand or docked molecule name using the input search pattern.
           lig_name = re.search("INPUT-LIGAND-PDBQT: REMARK  Name = ([\w-]+)", data)
           prot_name = re.search("Macromolecule file used to create Grid Maps =([\t\w-]+)", data)
           line.close()
          # The next part will use the infromation saved in "lig_name" to rename the file
           if lig_name:
             C = lig_name.group(1)
             B = prot_name.group(1).strip()
             #settig up new name variables
             newname = B+"_"+C+".dlg"
             newname_add = ""
             newname_add2 = ""
             #The next part searches for new names and old namaes and comes up with unique new names
             if os.path.isfile(filename):
                add = 0
                while True:
                    add += 1
                    newname_add = B+"_"+C+"_"+str(add)+".dlg"
                    if os.path.isfile(newname_add):
                      #This checks for the final name first. I know this is a silly order but this is what works
                      continue
                    if os.path.isfile(newname):
                      #This checks for the presence of a new name file if it is there it updates the name with a number in the name
                      os.renames(newname, newname_add)
                    # break
                    else:
                      #if the original file name is found it is given a new name that is derived from the chemical identifier or name
                      os.renames(filename, newname)
                      break
       except FileNotFoundError:
            pass
       return

    def find_dlgs(self, dpf_files, dlg_files):
        for dpf_file, subdirs, files in os.walk(dpf_files):
              for filename_dlg in files:
                 if filename_dlg.endswith(".dlg"):
                     orifile = os.path.join(dpf_files, filename_dlg)
                     renamefile = os.path.join(dlg_files, filename_dlg)
                     shutil.move(orifile, renamefile)
                     #self.rename_file(renamefile)
                     del renamefile
        return

    def run_binana(self, binana_list, install_dir, dlg_files, protein_files):
        #this is needed to solve a problem with multiprocessing and not starting all jobs.
        #os.chdir(dlg_file)
        #with get_context("spawn").Pool() as sum_pool:
        #sum_pool = multiprocessing.Pool()
        work_list = []
        binana_loc = os.path.join(install_dir, "binana.py")
        for get_dlg_conf in binana_list:
           get_dlg_conf = get_dlg_conf.strip('\n')
           if get_dlg_conf.endswith(".dlg"):
              try:
                  with open(get_dlg_conf) as current_dlg_file:
                       info = current_dlg_file.read()
                       outcome_good = re.search("Successful", info)
                       if outcome_good:
                            prot_line = re.search("Macromolecule file used to create Grid Maps =	(\w+.+)", info)
                            lig_line = re.search("INPUT-LIGAND-PDBQT: REMARK  Name = ([\w-]+)", info)
                            for prot_path, prot_dir, prot_file in os.walk(protein_files):
                                for prot_file_found in prot_file:
                                    if prot_file_found.startswith(prot_line.group(1)) and lig_line:
                                        prot_name = prot_line.group(1)
                                        full_prot_path = os.path.join(prot_path, prot_name)
                                        first_part = str(binana_loc)+" -receptor "+str(full_prot_path)+" -ligand "+get_dlg_conf.replace(".dlg",".pdbqt")
                                        second_part = get_dlg_conf.replace('.dlg','_full_log.txt')
                                        work_dict = {}
                                        work_dict[first_part] = second_part
                                        work_list.append(work_dict)
                                        del work_dict
                                        #work_list.append(" "+str(binana_loc)+" -receptor "+str(full_prot_path)+" -ligand "+get_dlg_conf.replace(".dlg",".pdbqt")+" > "+get_dlg_conf.replace('.dlg','_full_log.txt'))
                                        current_dlg_file.close()
                                        #sum_pool.apply_async(self.run_summarize_interactions, args=(install_dir,full_prot_path,get_dlg_conf, ))
                                    else:
                                        pass
              except IOError as exc:
                  if exc.errno != errno.EISDIR: # Do not fail if a directory is found, just ignore it.
                     raise # Propagate other kinds of IOError.
                     break
                  else:
                     pass
        binana_list_file = os.path.join(working_project, "list_for_binana.json")
        with open(binana_list_file, "w") as file:
            json.dump(work_list, file)
        file.close()

    def run_summarize_docking(self, pythonsh_loc,summarize_dock_path,lig_key,prot_def):
        subprocess.run([str(pythonsh_loc),str(summarize_dock_path),"-l",str(lig_key),"-o","docking_summary_log.txt", "-a", "-b","-r",str(prot_def)])
        return

    def make_results_log(self, all_dlg_undone, dlg_files, protein_files, install_dir):
        if platform == "win32":
            print("windows was found")
            summarize_dock_path = os.path.join(install_dir,"mgltools_win32_1.5.6\\AutoDockTools\\Utilities24\\summarize_docking.py")
            pythonsh_loc = os.path.join(install_dir,"python.exe")
        elif platform == "linux" or "linux2":
            print("linux was found")
            summarize_dock_path = os.path.join(install_dir,"mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/summarize_docking.py")
            pythonsh_loc = os.path.join(install_dir,"mgltools_x86_64Linux2_1.5.6/bin/pythonsh")
        os.chdir(dlg_files)
        work_dict = {}
        if not os.path.isfile("docking_summary_log.txt"):
            docking_file = pd.DataFrame(columns=['lowestEnergy_dlgfn','#runs','#cl','#LEC','LE','rmsd_LE','#hb','#ESTAT'])
            docking_file.to_csv("docking_summary_log.txt", index=False)
        else:
            pass
        for dlg_file_found in all_dlg_undone:
            dlg_file_found = dlg_file_found.strip('\n')
            print("this is dlg_file_found: ", dlg_file_found)
            if dlg_file_found.endswith(".dlg"):
               #try:
                   heatmaps_left = open(dlg_file_found, "r").read()
                   with open(dlg_file_found) as current_dlg_file: # No need to specify 'r': this is the default.
                        info = current_dlg_file.read()
                        prot_line = re.search("Macromolecule file used to create Grid Maps =	(\w+.+)", info)
                        lig_line = re.search("INPUT-LIGAND-PDBQT: REMARK  Name = ([\w-]+)", info)
                        current_dlg_file.close()
                        for prot_path, prot_dir, prot_file in os.walk(protein_files):
                            for prot_file_found in prot_file:
                                if prot_file_found.startswith(prot_line.group(1)) and lig_line:
                                   prot_name = prot_line.group(1)
                                   full_prot_path = os.path.join(prot_path, prot_name)
                                   work_dict[dlg_file_found] = full_prot_path
                                   #self.run_summarize_docking(pythonsh_loc,summarize_dock_path,dlg_file_found,full_prot_path,)
                                else:
                                   pass
               #except IOError as exc:
            #     if exc.errno != errno.EISDIR: # Do not fail if a directory is found, just ignore it.
            #        raise # Propagate other kinds of IOError.
            #        break
            #     else:
            #        pass
        for lig_key, prot_def in work_dict.items():

            self.run_summarize_docking(pythonsh_loc,summarize_dock_path,lig_key,prot_def)
        return

    def gen_mol_prot_key(self, dlg_files):
            print("Starting Protein Key Info Add")
            docking_info = pd.read_csv("docking_summary_log.txt", delimiter=',', header=0, engine='python')
            molecule_name = []
            receptor_name = []
            dlgf_filename = []
            molecule_mass = []
            for dlg_file_name in docking_info['lowestEnergy_dlgfn']: # 'file' is a builtin type, 'name' is a less-ambiguous variable name.
                new_dlg_name = dlg_file_name+'.dlg'
                try:
                    with open(new_dlg_name) as dlg_file: # No need to specify 'r': this is the default.
                         info = dlg_file.read()
                         prot_line = re.search("Macromolecule file used to create Grid Maps =	(\w+.+)", info)
                         #BUG HERE if there is not title in the PDBQT this will be blank... therefore skiped in the heatmap
                         lig_line = re.search("INPUT-LIGAND-PDBQT: REMARK  Name = ([\w-]+)", info)
                         dlg_file.close()
                         if prot_line and lig_line:
                            molecule_name.append(lig_line.group(1))
                            receptor_name.append(prot_line.group(1).strip(".pdbqt"))
                            obmol = OB.OBMol()
                            obconversion = OB.OBConversion()
                            obconversion.SetInFormat("pdbqt")
                            molecule = obconversion.ReadFile(obmol, new_dlg_name.replace(".dlg",".pdbqt"))
                            molwght = obmol.GetExactMass()
                            molecule_mass.append(molwght)
                            dlgf_filename.append(new_dlg_name.strip(".dlg"))
                         else:
                            pass
                except IOError as exc:
                    if exc.errno != errno.EISDIR: # Do not fail if a directory is found, just ignore it.
                       raise # Propagate other kinds of IOError.
                       break
                    else:
                       pass
            docking_info['mol_name'] = molecule_name
            docking_info['protein_id'] = receptor_name
            docking_info['exact_mass'] = molecule_mass
            docking_info['Score_per_Mass'] = docking_info['Score_per_Mass'] =docking_info['LE']/docking_info['exact_mass']
            docking_info.to_csv("docking_summary_protein_key.txt")
            dlg_file.close()

            return

    def get_conformations(self, get_conf, pdbqt_name_base, pythonsh_loc, write_lig_path):
            pdbqt_name = pdbqt_name_base+".pdbqt"
            subprocess.run([str(pythonsh_loc),str(write_lig_path),"-f",get_conf,"-o",pdbqt_name])
            #os.popen(str(pythonsh_loc)+" "+str(write_lig_path)+" -f "+ get_conf+" -o "+pdbqt_name)
            return

    def convert_pdbqt_to_pdb(self, pdbqt_name_base):
        obconvert = OB.OBConversion()
        obconvert.SetInAndOutFormats("pdbqt", "pdb")
        obmol = OB.OBMol()
        conver_mol = pdbqt_name_base+".pdbqt"
        obconvert.ReadFile(obmol, conver_mol)
        obmol.AddHydrogens()
        obmol.CorrectForPH()
        obconvert.WriteFile(obmol, pdbqt_name_base+'.pdb')
        obconvert.CloseOutFile()
        return

    def write_le_pose(self, all_dlg_undone, dlg_files, install_dir):
        if platform == "win32":
           print("windows was found")
           write_lig_path = os.path.join(install_dir,"mgltools_win32_1.5.6\\AutoDockTools\\Utilities24\\write_lowest_energy_ligand.py")
           pythonsh_loc = os.path.join(install_dir,"python.exe")
        elif platform == "linux" or "linux2":
           print("linux was found")
           write_lig_path = os.path.join(install_dir,"mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/write_lowest_energy_ligand.py")
           pythonsh_loc = os.path.join(install_dir,"mgltools_x86_64Linux2_1.5.6/bin/pythonsh")
        for get_dlg_conf in all_dlg_undone:
            get_dlg_conf = get_dlg_conf.strip('\n')
            if get_dlg_conf.endswith(".dlg"):
               pdbqt_name_base = join(dlg_files, (get_dlg_conf.split('.')[0]))
               get_conf = join(dlg_files, get_dlg_conf)
               self.get_conformations(get_conf, pdbqt_name_base, pythonsh_loc, write_lig_path)
            ###################################################################################
            ####### RIGHT NOW THIS PART DOESN"T WORK. It may be useful to have the pdbs #######
            #for get_pdbqt_conf in files:
            #    if get_pdbqt_conf.endswith(".pdbqt"):
            #       print("this is the pdbqt file name ", get_dlg_conf)
            #       pdbqt_name_base = join(dlg_file, (get_pdbqt_conf.split('.')[0]))
            #       get_conf = join(dlg_file, get_pdbqt_conf)
            #       self.convert_pdbqt_to_pdb(pdbqt_name_base)
            #    return

    def make_parsed_heatmap(self, data_index_rows, data_frame, data_frame2, get_iactdata):
            print("making the heatmaps now")
            for i in list(data_index_rows):
                for j in list(data_frame):
                    data_to_split = data_frame.at[i, j]
                    if type(data_to_split) == str:
                       list_element_one = data_to_split.split(':', 2)
                       col_name = (list_element_one[0].replace(" ",""))
                       data_point = (list_element_one[1].replace(" ",""))
                    if col_name != '0' and i == "Short-range contacts":
                       #the color here should be unique
                       data_frame2.at[i, col_name] = 1
                    elif col_name != '0' and i == "Hydrogen bonds" :
                       #the color here should have a lot of 1 in it
                       data_frame2.at[i, col_name] = 2
                    elif col_name != '0' and i == "Hydrophobic contacts":
                       #the color here should be unique and disctinct from hbonds and short-range contacts
                       data_frame2.at[i, col_name] = 7
                    elif col_name != '0' and i == "pi-pi stacking interactions":
                       #this color should be closer to hydrophobic contacts than hbonds
                       data_frame2.at[i, col_name] = 5
                    elif col_name != '0' and i == "T-stacking":
                       #this color also should be closter to hydrophobic contacts than hbonds
                       data_frame2.at[i, col_name] = 6
                    elif col_name != '0' and i == "cation-pi":
                       #this color and the color for salt bridge should be closer to hbonds than the others
                       data_frame2.at[i, col_name] = 4
                    elif col_name != '0' and i == "Salt Bridges":
                       data_frame2.at[i, col_name] = 4
                    elif col_name == '0':
                       #print("setting the column name to 0")
                       pass
                    else:
                       data_frame2.at[i, col_name] = 0
            data_frame2.fillna(0, inplace=True)
            for i in np.arange(0, data_frame2.shape[0]):           #for lop to cast data to numeric
                data_frame2.iloc[i, :] = pd.to_numeric(data_frame2.iloc[i, :], errors= 'ignore')
            data_frame2.fillna(0, inplace=True)
            #these colors were picked to be oposites on the color scale so that hydrophobic type interactions would be in diiferent channels then electrostatic interactions such as hbonds
            colorscale = [[0, '#ffffff'], [.2, '#ff9900'], [.4, '#ff0000'], [.6, '#ff9966'], [.8, '#00ffff'], [1, '#0000ff']]#
            fig = go.Figure(go.Heatmap(z=data_frame2.values,x=data_frame2.columns,y=None,zmax=data_frame2.values.max(),zmin=0,showscale=False,colorscale=colorscale))
            fig.update_layout(
                font=dict(
                size=9,
                 ),
                autosize=False,
                width=356, height=356,
                margin=dict(
                    l=2,
                    r=2,
                    b=2,
                    t=2,
                    pad=0
            )
            )
            #py.offline.plot(fig, auto_open=False, output_type='file', filename=input_file.replace(".txt","_heatmap.html"))
            py.io.write_image(fig, get_iactdata.replace("_log.txt","_heatmap.png"), format='png', scale=1.0, width=356, height=356)
            #py.io.orca.shutdown_server()
            return

    def generate_sankey_diagram(self, data_index_rows, data_frame, data_frame2, get_iactdata, **kwargs):
            for i in list(data_index_rows):
                for j in list(data_frame):
                  data_to_split = data_frame.at[i, j]
                  if type(data_to_split) == str:
                     list_element_one = data_to_split.split(':', 2)
                     col_name = (list_element_one[0].replace(" ",""))
                     data_point = (list_element_one[1].replace(" ",""))
                  if col_name != '0':
                     data_frame2.at[i, col_name] = data_point
                  elif col_name == '0':
                     pass
                  else:
                     data_frame2.at[i, col_name] = 0
            for i in np.arange(0, data_frame2.shape[0]):           #for loop to cast data to numeric
                data_frame2.iloc[i, :] = pd.to_numeric(data_frame2.iloc[i, :], errors= 'ignore')
            data_frame2.fillna(0, inplace=True)

            source_list = []
            target_list = []
            value_list = []
            all_labels = []
            short_col_names = ""

            col_aa = []

            for i in list(data_frame2):
                col_name_aa_no = i.split("(", 2)
                col_aa.append(col_name_aa_no[1].replace(")",""))
                col_aa2 =  ', '.join(col_aa)

            print("putting info into Plotly Sankey generator")

            row_names = list(data_frame2.index)
            column_names = list(data_frame2)
            short_col_names = list(col_aa2)
            all_labels = row_names + column_names

            for i in list(data_frame2.index):
                for j in column_names:
                    target_list.append(column_names.index(j)+len(row_names))
                    source_list.append(row_names.index(i))
                    value_list.append(data_frame2.at[i,j])

            print("Generating Sankey diagrams")

            fig = go.Figure(data=[go.Sankey(
            node = dict(
              pad = 15,
              thickness = 50,
              line = dict(color = "black" , width = 0.7),
              label = list(all_labels),
              color = ["mediumblue", "mediumorchid", "mediumpurple", "mediumseagreen", "mediumslateblue", "mediumspringgreen", "mediumturquoise", "mediumvioletred", "mediumblue", "mediumorchid", "mediumpurple", "mediumseagreen", "mediumslateblue", "mediumspringgreen", "mediumturquoise", "mediumvioletred", "mediumblue", "mediumorchid", "mediumpurple", "mediumseagreen", "mediumslateblue", "mediumspringgreen", "mediumturquoise", "mediumvioletred", "mediumblue", "mediumvioletred"]
              ),
              link = dict(
                 source = list(source_list), # indices correspond to labels, eg A1, A2, A2, B1, ...
                 target = list(target_list),
                 value = list(value_list)
               #color = ["mediumblue", "mediumorchid", "mediumpurple", "mediumseagreen", "mediumslateblue", "mediumspringgreen", "mediumturquoise", "mediumvioletred"]
              ))])
            if 'sankey-html' == kwargs.get('param'):
                py.offline.plot(fig, auto_open=False, output_type='file', image_width=800, image_height=600, filename=get_iactdata.replace("_parsed_log.txt","_sankey.html"))
            if 'sankey-png' == kwargs.get('param'):
                py.io.write_image(fig, get_iactdata.replace("_parsed_log.txt","_sankey.png"), format='png', scale=None, width=800, height=600)
                #py.offline.plot(fig, auto_open=False, output_type='file', image_width=800, image_height=600, filename=get_iactdata.replace("_parsed_log.txt","_sankey.html"))
                #py.io.write_image(fig, get_iactdata.replace("_parsed_log.txt","_sankey.png"), format='png', scale=None, width=800, height=600)
                #py.io.orca.shutdown_server()
            print("this should be the end of this one... starting the next one")
            return

    def build_sankey_and_other_logs(self, file_list, **kwargs):
            result_list_tqdm = []
            for get_iactdata in file_list:
                get_iactdata = get_iactdata.strip('\n')
                if get_iactdata.endswith("parsed_log.txt"):
                   data_frame = pd.read_csv(get_iactdata, header=0, index_col=0, sep=', ', engine='python')
                   data_frame.fillna(0, inplace=True)
                   col_names = list(data_frame)
                   data_index_rows = data_frame.index.values
                   rows_list = list(data_index_rows)
                   data_frame2 = pd.DataFrame(data='0', index=rows_list, columns=col_names)
                   if kwargs.get('param') == 'heatmaps':
                       self.make_parsed_heatmap(data_index_rows, data_frame, data_frame2, get_iactdata)
                       continue
                   elif kwargs.get('param') == 'sankey-png':
                       self.generate_sankey_diagram(data_index_rows, data_frame, data_frame2, get_iactdata, param="sankey-png")
                       continue
                   elif kwargs.get('param') == 'sankey-html':
                       self.generate_sankey_diagram(data_index_rows, data_frame, data_frame2, get_iactdata, param="sankey-html")
                       continue
                   else:
                       break
            return

class Worker(QRunnable):
    '''
    Worker thread

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.

    :param callback: The function callback to run on this worker thread. Supplied args and
                     kwargs will be passed through to the runner.
    :type callback: function
    :param args: Arguments to pass to the callback function
    :param kwargs: Keywords to pass to the callback function

    '''

    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()

        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

    @pyqtSlot()
    def run(self):

        try:
            result = self.fn(*self.args, **self.kwargs)
            #started = pyqtSignal()
        except:
            exctype, value = sys.exc_info()[:2]
        else:
            self.signals.started.emit()
            self.signals.result.emit(result)  # Return the result of the processing
        finally:
            self.signals.finished.emit()  # Done


class Worker_docking(QRunnable):

    def __init__(self, fn, *args, **kwargs):
        super(Worker_docking, self).__init__()

        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()


    @pyqtSlot()
    def run(self):
        # Retrieve args/kwargs here; and fire processing using them
        try:
            result = self.fn(*self.args, **self.kwargs)
        except:
            #traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            #self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.signals.result.emit(result)  # Return the result of the processing
        finally:
            self.signals.finished.emit()  # Done

class MainWindow(QMainWindow):
    message_status = pyqtSignal(str)
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.title = 'Molecular docking and deep-learning toxin prediction'
        self.left = 3
        self.top = 28
        self.width = 1300
        self.height = 640
        self.initUI()
        print("\n Molecular docking and deep-learning toxin prediction \n GUI was written by Michael J. McCarthy, 2019")

    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        self.status = QStatusBar()
        self.setStatusBar(self.status)

        self.createTable()
        self.createtopLeftGroupBox()
        self.createtopRightGroupBox()

        ##Using docked widg
        localleft_groupbox = QGroupBox()
        layoutleft = QVBoxLayout()
        layoutleft.addLayout(self.newlayoutleft)
        localleft_groupbox.setLayout(layoutleft)
        dockWidget_setupdocking = QDockWidget('Setup Docking', self)
        dockWidget_setupdocking.setWidget(localleft_groupbox)
        dockWidget_setupdocking.setFloating(False)

        localright_groupbox = QGroupBox()
        layoutright = QVBoxLayout()
        layoutright.addLayout(self.newlayout)
        localright_groupbox.setLayout(layoutright)
        dockWidget_dockingresults = QDockWidget('Docking Results', self)
        dockWidget_dockingresults.setWidget(localright_groupbox)
        dockWidget_dockingresults.setFloating(False)

        #self.editor = QTextEdit()

        self.setCentralWidget(self.tableWidget)
        self.addDockWidget(Qt.LeftDockWidgetArea, dockWidget_setupdocking)
        self.addDockWidget(Qt.TopDockWidgetArea, dockWidget_dockingresults)

        file_menu = self.menuBar().addMenu("&File")

        open_file_action = QAction("Open molecules info file", self)
        open_file_action.setStatusTip("Open file")
        open_file_action.triggered.connect(self.file_open)
        file_menu.addAction(open_file_action)

        start_molfile_action = QAction("Input molecules from table", self)
        start_molfile_action.setStatusTip("Input molecules from table, by name or in SMILE or InChI format")
        start_molfile_action.triggered.connect(self.molfile_start)
        file_menu.addAction(start_molfile_action)

        open_molfile_action = QAction("Input molecules by 3D SDF file", self)
        open_molfile_action.setStatusTip("Open a 2D or 3D molecules file")
        open_molfile_action.triggered.connect(self.molfile_open)
        file_menu.addAction(open_molfile_action)

        file_menu.addSeparator()

        ##This next section is action for the file menu
        add_column = QAction("Add Column", self)
        add_column.setStatusTip("Add a new column to the table")
        add_column.triggered.connect(self.addcolumn)
        file_menu.addAction(add_column)

        add_row = QAction("Add Row", self)
        add_row.setStatusTip("Add a new row to the table")
        add_row.triggered.connect(self.addrow)
        file_menu.addAction(add_row)

        add_project_setup = QAction("1. Add Docking Project", self)
        add_project_setup.setStatusTip("Create a new docking project")
        add_project_setup.triggered.connect(self.setup_docking_project)
        file_menu.addAction(add_project_setup)

        add_setup_docking = QAction("2. Setup Docking Project", self)
        add_setup_docking.setStatusTip("Generates all files for docking with AutoDock")
        add_setup_docking.triggered.connect(self.setup_docking_files)
        file_menu.addAction(add_setup_docking)

        add_run_docking = QAction("3. Run Docking", self)
        add_run_docking.setStatusTip("This will run docking with AutoDock")
        add_run_docking.triggered.connect(self.docking_button_clicked_run)
        file_menu.addAction(add_run_docking)

        add_process_docking = QAction("4. Process Docking", self)
        add_process_docking.setStatusTip("This will processs the docking logs and generate summary of binding scores")
        add_process_docking.triggered.connect(self.process_results)
        file_menu.addAction(add_process_docking)

        file_menu.addSeparator()

        add_existing_docking = QAction("First Check Restart Project", self)
        add_existing_docking.setStatusTip("This will check for where you left off with a project")
        add_existing_docking.triggered.connect(self.check_project_status)
        file_menu.addAction(add_existing_docking)

        add_process_docking = QAction("Second Restart Project", self)
        add_process_docking.setStatusTip("Pick up where you left off with a project")
        add_process_docking.triggered.connect(self.restart_process_results)
        file_menu.addAction(add_process_docking)

        file_menu.addSeparator()

        add_existing_docking = QAction("Open Existing Project", self)
        add_existing_docking.setStatusTip("Open an existing project")
        add_existing_docking.triggered.connect(self.open_exising_project)
        file_menu.addAction(add_existing_docking)

        add_project_remove = QAction("Delete Docking Project", self)
        add_project_remove.setStatusTip("Delete a docking project... This can not be undone!")
        add_project_remove.triggered.connect(self.remove_docking_project)
        file_menu.addAction(add_project_remove)

        file_menu.addSeparator()

        add_gen3dfrom2d = QAction("Generate 3D molecules", self)
        add_gen3dfrom2d.setStatusTip("This will attempt to generate 3D molecules from 2D input such as SMILES or InChI. It will look for a text file (\".txt\") in the molecule_files project direcorty")
        add_gen3dfrom2d.triggered.connect(self.gen3dfrom2d)
        file_menu.addAction(add_gen3dfrom2d)

        file_menu.addSeparator()

        exit_file_action = QAction("Exit", self)
        exit_file_action.setStatusTip("Exit Application")
        exit_file_action.triggered.connect(self.file_exit)
        file_menu.addAction(exit_file_action)
    def molfile_open(self):
        try:
            working_project
            fileName, _ = QFileDialog.getOpenFileName(self,'Open file', None,"All Files (*)")
            if fileName == '':
               return
            if fileName != '' and QFileInfo(fileName).suffix() == "":
               fileName += '.txt'
               shutil.copy(fileName, molecule_files)
               os.chdir(working_project)
               file = os.listdir(molecule_files)
               send_sdf_file = os.path.join(molecule_files, file.pop(0))
               generate_3Dligands_10_28_2020.update_sdf_file(sent_sdf_file)
            else:
               shutil.copy(fileName, molecule_files)
               os.chdir(working_project)
               file = os.listdir(molecule_files)
               sent_sdf_file = os.path.join(molecule_files, file.pop(0))
               generate_3Dligands_10_28_2020.update_sdf_file(sent_sdf_file)
               return
        except NameError:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText("Open or start a project")
            #msg.setInformativeText("To import molecule there needs to be a project directory with a molecules directory. These will be generated if a new project is started or found if an existing project is opened")
            msg.setWindowTitle("Working directory not found")
            msg.setDetailedText("To import molecules there needs to be a project directory with a molecule_files directory. This will be generated if a new project is started or found if an existing project is opened")
            msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
            msg.exec_()
    def file_open(self):#There is a problem that I need to fix. if there is an empty line pasteing will only work for the stuff after that line and if there is an empty line at the end it will only put in empy lines
        try:
            os.chdir(working_project)
            fileName, _ = QFileDialog.getOpenFileName(self,'Open molecules info file', None,"All Files (*)")
            if fileName == '':
               return
            if fileName != '':
               if QFileInfo(fileName).suffix() == "": fileName += '.txt'
            with open(fileName, 'r') as stream:
                file_data = csv.reader(stream, delimiter=':')
                header = next(file_data, None)
                self.tableWidget.setRowCount(0)
                self.tableWidget.setColumnCount(0)
                for item in header:
                    row = self.tableWidget.columnCount()
                    self.tableWidget.insertColumn(row)
                    self.tableWidget.setRowCount(len(item))
                    head = QTableWidgetItem(item)
                    self.tableWidget.setHorizontalHeaderItem(row, head)
                self.tableWidget.setRowCount(0)
                for rowdata in csv.reader(stream, delimiter=':'):
                     row = self.tableWidget.rowCount()
                     self.tableWidget.insertRow(row)
                     self.tableWidget.setColumnCount(len(rowdata))
                     for column, data in enumerate(rowdata):
                        if data != '':
                           png_file = re.search("([\w-]+).png", data)
                           if png_file:
                               pic  = png_file.group(0)
                               pic_path = os.path.join(os.getcwd(), pic)
                               picture = QPixmap(pic_path)
                               self.label = QLabel()
                               self.label.setPixmap(picture)
                               self.tableWidget.setCellWidget(row, column, self.label)
                           item = QTableWidgetItem(data)
                           self.tableWidget.setItem(row, column, item)
                        else:
                             continue
        except NameError:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText("Open or start a project")
            #msg.setInformativeText("To import molecule there needs to be a project directory with a molecules directory. These will be generated if a new project is started or found if an existing project is opened")
            msg.setWindowTitle("Working directory not found")
            msg.setDetailedText("To Open a molecules file there needs to be a project directory with a molecule_files directory. This will be generated if a new project is started or found if an existing project is opened")
            msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
            msg.exec_()
    def file_save(self):
        fileName, _ = QFileDialog.getSaveFileName(self,'Save file', None,"All Files (*)")
        if fileName == '':
           return
        if fileName != '':
           if QFileInfo(fileName).suffix() == "": fileName += '.txt'
        with open(fileName, 'w') as stream:
             writer = csv.writer(stream)
             for row in range(self.tableWidget.rowCount()):
                 rowdata = []
                 for column in range(self.tableWidget.columnCount()):
                     item = self.tableWidget.item(row, column)
                     if item:
                       rowdata.append(item.text())
                     else:
                       rowdata.append('')
                 writer.writerow(rowdata)
        shutil.copy(fileName, molecule_files)
    def molfile_start(self):
           try:
               working_project
           except NameError:
               msg = QMessageBox()
               msg.setIcon(QMessageBox.Information)
               msg.setText("Open or start a project")
               msg.setWindowTitle("Working directory not found")
               msg.setDetailedText("To import molecules there needs to be a project directory with a molecule_files directory. This will be generated if a new project is started or found if an existing project is opened")
               msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
               msg.exec_()
           else:
               os.chdir(working_project)
               fileName, _ = QFileDialog.getSaveFileName(self,'Save file', working_project,"All Files (*)")
               if fileName == '':
                  return
               if fileName != '':
                  if QFileInfo(fileName).suffix() == "": fileName += '_from_table.txt'
               with open(fileName, 'w') as stream:
                    writer = csv.writer(stream)
                    headers = []
                    for column in range(self.tableWidget.columnCount()):
                        header = self.tableWidget.horizontalHeaderItem(column)
                        if header is not None:
                           headers.append(header.text())
                        else:
                           continue
                    writer.writerow(headers)
                    for row in range(self.tableWidget.rowCount()):
                        rowdata = []
                        for column in range(self.tableWidget.columnCount()):
                            item = self.tableWidget.item(row, column)
                            if item:
                               rowdata.append(item.text())
                            else:
                               continue
                        writer.writerow(rowdata)
               shutil.copy(fileName, molecule_files)
               os.chdir(molecule_files)
               self.threadpoolmolstart = QThreadPool()
               self.start_mols = process_input()
               def interact_diagrams(fileName):
                   self.start_mols.start_the_search(fileName)
               worker_processmolstart = Worker(interact_diagrams, fileName)
               self.threadpoolmolstart.start(worker_processmolstart)
               self.start_mols.message_status.connect(self.update_message)
               worker_processmolstart.signals.finished.connect(self.thread_complete)
               self.threadpoolmolstart.clear()
    def file_exit(self):
        QApplication.quit()
    def setup_docking_project(self):
        ## the only purpose of this function is to make a tar ball of an existing directory
        def tardir(path, tar_name):
            with tarfile.open(tar_name, "w:gz") as tar_handle:
                 for root, dirs, files in os.walk(path):
                     for file in files:
                         tar_handle.add(os.path.join(root, file))
        ## This will open a dialog box to get the name of a new project and then a directory will be made with that name with all project directoies
        text_input, okPressed = QInputDialog.getText(self, "Project Name","Enter a New Project Name:", QLineEdit.Normal, "")
        ## Set the location of the install directory because it should contain helper software such as mgltools, autodock, binana and etc
        global install_dir
        install_dir = os.getcwd()
        ## this will put a date stamp on a existing directory for archiving if the user wants to use the same project name
        Current_Date = datetime.datetime.today().strftime ('%d-%b-%Y')
        home = os.getcwd()
        ## this will provide information of current directories to check for re-used names later
        check_existing = os.listdir(home)
        if okPressed and text_input != '':
            ## now we check for existing directories
            if text_input in check_existing:
                text = text_input
                ## if an existing project name is found ask the user again for a new name or if they want to overwite
                text_new, okPressed = QInputDialog.getText(self, "Project Already Used","Enter a New Name. \n If You Enter The Same Name Again \n The Old Project Will Be Overwritten", QLineEdit.Normal, "")
                ## if they want to overwirite we are actually goind to archive just-in-case
                if okPressed and text_new != '' and text_new in check_existing:
                    tar_name = text_new+str(Current_Date)
                    ## renanme the existing project with the current date
                    os.rename(text_new, tar_name)
                    path = home+"/"+tar_name
                    ## zip it up and then delete the renamed project so it only exists as an archive
                    tardir(path, tar_name+".tar.gz")
                    shutil.rmtree(tar_name)
                    text = text_input
            else:
                text = text_input
                pass
            working_directories = ["dpf_files", "molecule_files", "protein_files", "dlg_files", "results_files"]
            project_path = os.getcwd()
            try:
                os.mkdir(project_path+"/"+text)
                os.chdir(project_path+"/"+text)
                new_project_path = os.getcwd()
            except:
                return
            for directory in working_directories:
                pathvar = directory
                if directory == "protein_files":
                   pass
                else:
                   os.mkdir(new_project_path+"/"+directory)
                ## the next two lines makes a project source file which is a list of directories and their location
                ## tdpa is for ToxDockingProjectAutodock so for vina this would be tdpv
                if directory == "protein_files":
                   project_source_file = pathvar+" = "+"\""+project_path+"/"+directory+"\""
                else:
                   project_source_file = pathvar+" = "+"\""+new_project_path+"/"+directory+"\""
            os.chdir(project_path)
            # try to change these to private variables and then pass them into setup_docking_project function
            #shutil.copy("multithread_example.py", new_project_path)
            shutil.copy("parallel_autodock4.py", new_project_path)
            shutil.copy("collect_docking_interaction_info2.py", new_project_path)
            if platform == 'win32':
                shutil.copy("autodock4.exe", new_project_path)
            else:
                pass
            global dpf_files
            dpf_files = os.path.join(new_project_path,"dpf_files")
            global molecule_files
            molecule_files = os.path.join(new_project_path,"molecule_files")
            global protein_files
            protein_files = os.path.join(project_path,"protein_files")
            global dlg_files
            dlg_files = os.path.join(new_project_path,"dlg_files")
            global results_files
            results_files = os.path.join(new_project_path,"results_files")
            global working_project
            working_project = os.path.join(project_path, text)
        else:
            return
        self.label_loc.setText("Current working directory: "+str(working_project))
        self.tree.setRootIndex(self.model.index(working_project))
        return working_project

    def gen3dfrom2d(self):
        ## this will branch based on user input. For openbabel we need info on how to run simulations
        ## But rdit just uses the internal EDGK method
         try:
            molecule_files
            def getInteger( method, ffield, mol_file):
                steps, okPressed = QInputDialog.getInt(self, "Number of steps","Steps:", 500, 0, 1000, 1)
                if okPressed:
                   print(steps, method, ffield)
                if method == 'Openbabel':
                   fileName = os.path.join(molecule_files,mol_file)
                   generate_3Dligands_10_28_2020.gen_obabel_mols(fileName,ffield,steps)

            def getDouble():
                d, okPressed = QInputDialog.getDouble(self, "Get double","Value:", 10.50, 0, 100, 10)
                if okPressed:
                   print( d)

            def getChoice_forcefield(method,mol_file):
                items = ("MMFF94", "Ghemical", "GAFF", "MMFF94s", "UFF")
                ffield, okPressed = QInputDialog.getItem(self, "Choose RdKit or Openbabel to generate 3D structures","Method:", items, 0, False)
                if okPressed and ffield:
                   getInteger(method, ffield, mol_file)

            def getChoice(mol_file):
                items = ("RdKit", "Openbabel")
                method, okPressed = QInputDialog.getItem(self, "Choose Openbabel forcefield to generate 3D structures","Method:", items, 0, False)
                if okPressed and method == "Openbabel":
                   getChoice_forcefield(method,mol_file)
                elif okPressed and method == "RdKit":
                   smi_input_file = os.path.join(molecule_files,mol_file)
                   generate_3Dligands_10_28_2020.gen_rdkit_mols(smi_input_file)
                else:
                   pass

            #getChoice()
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText("Only use if no 3D structures are available")
            msg.setInformativeText("For further instructs on InChI and SMILES input read message")
            msg.setWindowTitle("Generate 3D structures")
            msg.setDetailedText("Warning: These methods are simple and may make mistakes, only use if you can't find 3D structures and can't generate them any other way. For InChI input put quotes around entire InChi stirng. Add SMILES and InChI, one per line and one column only with no header row")
            msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
            choice = msg.exec()
            if choice == QMessageBox.Cancel:
                pass
            elif choice == QMessageBox.Ok:
                for mol_file in os.listdir(molecule_files):
                   if mol_file.endswith(".txt"):
                      getChoice(mol_file)
                      for info_file in os.listdir(working_project):
                          if info_file.startswith('excel_molecules'):
                              with open(info_file, 'r') as stream:
                                  file_data = csv.reader(stream, delimiter=':')
                                  header = next(file_data, None)
                                  self.tableWidget.setRowCount(0)
                                  self.tableWidget.setColumnCount(0)
                                  for item in header:
                                      row = self.tableWidget.columnCount()
                                      self.tableWidget.insertColumn(row)
                                      self.tableWidget.setRowCount(len(item))
                                      head = QTableWidgetItem(item)
                                      self.tableWidget.setHorizontalHeaderItem(row, head)
                                  self.tableWidget.setRowCount(0)
                                  for rowdata in csv.reader(stream, delimiter=':'):
                                          row = self.tableWidget.rowCount()
                                          self.tableWidget.insertRow(row)
                                          self.tableWidget.setColumnCount(len(rowdata))
                                          for column, data in enumerate(rowdata):
                                              if data != '':
                                                 png_file = re.search("mol_([\w-]+).png", data)
                                                 if png_file:
                                                    pic  = png_file.group(0)
                                                    pic_path = os.path.join(os.getcwd(), pic)
                                                    picture = QPixmap(pic_path)
                                                    self.label = QLabel()
                                                    self.label.setPixmap(picture)
                                                    self.tableWidget.setCellWidget(row, column, self.label)
                                                 item = QTableWidgetItem(data)
                                                 self.tableWidget.setItem(row, column, item)
                                              else:
                                                    continue
         except NameError:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText("Open or start a project")
            msg.setWindowTitle("Working directory not found")
            msg.setDetailedText("To generate 3D molecules there needs to be a project directory with a molecule_files directory and a file containing a list of SMILES or InChI, one per line and only one column. This will be generated if a new project is started or found if an existing project is opened and there is a text file with the appropriate information")
            msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
            msg.exec_()

    def open_exising_project(self):
        ## Fist I am going to put launch path as install_dir I will probably have to change this later
        global install_dir
        try:
            install_dir
        except:
            install_dir = os.getcwd()
        ## This will open a dialog box to get the name of a new project and then a directory will be made with that name with all project directoies
        filename_expro = QFileDialog.getExistingDirectory(self,"Existing Project")
        #home = os.getcwd()
        ## this will provide information of current directories to check for re-used names later
        #check_existing = os.listdir(home)
        if filename_expro != '':
            ## now we check for existing directories
            os.chdir(filename_expro)
            working_directories = ["dpf_files", "molecule_files", "protein_files", "dlg_files", "results_files"]
            project_path = os.getcwd()
            try:
                new_project_path = os.getcwd()
            except:
                return
            # try to change these to private variables and then pass them into setup_docking_project function
            global dpf_files
            dpf_files = os.path.join(new_project_path,"dpf_files")
            global molecule_files
            molecule_files = os.path.join(new_project_path,"molecule_files")
            global protein_files
            protein_files = os.path.join(install_dir,"protein_files")
            global dlg_files
            dlg_files = os.path.join(new_project_path,"dlg_files")
            global results_files
            results_files = os.path.join(new_project_path,"results_files")
            global working_project
            working_project = new_project_path
            global sankey_dir
            sankey_dir = os.path.join(new_project_path,"dlg_files")
        else:
            return
        if self.formLayout.itemAt(0,1) != None:
            while self.formLayout.itemAt(0,1) != None:
                  self.formLayout.removeRow(0)
        else:
            pass
        if self.listwidget.item(0) != None:
            while self.listwidget.item(0) != None:
                #thing = self.listwidget.row(0)
                self.listwidget.clear()
        else:
            pass
        if self.listwidget2.item(0) != None:
            while self.listwidget2.item(0) != None:
                #thing = self.listwidget2.row(0)
                self.listwidget2.clear()
        else:
            pass
        if self.listwidget3.item(0) != None:
            while self.listwidget3.item(0) != None:
                #thing = self.listwidget3.row(0)
                self.listwidget3.clear()
        else:
            pass
        self.label_loc.setText("Current working directory: "+str(working_project))
        self.tree.setRootIndex(self.model.index(working_project))
        return working_project

    def setup_docking_files(self):
        ##This should utilize a threads to set up molecules for docking.
        try:
            os.chdir(working_project)
            def step6():
                ##This final step makes a docking list that will be used later for docking or could be moved to a cluster along with parallel_autodock4.py and the dpf_files dir and then run separately witha pbs script
                self.threadpoolls = QThreadPool()
                self.setup_docking = setup_Autodocking.setup_docking()
                def make_list(dpf_files, working_project):
                    self.setup_docking.make_list_to_dock(dpf_files, working_project)
                worker_processE = Worker(make_list, dpf_files, working_project)
                self.message_prog6 = process_setup_progress(param="writing_docking_list")
                self.threadpoolls.start(worker_processE)
                self.message_prog6.start()
                self.message_prog6.message_status.connect(self.update_message)
                worker_processE.signals.finished.connect(self.thread_complete)
                self.threadpoolsep.clear()
            def step5():
                ##Initially the dpf files are written in the proteins_files directory so they need to be moved to the dpf_files directory in the current project
                self.threadpoolmv = QThreadPool()
                self.setup_docking = setup_Autodocking.setup_docking()
                def move_dpf(protein_files, molecule_files, dpf_files):
                    self.setup_docking.move_files_for_docking(protein_files, molecule_files, dpf_files)
                worker_processD = Worker(move_dpf, protein_files, molecule_files, dpf_files)
                self.message_prog5 = process_setup_progress(param="moving_files")
                self.threadpoolmv.start(worker_processD)
                self.message_prog5.start()
                self.message_prog5.message_status.connect(self.update_message)
                worker_processD.signals.finished.connect(step6)
            def step4():
                ## this step will write the dpf files
                self.threadpooldpf = QThreadPool()
                self.setup_docking = setup_Autodocking.setup_docking()
                def setup_dpf(molecule_files, protein_files):
                    self.setup_docking.setup_dpf_files(molecule_files, protein_files)
                worker_processC = Worker(setup_dpf, molecule_files, protein_files)
                self.message_prog4 = process_setup_progress(param="dpf_files")
                self.threadpooldpf.start(worker_processC)
                self.message_prog4.start()
                self.message_prog4.message_status.connect(self.update_message)
                worker_processC.signals.finished.connect(step5)
                self.threadpooldpf.clear()
                return
            def step3():
                ## This step will use prepare_ligand4.py from AutoDockTools (within mgltools) to prepare pdbqt files for each molecule
                ## However, MACOS can't use mgltools anymore so we must use openbabel... results may not be the same when comparing pdbqt files
                self.threadpoolgenqt = QThreadPool()
                self.setup_docking = setup_Autodocking.setup_docking()
                def make_pdbqt(molecule_files,install_dir):
                    print("starting step3 generate pdbqt")
                    self.setup_docking.generate_pdbqt(molecule_files, install_dir)
                worker_processB2 = Worker(make_pdbqt, molecule_files, install_dir)
                self.message_prog3 = process_setup_progress(param="gen_pdbqt")
                self.message_prog3.start()
                self.threadpoolgenqt.start(worker_processB2)
                self.message_prog3.message_status.connect(self.update_message)
                worker_processB2.signals.finished.connect(step4)
                self.threadpoolgenqt.clear()

                return
            def step2():
                ## We must covert molecules to pdb format to use in step 3
                self.threadpoolgen = QThreadPool()
                self.setup_docking = setup_Autodocking.setup_docking()
                def make_pdb(molecule_files,install_dir):
                    print("starting step2 generate pdb")
                    self.setup_docking.generate_pdb(molecule_files, install_dir)
                worker_processB = Worker(make_pdb, molecule_files, install_dir)
                self.message_prog2 = process_setup_progress(param="gen_pdb")
                self.threadpoolgen.start(worker_processB)
                self.message_prog2.start()
                self.message_prog2.message_status.connect(self.update_message)
                self.threadpoolgen.clear()
                if platform != "darwin":
                   worker_processB.signals.finished.connect(step3)
                   self.threadpoolgen.clear()
                else:
                   worker_processB.signals.finished.connect(step4)
                   self.threadpoolgen.clear()
                return
            def step1():
                ## we need to work with individual molecule files so this will split up master files
                self.threadpoolsep = QThreadPool()
                self.setup_docking = setup_Autodocking.setup_docking()
                def separate_mols(molecule_files, install_dir):
                    print("starting separate_mols in setup")
                    self.setup_docking.setup_mol_file_for_docking(molecule_files, install_dir)
                worker_processA = Worker(separate_mols, molecule_files, install_dir)
                self.message_prog = process_setup_progress(param="separate_mols")
                self.threadpoolsep.start(worker_processA)
                self.message_prog.start()
                self.message_prog.message_status.connect(self.update_message)
                worker_processA.signals.finished.connect(step2)
                self.threadpoolsep.clear()
                return
            #return
            step1()
        except:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText("Open or start a project")
            msg.setWindowTitle("Working directory not found")
            msg.setDetailedText("To setup a docking project there needs to be a project directory with a molecule_files directory and a file containing a list of SMILES or InChI, one per line and only one column. This will be generated if a new project is started or found if an existing project is opened and there is a text file with the appropriate information")
            msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
            msg.exec_()

    def execute_this_fn(self):  #I did try subprocess.run but it really slows this down. best results with os.system
        ##This won't be recgonized on windows I need to make the platform dependent
        if not platform == "win32":
            version_check = subprocess.check_output("mpirun --version", shell=True)
            if "(Open MPI) 4" in str(version_check):
                if platform == 'darwin':
                    os.system("mpirun --use-hwthread-cpus -n "+str(self.text)+" python parallel_autodock4.py")
                elif platform == "linux" or "linux2":
                    os.system("mpirun --use-hwthread-cpus -n "+str(self.text)+" python parallel_autodock4.py") #subprocess.run(["mpirun","--use-hwthread-cpus",str(self.text),"python parallel_autodock4.py",">>","list_of_docked.txt","2>&1"])
            else:
                if platform == 'darwin':
                    os.system("mpirun -n "+str(self.text)+" python parallel_autodock4.py")
                elif platform == "linux" or "linux2":
                    os.system("mpirun  -n "+str(self.text)+" python parallel_autodock4.py")

        elif platform == "win32":
             print("starting on windows")
             os.system("mpiexec -n "+str(self.text)+" python parallel_autodock4.py")

    def execute_binana(self):
        if not platform == "win32":
            version_check = subprocess.check_output("mpirun --version", shell=True)
            if "(Open MPI) 4" in str(version_check):
                if platform == 'darwin':
                    os.system("mpirun --use-hwthread-cpus -n "+str(self.text)+" python collect_docking_interaction_info2.py")
                elif platform == "linux" or "linux2":
                    os.system("mpirun --use-hwthread-cpus -n "+str(self.text)+" python collect_docking_interaction_info2.py")
            else:
                if platform == 'darwin':
                    os.system("mpirun -n "+str(self.text)+" python collect_docking_interaction_info2.py")
                elif platform == "linux" or "linux2":
                    os.system("mpirun -n "+str(self.text)+" python collect_docking_interaction_info2.py")
        elif platform == "win32":
             os.system("mpiexec -n "+str(self.text)+" python collect_docking_interaction_info2.py")

    def thread_start(self):
        self.label_update.setText("Job Started")
        return

    def thread_complete(self):
        print("job done")
        self.label_update.setText("Job Done")
        return

    def createProgressBar(self):
        self.progressBar = QProgressBar()
        self.progressBar.setRange(0, 100)
        return

    def update_status(self, prog):
        self.progressBar.setValue(prog)
        return

    def update_message(self, message):
        self.label_update.setText(message)
        return

    def docking_button_clicked_run(self):
        try:
            os.chdir(working_project)
            self.threadpool = QThreadPool()
            maxnum = self.threadpool.maxThreadCount()
            self.text, okPressed = QInputDialog.getText(self, "Choose number of processor to use","Number of processors to use \n or click OK to use the maximum \n Max detected: "+str(maxnum)+ " :", QLineEdit.Normal, "")
            if okPressed and self.text == "":
              self.text = maxnum - 1
            elif okPressed and int(self.text) != maxnum and not int(self.text) == 0:
              print("number other than 0 or max: ",self.text)
            elif okPressed and int(self.text) == '0':
              self.text = 1
            elif not okPressed:
              return
            else:
              self.text = maxnum - 1
            self.label_update.setText("Starting")
            worker = Worker_docking(self.execute_this_fn) # Any other args, kwargs are passed to the run function
            self.prog_check = check_progres()
            self.prog_check.prog_status.connect(self.update_status)
            self.message_prog = check_progres()
            self.message_prog.message_status.connect(self.update_message)
            # Execute
            self.threadpool.start(worker)
            self.prog_check.start()
            self.message_prog.start()
            worker.signals.finished.connect(self.thread_complete)
        except:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText("Open or start a project")
            msg.setWindowTitle("Working directory not found")
            msg.setDetailedText("To run docking there needs to be a project directory with a molecule_files, dpf files etc. This will be generated if a new project is started or found if an existing project is opened and there is a text file with the appropriate information")
            msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
            msg.exec_()
        return

    def process_results(self):
        try:
            os.chdir(dlg_files)
            def summarize_and_plot():
                os.chdir(results_files)
                self.threadpoolsum = QThreadPool()
                self.threadpoolplot = QThreadPool()
                self.start_log = process_results()
                scores_info_file = "docking_summary_protein_key.txt"
                prediction_info_file = "table_of_prediction.csv"
                def plot_final(scores_info_file, prediction_info_file):
                    self.start_log.plot_merg_scores_and_prediction(scores_info_file, prediction_info_file)
                def summarize_interact(dlg_files, protein_files):
                    self.start_log.sumarize_results_and_generate_tables(dlg_files, protein_files)
                worker_processA = Worker(summarize_interact, dlg_files, protein_files)
                worker_processB = Worker(plot_final, scores_info_file, prediction_info_file)
                self.prog_check = process_results_progress(param="final_work")
                self.message_prog = process_results_progress(param="final_work")
                self.prog_check.start()
                self.message_prog.start()
                self.threadpoolsum.start(worker_processA)
                self.threadpoolplot.start(worker_processB)
                self.prog_check.prog_status.connect(self.update_status)
                self.message_prog.message_status.connect(self.update_message)
                worker_processB.signals.finished.connect(self.thread_complete)
                self.threadpoolplot.clear()
                worker_processA.signals.finished.connect(self.thread_complete)
                self.threadpoolsum.clear()
                return
            def run_model():
                self.threadpoolbina = QThreadPool()
                self.start_log = process_results()
                def interact_diagrams(install_dir, results_files, dlg_files):
                    self.start_log.start_predictions(install_dir, results_files, dlg_files)
             #       self.start_log.run_binana(install_dir, dlg_files, protein_files)
                worker_process2 = Worker(interact_diagrams, install_dir, results_files, dlg_files)
                self.prog_check = process_results_progress(param="model_run")
                self.message_prog = process_results_progress(param="model_run")
                self.prog_check.start()
                self.message_prog.start()
                self.threadpoolbina.start(worker_process2)
                self.prog_check.prog_status.connect(self.update_status)
                self.message_prog.message_status.connect(self.update_message)
                worker_process2.signals.finished.connect(self.thread_complete)
                worker_process2.signals.finished.connect(summarize_and_plot)
                self.threadpoolbina.clear()
                return
            def model_input():
                self.threadpoolbina = QThreadPool()
                self.start_log = process_results()
                def deep_learning_input(dlg_files, results_files):
                     self.start_log.generate_deep_learning_input(dlg_files, results_files)
                worker_process2 = Worker(deep_learning_input, dlg_files, results_files)
                self.prog_check = process_results_progress(param="find_results")
                self.message_prog = process_results_progress(param="find_results")
                self.prog_check.start()
                self.message_prog.start()
                self.threadpoolbina.start(worker_process2)
                self.prog_check.prog_status.connect(self.update_status)
                self.message_prog.message_status.connect(self.update_message)
                worker_process2.signals.finished.connect(self.thread_complete)
                #time.sleep(2)
                worker_process2.signals.finished.connect(run_model)
                self.threadpoolbina.clear()
                return
            def sankey_html():
                list_from_scratch = []
                for i in os.listdir(dlg_files):
                    if i.endswith("parsed_log.txt"):
                       list_from_scratch.append(str(dlg_files)+"/"+i)
                self.threadpoolprog = QThreadPool()
                self.start_log = process_results()
                def interact_diagrams(working_project):
                    self.start_log.build_sankey_and_other_logs(list_from_scratch, param="sankey-html")
                worker_process2 = Worker(interact_diagrams, working_project)
                self.prog_check = process_results_progress(param="sankey-html")
                self.message_prog = process_results_progress(param="sankey-html")
                self.prog_check.start()
                self.message_prog.start()
                self.threadpoolprog.start(worker_process2)
                self.prog_check.prog_status.connect(self.update_status)
                self.message_prog.message_status.connect(self.update_message)
                worker_process2.signals.finished.connect(self.thread_complete)
                worker_process2.signals.finished.connect(model_input)
                self.threadpoolprog.clear()
                return
            def sankey_png():
                list_from_scratch = []
                for i in os.listdir(dlg_files):
                    if i.endswith("parsed_log.txt"):
                       list_from_scratch.append(str(dlg_files)+"/"+i)
                self.threadpoolprog = QThreadPool()
                self.start_log = process_results()
                def start_png_list(working_project):
                    self.start_log.build_sankey_and_other_logs(list_from_scratch, param="sankey-png")
                worker_process2 = Worker(start_png_list, working_project)
                self.prog_check = process_results_progress(param="sankey-png")
                self.message_prog = process_results_progress(param="sankey-png")
                self.prog_check.start()
                self.message_prog.start()
                self.threadpoolprog.start(worker_process2)
                self.prog_check.prog_status.connect(self.update_status)
                self.message_prog.message_status.connect(self.update_message)
                worker_process2.signals.finished.connect(self.thread_complete)
                worker_process2.signals.finished.connect(sankey_html)
                self.threadpoolprog.clear()
                return
            def heatmap_png():
                list_from_scratch = []
                for i in os.listdir(dlg_files):
                    if i.endswith("parsed_log.txt"):
                       list_from_scratch.append(str(dlg_files)+"/"+i)
                self.threadpoolprog = QThreadPool()
                self.start_logs = process_results()
                def start_list(working_project):
                    self.start_logs.build_sankey_and_other_logs(list_from_scratch, param="heatmaps")
                worker_process = Worker(start_list, working_project)
                self.prog_check = process_results_progress(param="heatmaps")
                self.message_prog = process_results_progress(param="heatmaps")
                self.prog_check.start()
                self.message_prog.start()
                self.threadpoolprog.start(worker_process)
                self.prog_check.prog_status.connect(self.update_status)
                self.message_prog.message_status.connect(self.update_message)
                worker_process.signals.finished.connect(self.thread_complete)
                worker_process.signals.finished.connect(sankey_png)
                #self.threadpoolprog.clear()
                return
            def count_contacts():
                os.chdir(dlg_files)
                self.threadpoolbina = QThreadPool()
                self.start_log = count_contacts_bina.count()
                #print("THIS IS THE RESULT OF PRINTING A CLASS INSTANCE : ", self.start_log )
                def collect_info(dlg_files):
                     self.start_log.do_the_deed(dlg_files)
                worker_process2 = Worker(collect_info, dlg_files)
                self.prog_check = process_results_progress(param="collect_binana")
                self.message_prog = process_results_progress(param="collect_binana")
                self.prog_check.start()
                self.message_prog.start()
                self.threadpoolbina.start(worker_process2)
                print(self.threadpoolbina.activeThreadCount())
                self.prog_check.prog_status.connect(self.update_status)
                self.message_prog.message_status.connect(self.update_message)
                worker_process2.signals.finished.connect(self.thread_complete)
                worker_process2.signals.finished.connect(heatmap_png)
                self.threadpoolbina.clear()
                return
            def run_binana():
                os.chdir(working_project)
                self.threadpool = QThreadPool()
                maxnum = self.threadpool.maxThreadCount()
                ## I was going to have the user choose how many cores to use with Bianana but I kept having problems with that.
                ## So I changed to 10 or less because it seems to work better.
                #self.text, okPressed = QInputDialog.getText(self, "Choose number of processor to use","Number of processors to use \n or click OK to use the maximum \n Max detected: "+str(maxnum)+ " :", QLineEdit.Normal, "")
                #if okPressed and self.text == "":
                #  self.text = maxnum - 1
                #  print(self.text)
                #elif okPressed and int(self.text) != maxnum and not int(self.text) == 0:
                #  print("number other than 0 or max: ",self.text)
                #elif okPressed and int(self.text) == '0':
                #  self.text = 1
                #else:
                #  self.text = maxnum - 1
                #  print("number max - 1: ",self.text)
                if maxnum > 10:
                    self.text = 10
                else:
                    self.text = maxnum - 1
                self.label_update.setText("Starting")
                worker = Worker_docking(self.execute_binana) # Any other args, kwargs are passed to the run function
                self.prog_check = process_results_progress(param="run_binana")
                self.prog_check.prog_status.connect(self.update_status)
                self.message_prog = process_results_progress(param="run_binana")
                self.message_prog.message_status.connect(self.update_message)
                self.threadpool.start(worker)
                self.prog_check.start()
                self.message_prog.start()
                worker.signals.finished.connect(self.thread_complete)
                worker.signals.finished.connect(count_contacts)
            def start_binana():
                os.chdir(dlg_files)
                binana_list = []
                for dlg_file, subdirs, files in os.walk(dlg_files):
                    for get_dlg_conf in files:
                        if get_dlg_conf.endswith(".dlg") and get_dlg_conf not in binana_list:
                           binana_list.append(os.path.join(dlg_files, get_dlg_conf))
                self.threadpoolbina = QThreadPool()
                self.start_log = process_results()
                def run_binana_list(binana_list, install_dir, dlg_files, protein_files):
                    self.start_log.run_binana(binana_list, install_dir, dlg_files, protein_files)
                worker_process2 = Worker(run_binana_list, binana_list, install_dir, dlg_files, protein_files)
                self.prog_check = process_results_progress(param="run_binana")
                self.message_prog = process_results_progress(param="run_binana")
                self.prog_check.start()
                self.message_prog.start()
                self.threadpoolbina.start(worker_process2)
                self.prog_check.prog_status.connect(self.update_status)
                self.message_prog.message_status.connect(self.update_message)
                worker_process2.signals.finished.connect(self.thread_complete)
                worker_process2.signals.finished.connect(run_binana)
                self.threadpoolbina.clear()
                return
            def prot_mol_key():
                self.threadpoolprotkey = QThreadPool()
                self.start_log = process_results()
                def mol_prot_key(dlg_files):
                     self.start_log.gen_mol_prot_key(dlg_files)
                worker_process2 = Worker(mol_prot_key, dlg_files)
                self.message_prog = process_results_progress(param="protein_key")
                self.message_prog.start()
                self.threadpoolprotkey.start(worker_process2)
                self.message_prog.message_status.connect(self.update_message)
                worker_process2.signals.finished.connect(self.thread_complete)
                worker_process2.signals.finished.connect(start_binana)
                self.threadpoolprotkey.clear()
                self.message_prog.quit()
                return
            def docking_summary():
                all_dlg_undone = []
                for dlg_file_found in os.listdir(dlg_files):
                    if dlg_file_found.endswith(".dlg") and dlg_file_found not in all_dlg_undone:
                       all_dlg_undone.append(dlg_file_found)
                self.threadpooldocksum = QThreadPool()
                self.start_log = process_results()
                def dlg_undone(all_dlg_undone, dlg_files, protein_files, install_dir):
                     self.start_log.make_results_log(all_dlg_undone, dlg_files, protein_files, install_dir)
                worker_process2 = Worker(dlg_undone, all_dlg_undone, dlg_files, protein_files, install_dir)
                self.prog_check = process_results_progress(param="get_score")
                self.message_prog = process_results_progress(param="get_score")
                self.prog_check.start()
                self.message_prog.start()
                self.threadpooldocksum.start(worker_process2)
                self.prog_check.prog_status.connect(self.update_status)
                self.message_prog.message_status.connect(self.update_message)
                worker_process2.signals.finished.connect(self.thread_complete)
                #time.sleep(2)
                worker_process2.signals.finished.connect(prot_mol_key)
                self.threadpooldocksum.clear()
                self.message_prog.quit()
                self.prog_check.quit()
                #self.message_prog.terminate()
                return
            def extract_LE_molecule():
                all_dlg_undone = []
                if platform == "darwin":
                    print("found MAC")
                    list_of_dlg =[]
                    for file in os.listdir(dlg_files):
                        if file.endswith(".dlg"):
                           list_of_dlg.append(file)
                    self.threadpoolextract = QThreadPool()
                    self.start_log = parse_dlg.get_le_mols()
                    def collect_logs(dlg_files, list_of_dlg):
                        self.start_log.collect_and_log(dlg_files, list_of_dlg)
                    worker_process2 = Worker(collect_logs, dlg_files, list_of_dlg)
                    self.prog_check = process_results_progress(param="write_le_pose")
                    self.message_prog = process_results_progress(param="write_le_pose")
                    self.prog_check.start()
                    self.message_prog.start()
                    self.threadpoolextract.start(worker_process2)
                    self.prog_check.prog_status.connect(self.update_status)
                    self.message_prog.message_status.connect(self.update_message)
                    worker_process2.signals.finished.connect(self.thread_complete)
                    worker_process2.signals.finished.connect(prot_mol_key)
                    self.threadpoolextract.clear()
                else:
                    for dlg_file_found in os.listdir(dlg_files):
                        if dlg_file_found.endswith(".dlg") and dlg_file_found not in all_dlg_undone:
                           all_dlg_undone.append(dlg_file_found)
                    self.threadpoolextract = QThreadPool()
                    self.start_log = process_results()
                    def get_poses(all_dlg_undone, dlg_files, install_dir):
                        self.start_log.write_le_pose(all_dlg_undone, dlg_files, install_dir)
                    worker_process2 = Worker(get_poses, all_dlg_undone, dlg_files, install_dir)
                    self.prog_check = process_results_progress(param="write_le_pose")
                    self.message_prog = process_results_progress(param="write_le_pose")
                    self.prog_check.start()
                    self.message_prog.start()
                    self.threadpoolextract.start(worker_process2)
                    self.prog_check.prog_status.connect(self.update_status)
                    self.message_prog.message_status.connect(self.update_message)
                    worker_process2.signals.finished.connect(self.thread_complete)
                    worker_process2.signals.finished.connect(docking_summary)
                    self.threadpoolextract.clear()
                    return
            def start_process_result():
                os.chdir(dlg_files)
                self.threadpooldlgmove = QThreadPool()
                self.start_log = process_results()
                def start_dlgs(dpf_files, dlg_files):
                     self.start_log.find_dlgs(dpf_files, dlg_files)
                worker_process2 = Worker(start_dlgs, dpf_files, dlg_files)
                self.message_prog = process_results_progress(param="rename")
                self.message_prog.start()
                self.threadpooldlgmove.start(worker_process2)
                self.message_prog.message_status.connect(self.update_message)
                worker_process2.signals.finished.connect(self.thread_complete)
                worker_process2.signals.finished.connect(extract_LE_molecule)
                self.threadpooldlgmove.clear()
            start_process_result()
        except:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText("Open or start a project")
            msg.setWindowTitle("Working directory not found")
            msg.setDetailedText("To process docking there needs to be a project directory with a molecule_files, dpf files and dlg files from docking etc. If processing docking stopped for some reason it might be best to try restart processing. This will be generated if a new project is started or found if an existing project is opened and there is a text file with the appropriate information")
            msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
            msg.exec_()
        return

    def restart_process_results(self):
        try:
            restart_lists = []
            for file in os.listdir(working_project):
                if file.startswith("list") and file.endswith(".txt") and file not in restart_lists:
                    restart_lists.append(file)
            def summarize_and_plot():
                    os.chdir(results_files)
                    self.threadpoolsum = QThreadPool()
                    self.threadpoolplot = QThreadPool()
                    self.start_log = process_results()
                    scores_info_file = "docking_summary_protein_key.txt"
                    prediction_info_file = "table_of_prediction.csv"
                    def plot_final(scores_info_file, prediction_info_file):
                        self.start_log.plot_merg_scores_and_prediction(scores_info_file, prediction_info_file)
                    def summarize_interact(dlg_files, protein_files):
                        self.start_log.sumarize_results_and_generate_tables(dlg_files, protein_files)
                    worker_processA = Worker(summarize_interact, dlg_files, protein_files)
                    worker_processB = Worker(plot_final, scores_info_file, prediction_info_file)
                    self.prog_check = process_results_progress(param="final_work")
                    self.message_prog = process_results_progress(param="final_work")
                    self.prog_check.start()
                    self.message_prog.start()
                    self.threadpoolsum.start(worker_processA)
                    self.threadpoolplot.start(worker_processB)
                    self.prog_check.prog_status.connect(self.update_status)
                    self.message_prog.message_status.connect(self.update_message)
                    worker_processB.signals.finished.connect(self.thread_complete)
                    self.threadpoolplot.clear()
                    worker_processA.signals.finished.connect(self.thread_complete)
                    self.threadpoolsum.clear()
                    return
            def run_model():
                    self.threadpoolbina = QThreadPool()
                    self.start_log = process_results()
                    def interact_diagrams(install_dir, results_files, dlg_files):
                        self.start_log.start_predictions(install_dir, results_files, dlg_files)
                 #       self.start_log.run_binana(install_dir, dlg_files, protein_files)
                    worker_process2 = Worker(interact_diagrams, install_dir, results_files, dlg_files)
                    self.prog_check = process_results_progress(param="model_run")
                    self.message_prog = process_results_progress(param="model_run")
                    self.prog_check.start()
                    self.message_prog.start()
                    self.threadpoolbina.start(worker_process2)
                    self.prog_check.prog_status.connect(self.update_status)
                    self.message_prog.message_status.connect(self.update_message)
                    worker_process2.signals.finished.connect(self.thread_complete)
                    worker_process2.signals.finished.connect(summarize_and_plot)
                    self.threadpoolbina.clear()
                    return
            def model_input():
                    for filefound in os.listdir(results_files):
                        file = os.path.join(results_files, filefound)
                        if file.endswith("prediction_key.txt"):
                            os.remove(file)
                        elif file.endswith("interaction_info.csv"):
                            os.remove(file)
                        elif file.endswith("prediction.csv"):
                            os.remove(file)
                        elif file.endswith("image_dataframe.csv"):
                            os.remove(file)
                        else:
                            continue
                    os.chdir(dlg_files)
                    self.threadpoolbina = QThreadPool()
                    self.start_log = process_results()
                    def interact_diagrams(dlg_files, results_files):
                         self.start_log.generate_deep_learning_input(dlg_files, results_files)
                    worker_process2 = Worker(interact_diagrams, dlg_files, results_files)
                    self.prog_check = process_results_progress(param="find_results")
                    self.message_prog = process_results_progress(param="find_results")
                    self.prog_check.start()
                    self.message_prog.start()
                    self.threadpoolbina.start(worker_process2)
                    self.prog_check.prog_status.connect(self.update_status)
                    self.message_prog.message_status.connect(self.update_message)
                    worker_process2.signals.finished.connect(self.thread_complete)
                    worker_process2.signals.finished.connect(run_model)
                    self.threadpoolbina.clear()
                    return
            def sankey_html():
                os.chdir(working_project)
                if "list-left-to-process_sankey-html.txt" in restart_lists:
                    sankey_html_left = open("list-left-to-process_sankey-html.txt", "r").readlines()
                    self.threadpoolprog = QThreadPool()
                    self.start_log = process_results()
                    def interact_diagrams(working_project):
                        self.start_log.build_sankey_and_other_logs(sankey_html_left, param="sankey-html")
                    worker_process2 = Worker(interact_diagrams, working_project)
                    self.prog_check = process_results_progress(param="sankey-html")
                    self.message_prog = process_results_progress(param="sankey-html")
                    self.prog_check.start()
                    self.message_prog.start()
                    self.threadpoolprog.start(worker_process2)
                    self.prog_check.prog_status.connect(self.update_status)
                    self.message_prog.message_status.connect(self.update_message)
                    worker_process2.signals.finished.connect(self.thread_complete)
                    worker_process2.signals.finished.connect(model_input)
                    self.threadpoolprog.clear()
                else:
                    model_input()
                    return
            def sankey_png():
                os.chdir(working_project)
                if "list-left-to-process_sankey-png.txt" in restart_lists:
                    sankey_png_left = open("list-left-to-process_sankey-png.txt", "r").readlines()
                    self.threadpoolprog = QThreadPool()
                    self.start_log = process_results()
                    def interact_diagrams(working_project):
                        self.start_log.build_sankey_and_other_logs(sankey_png_left, param="sankey-png")
                    worker_process2 = Worker(interact_diagrams, working_project)
                    self.prog_check = process_results_progress(param="sankey-png")
                    self.message_prog = process_results_progress(param="sankey-png")
                    self.prog_check.start()
                    self.message_prog.start()
                    self.threadpoolprog.start(worker_process2)
                    self.prog_check.prog_status.connect(self.update_status)
                    self.message_prog.message_status.connect(self.update_message)
                    worker_process2.signals.finished.connect(self.thread_complete)
                    worker_process2.signals.finished.connect(sankey_html)
                    self.threadpoolprog.clear()
                else:
                    sankey_html()
                    return
            def heatmap_png():
                os.chdir(working_project)
                if "list-left-to-process_heatmaps.txt" in restart_lists:
                    heatmaps_left = open("list-left-to-process_heatmaps.txt", "r").readlines()
                    self.threadpoolprog = QThreadPool()
                    self.start_logs = process_results()
                    def interact_diagrams(working_project):
                        self.start_logs.build_sankey_and_other_logs(heatmaps_left, param="heatmaps")
                    worker_process = Worker(interact_diagrams, working_project)
                    self.prog_check = process_results_progress(param="heatmaps")
                    self.message_prog = process_results_progress(param="heatmaps")
                    self.prog_check.start()
                    self.message_prog.start()
                    self.threadpoolprog.start(worker_process)
                    self.prog_check.prog_status.connect(self.update_status)
                    self.message_prog.message_status.connect(self.update_message)
                    worker_process.signals.finished.connect(self.thread_complete)
                    worker_process.signals.finished.connect(sankey_png)
                    self.threadpoolprog.clear()
                else:
                    sankey_png()
                    return
            def count_contacts():
                os.chdir(dlg_files)
                self.threadpoolbina = QThreadPool()
                self.start_log = count_contacts_bina.count()
                #print("THIS IS THE RESULT OF PRINTING A CLASS INSTANCE : ", self.start_log )
                def interact_diagrams(dlg_files):
                     self.start_log.do_the_deed(dlg_files)
                worker_process2 = Worker(interact_diagrams, dlg_files)
                self.prog_check = process_results_progress(param="collect_binana")
                self.message_prog = process_results_progress(param="collect_binana")
                self.prog_check.start()
                self.message_prog.start()
                self.threadpoolbina.start(worker_process2)
                self.prog_check.prog_status.connect(self.update_status)
                self.message_prog.message_status.connect(self.update_message)
                worker_process2.signals.finished.connect(self.thread_complete)
                worker_process2.signals.finished.connect(heatmap_png)
                self.threadpoolbina.clear()
                return
            def run_binana():
                os.chdir(working_project)
                self.threadpool = QThreadPool()
                maxnum = self.threadpool.maxThreadCount()
                #self.text, okPressed = QInputDialog.getText(self, "Choose number of processor to use","Number of processors to use \n or click OK to use the maximum \n Max detected: "+str(maxnum)+ " :", QLineEdit.Normal, "")
                ## I was going to have the user choose how many cores to use with Bianana but I kept having problems with that.
                ## So I changed to 10 or less because it seems to work better.
                #self.text, okPressed = QInputDialog.getText(self, "Choose number of processor to use","Number of processors to use \n or click OK to use the maximum \n Max detected: "+str(maxnum)+ " :", QLineEdit.Normal, "")
                #if okPressed and self.text == "":
                #  self.text = maxnum - 1
                #  print(self.text)
                #elif okPressed and int(self.text) != maxnum and not int(self.text) == 0:
                #  print("number other than 0 or max: ",self.text)
                #elif okPressed and int(self.text) == '0':
                #  self.text = 1
                #else:
                #  self.text = maxnum - 1
                #  print("number max - 1: ",self.text)
                if maxnum > 10:
                    self.text = 10
                else:
                    self.text = maxnum - 1
                self.label_update.setText("Starting")
                worker = Worker_docking(self.execute_binana) # Any other args, kwargs are passed to the run function
                self.prog_check = process_results_progress(param="run_binana")
                self.prog_check.prog_status.connect(self.update_status)
                self.message_prog = process_results_progress(param="run_binana")
                self.message_prog.message_status.connect(self.update_message)
                self.threadpool.start(worker)
                self.prog_check.start()
                self.message_prog.start()
                worker.signals.finished.connect(self.thread_complete)
                worker.signals.finished.connect(count_contacts)
                return
            def start_binana():
                os.chdir(working_project)
                if "list-mols_left_for_binana.txt" in restart_lists:
                    self.threadpoolbina = QThreadPool()
                    self.start_log = process_results()
                    binana_list = open("list-mols_left_for_binana.txt").readlines()
                    os.chdir(dlg_files)
                    def interact_diagrams(binana_list, install_dir, dlg_files, protein_files):
                        self.start_log.run_binana(binana_list, install_dir, dlg_files, protein_files)
                    worker_process2 = Worker(interact_diagrams, binana_list, install_dir, dlg_files, protein_files)
                    self.prog_check = process_results_progress(param="run_binana")
                    self.message_prog = process_results_progress(param="run_binana")
                    self.prog_check.start()
                    self.message_prog.start()
                    self.threadpoolbina.start(worker_process2)
                    self.prog_check.prog_status.connect(self.update_status)
                    self.message_prog.message_status.connect(self.update_message)
                    worker_process2.signals.finished.connect(self.thread_complete)
                    worker_process2.signals.finished.connect(run_binana)
                    self.threadpoolbina.clear()
                elif "list-mols_left_for_counting.txt" in restart_lists and not "list-mols_left_for_binana.txt" in restart_lists:
                    count_contacts()
                else:
                    heatmap_png()
                    return
            def prot_mol_key():
                    print("running prot_mol_key from restart")
                    self.threadpoolprotkey = QThreadPool()
                    self.start_log = process_results()
                    def interact_diagrams(dlg_files):
                         self.start_log.gen_mol_prot_key(dlg_files)
                    worker_process2 = Worker(interact_diagrams, dlg_files)
                    #self.prog_check = process_results_progress(param="get_score")
                    self.message_prog = process_results_progress(param="protein_key")
                    #self.prog_check.start()
                    self.message_prog.start()
                    self.threadpoolprotkey.start(worker_process2)
                    #self.prog_check.prog_status.connect(self.update_status)
                    self.message_prog.message_status.connect(self.update_message)
                    worker_process2.signals.finished.connect(self.thread_complete)
                    worker_process2.signals.finished.connect(run_binana)
                    self.threadpoolprotkey.clear()
                    return
            def docking_summary():
                if "list-mols_left_to_get-info-on.txt" in restart_lists:
                    all_dlg_undone = open("list-mols_left_to_get-info-on.txt", "r").readlines()
                    self.threadpooldocksum = QThreadPool()
                    self.start_log = process_results()
                    def interact_diagrams(all_dlg_undone, dlg_files, protein_files, install_dir):
                         self.start_log.make_results_log(all_dlg_undone, dlg_files, protein_files, install_dir)
                    worker_process2 = Worker(interact_diagrams, all_dlg_undone, dlg_files, protein_files, install_dir)
                    self.prog_check = process_results_progress(param="get_score_redo")
                    self.message_prog = process_results_progress(param="get_score_redo")
                    self.prog_check.start()
                    self.message_prog.start()
                    self.threadpooldocksum.start(worker_process2)
                    self.prog_check.prog_status.connect(self.update_status)
                    self.message_prog.message_status.connect(self.update_message)
                    worker_process2.signals.finished.connect(self.thread_complete)
                    worker_process2.signals.finished.connect(prot_mol_key)
                else:
                    start_binana()
                    return
            def extract_LE_molecule():
                if "list-mols_left_to_extract.txt" in restart_lists:
                    if platform == "darwin":
                        print("found MAC")
                        list_of_dlg = "list-mols_left_to_extract.txt"
                        self.threadpoolextract = QThreadPool()
                        self.start_log = parse_dlg.get_le_mols()
                        def interact_diagrams(dlg_files, list_of_dlg):
                            #print("this is list_of_dlg: ", list_of_dlg)
                            self.start_log.collect_and_log(dlg_files, list_of_dlg)
                        worker_process2 = Worker(interact_diagrams, dlg_files, list_of_dlg)
                        self.prog_check = process_results_progress(param="write_le_pose")
                        self.message_prog = process_results_progress(param="write_le_pose")
                        self.prog_check.start()
                        self.message_prog.start()
                        self.threadpoolextract.start(worker_process2)
                        self.prog_check.prog_status.connect(self.update_status)
                        self.message_prog.message_status.connect(self.update_message)
                        worker_process2.signals.finished.connect(self.thread_complete)
                        worker_process2.signals.finished.connect(prot_mol_key)
                        self.threadpoolextract.clear()
                    else:
                        all_dlg_undone = open("list-mols_left_to_extract.txt", "r").readlines()
                        print("found info left to get ")
                        self.threadpoolextract = QThreadPool()
                        self.start_log = process_results()
                        def interact_diagrams(all_dlg_undone, dlg_files, install_dir):
                            self.start_log.write_le_pose(all_dlg_undone, dlg_files, install_dir)
                        worker_process2 = Worker(interact_diagrams, all_dlg_undone,  dlg_files, install_dir)
                        self.prog_check = process_results_progress(param="write_le_pose")
                        self.message_prog = process_results_progress(param="write_le_pose")
                        self.prog_check.start()
                        self.message_prog.start()
                        self.threadpoolextract.start(worker_process2)
                        self.prog_check.prog_status.connect(self.update_status)
                        self.message_prog.message_status.connect(self.update_message)
                        worker_process2.signals.finished.connect(self.thread_complete)
                        #worker_process2.signals.finished.connect(docking_summary)
                        if "list-mols_left_to_get-info-on.txt" in restart_lists:
                            worker_process2.signals.finished.connect(docking_summary)
                            self.threadpoolextract.clear()
                else:
                    docking_summary()
                return
            extract_LE_molecule()
        except:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText("Open or start a project")
            msg.setWindowTitle("Working directory not found")
            msg.setDetailedText("To restart ""process docking"" there needs to be a project directory with a molecule_files, dpf files and dlg files from docking etc. This will be generated if a new project is started or found if an existing project is opened and there is a text file with the appropriate information")
            msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
            msg.exec_()

    def check_project_status(self):
        try:
             self.threadpoolprog = QThreadPool()
             self.start_log = restart_progress_check()
             def check_proj_status(working_project, install_dir):
                  self.start_log.check_status(working_project, install_dir)
             worker_process2 = Worker(check_proj_status, working_project, install_dir)
             self.prog_check = process_results_progress(param="check_restart")
             self.message_prog = process_results_progress(param="check_restart")
             self.prog_check.start()
             self.message_prog.start()
             self.threadpoolprog.start(worker_process2)
             self.prog_check.prog_status.connect(self.update_status)
             self.message_prog.message_status.connect(self.update_message)
             worker_process2.signals.finished.connect(self.thread_complete)
             self.threadpoolprog.clear()
        except:
             msg = QMessageBox()
             msg.setIcon(QMessageBox.Information)
             msg.setText("Open or start a project")
             msg.setWindowTitle("Working directory not found")
             msg.setDetailedText("To restart ""process docking"" there needs to be a project directory with a molecule_files, dpf files and dlg files from docking etc. This will be generated if a new project is started or found if an existing project is opened and there is a text file with the appropriate information")
             msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
             msg.exec_()

    def createList(self):
        self.listwidget = QListWidget()
        self.listwidget.setWindowTitle("Interaction Diagram")
        self.listwidget.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.listwidget.setMaximumWidth(325)
        self.listwidget2 = QListWidget()
        self.listwidget2.setWindowTitle("Receptor Site")
        self.listwidget2.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.listwidget2.setMaximumWidth(325)
        self.listwidget3 = QListWidget()
        self.listwidget3.setWindowTitle("Selected diagrams")
        self.listwidget3.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.listwidget3.setMaximumWidth(325)

    def Load_list(self, item):
        ###I am going to need to if statements here becuase this will keep loading the diagrams if the button is reclicked
        for item in range(self.listwidget3.count()):
            file_name = self.listwidget3.item(item).text().replace(" with ", "_")+"_sankey.html"
            diag_name = " "+self.listwidget3.item(item).text()
            total_file = os.path.join(dlg_files, file_name)
            self.createbottomRightGroupbox(total_file, diag_name)
            self.formLayout.addRow(self.bottomRightGroupbox)

    def Report_Print(self, item):
        try:
            print_data_mol = []
            print_data_prot = []
            print_data_file = []
            print_data_cap = []
            for item in range(self.listwidget3.count()):
                file_name = self.listwidget3.item(item).text().replace(" with ", "_")+"_sankey.png"
                diag_name = " "+self.listwidget3.item(item).text()
                total_file = os.path.join(dlg_files, file_name)
                split_name =  self.listwidget3.item(item).text().split(" with ")
                print_data_prot.append(split_name[0].replace(" ",""))
                print_data_mol.append(split_name[1].replace(" ",""))
                print_data_file.append(total_file)
                print_data_cap.append(diag_name)
            report_select = pd.DataFrame({
                        "Protein Name" : print_data_prot,
                        "Report Molecule Name" : print_data_mol,
                        #"Original log file name" : ori_log,
                        "File location" : print_data_file,
                        #"Hydrophobic interactions" : hydophob,
                        "Caption" : print_data_cap,
                        })
            ###I am going to have to open up a save dialog to get a unique report name... cause this will only work once as is
            for heat_file in os.listdir(results_files):
                if heat_file.startswith("docking_summary_protein_key_heatmap") and heat_file.endswith(".png"):
                    file_heatmap = os.path.join(results_files, heat_file)
                elif heat_file.startswith("docking_summary_protein_key_score-per-mass_heatmap") and heat_file.endswith(".png"):
                    file_heatmap_mass = os.path.join(results_files, heat_file)
                elif heat_file.startswith("predictions_heatmap") and heat_file.endswith(".png"):
                    file_heatpredmap = os.path.join(results_files, heat_file)
                elif heat_file.startswith("score_predictions_heatmap") and heat_file.endswith(".png"):
                    file_heatmap_scoreden = os.path.join(results_files, heat_file)

            fileName, _ = QFileDialog.getSaveFileName(self,'Save file', None,"All Files (*)")
            if fileName == '':
                pass
            elif fileName != '':
                report_select.to_csv(fileName, mode='w')
                list_len = self.listwidget3.count()
                dlg_temp = dlg_files
                protein_temp = protein_files
                report_of_results_pylatex.generate_report(report_select, list_len, file_heatmap, file_heatmap_mass, file_heatpredmap, file_heatmap_scoreden, dlg_files, protein_files, fileName, working_project)
                del report_select
        except:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText("Open or start a project")
            msg.setWindowTitle("Working directory not found")
            msg.setDetailedText("To print results there needs to be a project directory with a molecule_files, dpf files and dlg files from docking etc. This will be generated if a new project is started or found if an existing project is opened and there is a text file with the appropriate information")
            msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
            msg.exec_()

    def clear_list(self, item):
        for element in self.listwidget3.selectedItems():
            thing = self.listwidget3.row(element)
            self.listwidget3.takeItem(thing)
            self.formLayout.removeRow(thing + 4)

    def populate_lists(self):
        try:
            os.chdir(dlg_files)
            button_list = []
            receptor_list = []
            #print("current working director from populate list function: ", os.getcwd())
            #print("the length or receptor_list: ", len(receptor_list))
            for dir, subdirs, files in os.walk(protein_files):
                for name in files:
                    protein_name = os.fsdecode(name)
                    if protein_name.endswith(".pdbqt"):
                        receptor_list.append(protein_name.replace(".pdbqt",""))
            #print("the length or receptor_list after loading: ", len(receptor_list))
            for file in os.listdir(os.getcwd()):
                if file.endswith(".html"):
                   name = file.split(".html")
                   for file_2ndpart in receptor_list:
                       if file_2ndpart in str(name[0]):
                          step_one = str(name[0]).replace("_sankey", "")
                          name2 =str(step_one).replace(file_2ndpart+"_","")
                          button_list.append(name2)
            local_path = os.getcwd()
            os.chdir(results_files)
            for heat_file in os.listdir(results_files):
                if heat_file.startswith("docking_summary_protein_key_heatmap") and heat_file.endswith(".html"):
                    file_heatmap = os.path.join(results_files, heat_file)
                    self.createtopRightGroupbox2(file_heatmap)
                    self.formLayout.addRow(self.topRightGroupbox2)
                elif heat_file.startswith("docking_summary_protein_key_score-per-mass_heatmap") and heat_file.endswith(".html"):
                    file_heatmap_mass = os.path.join(results_files, heat_file)
                    self.createtopRightGroupbox3(file_heatmap_mass)
                    self.formLayout.addRow(self.topRightGroupbox3)
            for heat_file in os.listdir(results_files):
                if heat_file.startswith("predictions_heatmap") and heat_file.endswith(".html"):
                    file_heatpredmap = os.path.join(results_files, heat_file)
                    self.createtopRightGroupboxpred(file_heatpredmap)
                    self.formLayout.addRow(self.topRightGroupboxpred)
                elif heat_file.startswith("score_predictions_heatmap") and heat_file.endswith(".html"):
                    file_heatmap_scoreden = os.path.join(results_files, heat_file)
                    self.createtopRightGroupboxpred2(file_heatmap_scoreden)
                    self.formLayout.addRow(self.topRightGroupboxpred2)
            for receptor in receptor_list:
                #QListWidgetItem(receptor, self.listwidget2)
                #print("list of receptors from populate list fucntion: ", receptor)
                self.listwidget2.addItem(receptor)
            for file in set(button_list):
                #QListWidgetItem(file, self.listwidget)
                self.listwidget.addItem(file)

        except:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText("Open or start a project")
            msg.setWindowTitle("Working directory not found")
            msg.setDetailedText("Get the results there needs to be a project directory with a molecule_files, dpf files and dlg files from docking etc. This will be generated if a new project is started or found if an existing project is opened and there is a text file with the appropriate information")
            msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
            msg.exec_()

        def prot_list_thing(item, mol_name):
           prot_name = self.listwidget2.currentItem().text()
           file_name = prot_name+"_"+mol_name+"_sankey.html"
           diag_name = " "+prot_name+" with "+mol_name
           total_file = os.path.join(dlg_files, file_name)
           QListWidgetItem(prot_name+" with "+mol_name, self.listwidget3)
           del item, mol_name
           return
        def mol_list_thing(item):
            new_var = self.listwidget.currentItem()
            self.listwidget2.itemClicked.connect(lambda item: prot_list_thing(item, new_var.text()))
            if new_var == True:
                del new_var
            else:
                pass
            return
        self.listwidget.itemClicked.connect(mol_list_thing)

    def remove_docking_project(self):
        fileName = QFileDialog.getExistingDirectory(self,"Remove Project")
        if fileName != '':
            project_path = os.getcwd()
            path, file = os.path.split(fileName)
            try:
                os.remove(path+"/tox_projects/"+file+".py")
            except:
                pass
            shutil.rmtree(fileName)
        else:
            return

    def addcolumn(self):
        text, okPressed = QInputDialog.getText(self, "Get text","Text:", QLineEdit.Normal, "")
        if okPressed and text != '':
            new_head = QTableWidgetItem(text)
        else:
            return
        columnPosition = self.tableWidget.columnCount()
        self.tableWidget.insertColumn(columnPosition)

        self.tableWidget.setHorizontalHeaderItem(columnPosition, new_head)

    def addrow(self):
        rowPosition = self.tableWidget.rowCount()
        self.tableWidget.insertRow(rowPosition)

    def edit_toggle_wrap(self):
        self.editor.setLineWrapMode( 1 if self.editor.lineWrapMode() == 0 else 0 )

    def keyPressEvent(self, event):
        if event.matches(QKeySequence.Paste):
            self.paste_data()
        else:
            return

    def keyPressEven(self, event):
        if event.matches(QKeySequence.Copy):
            self.copy_data()
        else:
            return
    @pyqtSlot()
    def on_click(self):
        for currentQTableWidgetItem in self.tableWidget.selectedItems():
            print(currentQTableWidgetItem.row(), currentQTableWidgetItem.column(), currentQTableWidgetItem.text())
    def paste_data(self):
        #There is a problem that I need to fix. if there is an empty line pasteing will only work for the stuff after that line and if there is an empty line at the end it will only put in empy lines
        #selection = self.tableWidget.selectedIndexes()
        test = ""
        for element in self.tableWidget.selectedItems():
            test = element.text()
        #return(test)
        selection = self.tableWidget.selectedIndexes()
        if selection and test == '':
          self.tableWidget.setRowCount(0)
          self.tableWidget.setColumnCount(0)
          paste_info = open("temp_data.txt", 'w')
          paste_info.write(QApplication.clipboard().text())
          paste_info.close()
          for currentQTableWidgetItem in self.tableWidget.selectedItems():
              print(currentQTableWidgetItem.row(), currentQTableWidgetItem.column(), currentQTableWidgetItem.text())

          with open("temp_data.txt", 'r') as stream:
               file_data = csv.reader(stream, delimiter=',')
               self.tableWidget.setRowCount(0)
               self.tableWidget.setColumnCount(0)
               for rowdata in csv.reader(stream, delimiter=','):
                   row = self.tableWidget.rowCount()
                   self.tableWidget.insertRow(row)
                   self.tableWidget.setColumnCount(len(rowdata))
                   for column, data in enumerate(rowdata):
                        item = QTableWidgetItem(data)
                        self.tableWidget.setItem(row, column, item)
                   self.tableWidget.setHorizontalHeaderLabels(["Compound", "Name"])
                   self.tableWidget.setColumnCount(10)

        return

    def copy_data(self):
         ##right now copy_data only captures the last line
         rowdata = []
         first = True
         i = 0
         item_copied = ""
         for item in self.tableWidget.selectedItems():
             item_copied += item.text()
             item_copied += '\n'
             QApplication.clipboard().setText(item_copied)
             i = i + 1
             #print(item_copied, file=open("temp_data.txt", "w"))
             if first:
                copied_info = open("temp_data.txt", "w")
                copied_info.write(item_copied+'\n')
                first = False
                #copied_info.close()
             #QApplication.clipboard().setText(item_copied)
             #if first:
            #    copied_info = open("temp_data.txt", "w")
            #    first = False
            #    copied_info.write(item_copied + '\n')
                #copied_info.close()
                #print(item_copied)
             else:
                copied_info = open("temp_data.txt", "a")
                copied_info.write(item_copied + '\n')
               #QApplication.clipboard().setText(item_copied)
                copied_info.close()

    def createtopRightGroupBox(self):
        self.topRightGroupBox = QGroupBox()
        self.formLayout = QFormLayout()
        self.topRightGroupBox.setLayout(self.formLayout)
        self.scroll = QScrollArea()
        self.scroll.setWidgetResizable(True)
        self.scroll.setWidget(self.topRightGroupBox)
        self.newlayout = QVBoxLayout()
        self.newlayout.addWidget(self.scroll)
        self.show()

    def createtopLeftGroupBox(self):
        self.topLeftGroupBox = QGroupBox()

        defaultPushButton3 = QPushButton("Process docking logs")
        defaultPushButton3.setDefault(False)
        defaultPushButton4 = QPushButton("Populate List")
        defaultPushButton4.setDefault(False)
        clearListPushButton = QPushButton("Remove From List")
        clearListPushButton.setDefault(False)
        loadListPushButton = QPushButton("Load listed interaction diagrams")
        loadListPushButton.setDefault(False)
        printreportPushButton = QPushButton("Print report of display molecules")
        printreportPushButton.setDefault(False)


        self.verticalLayoutWidget = QtWidgets.QWidget()
        #self.verticalLayoutWidget.setGeometry(QRect(0, 450, 741, 131))
        #self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        #self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        #self.verticalLayout.setObjectName("verticalLayout")

        self.label_loc = QtWidgets.QLabel()#self.verticalLayoutWidget)
        self.label_loc.setObjectName("label")
        self.label_loc.setText("No project open")
        self.label_loc.setWordWrap(True)
        self.verticalLayout.addWidget(self.label_loc)

        self.label_update = QtWidgets.QLabel()#self.verticalLayoutWidget)
        self.label_update.setObjectName("label")
        self.label_update.setText("-")
        self.verticalLayout.addWidget(self.label_update)

        self.progressBar = QtWidgets.QProgressBar()#self.verticalLayoutWidget)
        self.progressBar.setProperty("value", 0)
        self.progressBar.setObjectName("progressBar")
        self.verticalLayout.addWidget(self.progressBar)

        defaultPushButton3.clicked.connect(self.process_results)
        defaultPushButton4.clicked.connect(self.populate_lists)
        clearListPushButton.clicked.connect(self.clear_list)
        loadListPushButton.clicked.connect(self.Load_list)
        printreportPushButton.clicked.connect(self.Report_Print)

        self.createTreeView()
        self.createTabs()
        self.createList()

        layout_tabs = QVBoxLayout()
        ##I took self out of these three and solve and the warning below
        ##QLayout: Attempting to add QLayout "" to MainWindow "", which already has a layout
        ##I don't expect this to cause any problems in the future but I am going to keep the note and code for a while until I am sure
        layout3 = QGridLayout()#self)
        layout2 = QVBoxLayout()#self)
        layout_buttons = QVBoxLayout()#self)

        self.tab1.setLayout(layout2)
        self.tab2.setLayout(layout3)

        layout3.addLayout(layout_buttons, 1,1)
        layout3.addWidget(self.listwidget,0,0)
        layout3.addWidget(self.listwidget2,1,0)
        layout3.addWidget(self.listwidget3,0,1)

        layout_buttons.addWidget(defaultPushButton4)
        layout_buttons.addWidget(clearListPushButton)
        layout_buttons.addWidget(loadListPushButton)
        layout_buttons.addWidget(printreportPushButton)

        layout2.addWidget(self.tree)

        layout2.addWidget(self.verticalLayoutWidget)

        layout_tabs.addWidget(self.tabs)
        self.topLeftGroupBox.setLayout(layout_tabs)

        self.newlayoutleft = QVBoxLayout()
        self.newlayoutleft.addWidget(self.topLeftGroupBox)
        return

    def createTable(self):
       # Create table
        self.tableWidget = QTableWidget()
        self.tableWidget.setRowCount(100)
        self.tableWidget.setColumnCount(100)
        self.tableWidget.setHorizontalHeaderLabels(["Compound", "Name"])

    def createTabs(self):
        self.tabs = QTabWidget()
        self.tab1 = QWidget()
        self.tab2 = QWidget()
        self.tabs.resize(100,200)
        self.tabs.addTab(self.tab1,"Project info")
        self.tabs.addTab(self.tab2,"Docking Results and Maching Learning Predicitons")

    def createTreeView(self):
        path = os.getcwd()
        self.model = QFileSystemModel()
        self.model.setRootPath(os.getcwd())

        self.tree = QTreeView()
        self.tree.setGeometry(QRect(0, 10, 721, 431))
        self.tree.setModel(self.model)
        self.tree.setRootIndex(self.model.index(path))
        item_1 = self.tree.selectionModel().currentIndex()

        def start_editor(signal):
            ##find the index of the item clicked
            item_1 = self.tree.selectionModel().currentIndex()
            ##Use the index to get more information using pqt5 functions
            filePath = self.model.filePath(item_1)
            ##open the files with system available programs
            if platform == "win32":
                print("windows was found")
                os.startfile(filePath)
            elif platform == "darwin":
                            opener = "open"
                            subprocess.call([opener, filePath])
            elif platform == "linux" or "linux2":
                print("linux was found")
                subprocess.call(["xdg-open", filePath])
        self.tree.doubleClicked.connect(start_editor)

    def createtopRightGroupbox2(self, file_heatmap):
        self.topRightGroupbox2 = QGroupBox("Heatmap of docking energies")
        view_heatmap = QtWebEngineWidgets.QWebEngineView()
        view_heatmap.load(QUrl.fromLocalFile(file_heatmap))
        self.topRightGroupbox2.setMinimumSize(QSize(0, 650))
        layout = QGridLayout()
        layout.addWidget(view_heatmap)
        self.topRightGroupbox2.setLayout(layout)

    def createtopRightGroupbox3(self, file_heatmap_mass):
        self.topRightGroupbox3 = QGroupBox("Heatmap of docking energies per exact molecular mass")
        view_heatmap_mass = QtWebEngineWidgets.QWebEngineView()
        view_heatmap_mass.load(QUrl.fromLocalFile(file_heatmap_mass))
        self.topRightGroupbox3.setMinimumSize(QSize(0, 650))
        layout = QGridLayout()
        layout.addWidget(view_heatmap_mass)
        self.topRightGroupbox3.setLayout(layout)

    def createtopRightGroupboxpred(self, file_heatpredmap):
        self.topRightGroupboxpred = QGroupBox("Activity predictions")
        view_pred_heatmap = QtWebEngineWidgets.QWebEngineView()
        view_pred_heatmap.load(QUrl.fromLocalFile(file_heatpredmap))
        self.topRightGroupboxpred.setMinimumSize(QSize(0, 650))
        layout = QGridLayout()
        layout.addWidget(view_pred_heatmap)
        self.topRightGroupboxpred.setLayout(layout)

    def createtopRightGroupboxpred2(self, file_heatmap_scoreden):
        self.topRightGroupboxpred2 = QGroupBox("Activity perdictions per score density")
        view_heatmap_mass = QtWebEngineWidgets.QWebEngineView()
        view_heatmap_mass.load(QUrl.fromLocalFile(file_heatmap_scoreden))
        self.topRightGroupboxpred2.setMinimumSize(QSize(0, 650))
        layout = QGridLayout()
        layout.addWidget(view_heatmap_mass)
        self.topRightGroupboxpred2.setLayout(layout)

    def createbottomRightGroupbox(self, file, diag_name):
        self.bottomRightGroupbox = QGroupBox("Interaction Diagram"+diag_name)
        view_sankey = QtWebEngineWidgets.QWebEngineView()
        view_sankey.load(QUrl.fromLocalFile(file))
        self.bottomRightGroupbox.setMinimumSize(QSize(0, 650))
        layout = QGridLayout()
        layout.addWidget(view_sankey)
        self.bottomRightGroupbox.setLayout(layout)




if __name__ == '__main__':
    app = QApplication(sys.argv)
    ##ex = App()
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
