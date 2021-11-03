#!/usr/bin/env python

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

import re
import os
import sys
from os import path
import glob
import errno
from collections import Counter


    #dlg_path = "/media/mike/disk1/pangburn_lab/making_pipeline/neuro_sites_tox_pipline/results_files/dlg_files/"
    #log_path = "./*_full_log.txt"
class count():
    def __init__(self, *args, **kwargs):
        #super.__init__(self)
        #self.working_project = working_project
        self.args = args
        self.kwargs = kwargs

        inRecordingMode = False
        inRecordingMode2 = False
        inRecordingMode3 = False
        inRecordingMode4 = False
        inRecordingMode5 = False
        inRecordingMode6 = False
        inRecordingMode7 = False

        myline_dist1 = []
        myline_dist2 = []
        myline_dist3 = []
        myline_dist4 = []
        myline_dist5 = []
        myline_dist6 = []
        myline_dist7 = []

        long_output1 = ""
        long_output2 = ""
        long_output3 = ""
        long_output4 = ""
        long_output5 = ""
        long_output6 = ""
        long_output7 = ""

        #files = glob.glob(dlg_path)

    def find_info_in_file(self, get_dlg_conf):

        inRecordingMode = False
        inRecordingMode2 = False
        inRecordingMode3 = False
        inRecordingMode4 = False
        inRecordingMode5 = False
        inRecordingMode6 = False
        inRecordingMode7 = False
        inRecordingMode8 = False

        myline_dist1 = []
        myline_dist2 = []
        myline_dist3 = []
        myline_dist4 = []
        myline_dist5 = []
        myline_dist6 = []
        myline_dist7 = []
        myline_dist8 = []
        myline_dist9 = myline_dist4

        long_output1 = ""
        long_output2 = ""
        long_output3 = ""
        long_output4 = ""
        long_output5 = ""
        long_output6 = ""
        long_output7 = ""
        long_output8 = ""
        long_output9 = ""

        try:

               with open(get_dlg_conf) as dlg_log: # No need to specify 'r': this is the default.
                    for line in dlg_log:
                      if line.startswith('Atom-type pair counts within 2.5 angstroms:'):
                         inRecordingMode = True
                      elif line.startswith('Atom-type pair counts within 4.0 angstroms:'):
                         inRecordingMode = False
                      if inRecordingMode:
                         myline_dist1.append(line)
                      if line.startswith('Atom-type pair counts within 4.0 angstroms:'):
                         inRecordingMode2 = True
                      elif line.startswith('Summed electrostatic energy by atom-type pair, in J/mol:'):
                         inRecordingMode2 = False
                      if inRecordingMode2:
                         myline_dist2.append(line)
                      if line.startswith('Hydrogen bonds:'):
                         inRecordingMode3 = True
                      elif line.startswith('Hydrophobic contacts (C-C):'):
                         inRecordingMode3 = False
                      if inRecordingMode3:
                         myline_dist3.append(line)
                      if line.startswith('Hydrophobic contacts (C-C):'):
                         inRecordingMode4 = True
                      elif line.startswith('pi-pi stacking interactions:'):
                         inRecordingMode4 = False
                      if inRecordingMode4:
                         myline_dist4.append(line)
                      if line.startswith('pi-pi stacking interactions:'):
                         inRecordingMode5 = True
                      elif line.startswith('T-stacking (face-to-edge) interactions:'):
                         inRecordingMode5 = False
                      if inRecordingMode5:
                         myline_dist5.append(line)
                      if line.startswith('T-stacking (face-to-edge) interactions:'):
                         inRecordingMode6 = True
                      elif line.startswith('Cation-pi interactions:'):
                         inRecordingMode6 = False
                      if inRecordingMode6:
                         myline_dist6.append(line)
                      if line.startswith('Cation-pi interactions:'):
                         inRecordingMode7 = True
                      elif line.startswith('Salt Bridges:'):
                         inRecordingMode7 = False
                      if inRecordingMode7:
                         myline_dist7.append(line)
                      if line.startswith('Salt Bridges:'):
                         inRecordingMode8 = True
                      if inRecordingMode8:
                         myline_dist8.append(line)

                    for AAR in ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']:
                       prot_line9 = re.findall(AAR+"\W+\w+\W", str(myline_dist2))
                       if prot_line9:
                          amino_residue9 = prot_line9
                          long_output9 += (str(Counter(amino_residue9)).strip('Counter({})').replace("'","")+", ")
                    if long_output9 != "":
                       long_output10 = re.sub(": \d+","", str(long_output9))
                       print(("binding-interactions, %s" % (long_output10.rstrip(', '))), file=open(get_dlg_conf.replace("full_log.txt","parsed_log.txt"), "w")) #I specify write "w" here because if this script is run twice it will generate new files. #This avoid an error cause by adding data to the old file that the next script can't handle.
                    else:
                       print(("binding-interactions, %s" % '0: 0'), file=open(get_dlg_conf.replace("full_log.txt","parsed_log.txt"), "a"))
                    for AAR in ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']:
                        prot_line = re.findall(AAR+"\W+\w+\W", str(myline_dist1))
                        if prot_line:
                           amino_residue = prot_line
                           long_output1 += (str(Counter(amino_residue)).strip('Counter({})').replace("'","")+", ")
                    if long_output1 != "":
                       print(("Short-range contacts, %s" % (long_output1.rstrip(', '))), file=open(get_dlg_conf.replace("full_log.txt","parsed_log.txt"), "a"))
                    else:
                       print(("Short-range contacts, %s" % '0: 0'), file=open(get_dlg_conf.replace("full_log.txt","parsed_log.txt"), "a"))
                    for AAR in ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']:
                        prot_line3 = re.findall(AAR+"\W+\w+\W", str(myline_dist3))
                        if prot_line3:
                           amino_residue3 = prot_line3
                           long_output3 += (str(Counter(amino_residue3)).strip('Counter({})').replace("'","")+", ")
                    if long_output3 != "":
                       print(("Hydrogen bonds, %s" % (long_output3.rstrip(', '))), file=open(get_dlg_conf.replace("full_log.txt","parsed_log.txt"), "a"))
                    else:
                       print(("Hydrogen bonds, %s" %  '0: 0'), file=open(get_dlg_conf.replace("full_log.txt","parsed_log.txt"), "a"))
                    for AAR in ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']:
                       prot_line4 = re.findall(AAR+"\W+\w+\W", str(myline_dist4))
                       if prot_line4:
                          amino_residue4 = prot_line4
                          long_output4 += (str(Counter(amino_residue4)).strip('Counter({})').replace("'","")+", ")
                    if long_output4 != "":
                       print(("Hydrophobic contacts, %s" % (long_output4.rstrip(', '))), file=open(get_dlg_conf.replace("full_log.txt","parsed_log.txt"), "a"))
                    else:
                       print(("Hydrophobic contacts, %s" % '0: 0'), file=open(get_dlg_conf.replace("full_log.txt","parsed_log.txt"), "a"))
                    for AAR in ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']:
                        prot_line5 = re.findall(AAR+"\W+\w+\W", str(myline_dist5))
                        if prot_line5:
                           amino_residue5 = prot_line5
                           long_output5 += (str(Counter(amino_residue5)).strip('Counter.({})').replace("'","")+", ")
                    if long_output5 != "":
                       print(("pi-pi stacking interactions, %s" % (long_output5.rstrip(', '))), file=open(get_dlg_conf.replace("full_log.txt","parsed_log.txt"), "a"))
                    else:
                       print(("pi-pi stacking interactions, %s" % '0: 0'), file=open(get_dlg_conf.replace("full_log.txt","parsed_log.txt"), "a"))
                    for AAR in ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']:
                        prot_line6 = re.findall(AAR+"\W+\w+\W", str(myline_dist6))
                        if prot_line6:
                           amino_residue6 = prot_line6
                           long_output6 += (str(Counter(amino_residue6)).strip('Counter({})').replace("'","").replace("]","").replace("[","")+", ")
                    if long_output6 != "":
                       print(("T-stacking, %s" % (long_output6.rstrip(', '))), file=open(get_dlg_conf.replace("full_log.txt","parsed_log.txt"), "a"))
                    else:
                         print(("T-stacking, %s" % '0: 0'), file=open(get_dlg_conf.replace("full_log.txt","parsed_log.txt"), "a"))
                    for AAR in ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']:
                        prot_line7 = re.findall(AAR+"\W+\w+\W", str(myline_dist7))
                        if prot_line7:
                           amino_residue7 = prot_line7
                           long_output7 += (str(Counter(amino_residue7)).strip('Counter({})').replace("'","").replace("]","").replace("[","")+", ")
                    if long_output7 != "":
                       print(("cation-pi, %s" % (long_output7.rstrip(', '))), file=open(get_dlg_conf.replace("full_log.txt","parsed_log.txt"), "a"))
                         #del dlg_log
                    else:
                       print(("cation-pi, %s" % '0: 0'), file=open(get_dlg_conf.replace("full_log.txt","parsed_log.txt"), "a"))
                    for AAR in ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']:
                        prot_line8 = re.findall(AAR+"\W+\w+\W", str(myline_dist8))
                        if prot_line8:
                           amino_residue8 = prot_line8
                           long_output8 += (str(Counter(amino_residue8)).strip('Counter({})').replace("'","").replace("]","").replace("[","")+", ")
                    if long_output8 != "":
                       print(("Salt Bridges, %s" % (long_output8.rstrip(', '))), file=open(get_dlg_conf.replace("full_log.txt","parsed_log.txt"), "a"))
                         #del dlg_log
                    else:
                       print(("Salt Bridges, %s" % '0: 0'), file=open(get_dlg_conf.replace("full_log.txt","parsed_log.txt"), "a"))
        except IOError as exc:
           if exc.errno != errno.EISDIR: # Do not fail if a directory is found, just ignore it.
              raise # Propagate other kinds of IOError.
           else:
              del dlg_log
              pass
           return
    def do_the_deed(self, dlg_files):
        for dlg_file, subdirs, files in os.walk(dlg_files):
            #os.chdir(dlg_file)
            for get_dlg_conf in files: #for each file in files found
                if get_dlg_conf.endswith("_full_log.txt"): # for each file in the directory, if the file ends with _full_log.txt
                   self.find_info_in_file(get_dlg_conf)
        return



    #['Atom-type pair counts within 2.5 angstroms', 'Atom-type pair counts within 4.0 angstroms', 'Hydrogen bonds', 'Hydrophobic contacts', 'pi-pi stacking interactions', 'T-stacking (face-to-edge) interactions', 'Cation-pi interactions', 'Salt Bridges']
