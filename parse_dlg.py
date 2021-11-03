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
import pandas as pd

#input = sys.argv[-1]

#print("This is input: ", input)
class get_le_mols():
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
        print("Print anything from find_info_in_file")
        print("This is get_dlg_conf: ", get_dlg_conf)
        inRecordingMode = False
        inRecordingMode2 = False
        inRecordingMode3 = False
        inRecordingMode4 = False
        inRecordingMode5 = False
        inRecordingMode6 = False
        inRecordingMode7 = False
        inRecordingMode8 = False
        inRecordingMode9 = False
        inRecordingMode10 = False

        myline_dist1 = []
        myline_dist2 = []
        myline_dist3 = []
        myline_dist4 = []
        myline_dist5 = []
        myline_dist6 = []
        myline_dist7 = []
        myline_dist8 = []
        myline_dist9 = []
        myline_dist10 = []

        long_output1 = ""
        long_output2 = ""
        long_output3 = ""
        long_output4 = ""
        long_output5 = ""
        long_output6 = ""
        long_output7 = ""
        long_output8 = ""
        long_output9 = ""
        long_output10 = ""

        #ry:

        with open(get_dlg_conf) as dlg_log: # No need to specify 'r': this is the default.
            for line in dlg_log:
              if line.startswith('DOCKED: USER    Run = 1') and not line.startswith('DOCKED: USER    Run = 10'):
                 inRecordingMode = True
              elif line.startswith('DOCKED: ENDMDL'):
                 inRecordingMode = False
              if inRecordingMode:
                 myline_dist1.append(line)
              if line.startswith('DOCKED: USER    Run = 2'):
                 inRecordingMode2 = True
              elif line.startswith('DOCKED: ENDMDL'):
                 inRecordingMode2 = False
              if inRecordingMode2:
                 myline_dist2.append(line)
              if line.startswith('DOCKED: USER    Run = 3'):
                 inRecordingMode3 = True
              elif line.startswith('DOCKED: ENDMDL'):
                 inRecordingMode3 = False
              if inRecordingMode3:
                 myline_dist3.append(line)
              if line.startswith('DOCKED: USER    Run = 4'):
                 inRecordingMode4 = True
              elif line.startswith('DOCKED: ENDMDL'):
                 inRecordingMode4 = False
              if inRecordingMode4:
                 myline_dist4.append(line)
              if line.startswith('DOCKED: USER    Run = 5'):
                 inRecordingMode5 = True
              elif line.startswith('DOCKED: ENDMDL'):
                 inRecordingMode5 = False
              if inRecordingMode5:
                 myline_dist5.append(line)
              if line.startswith('DOCKED: USER    Run = 6'):
                 inRecordingMode6 = True
              elif line.startswith('DOCKED: ENDMDL'):
                 inRecordingMode6 = False
              if inRecordingMode6:
                 myline_dist6.append(line)
              if line.startswith('DOCKED: USER    Run = 7'):
                 inRecordingMode7 = True
              elif line.startswith('DOCKED: ENDMDL'):
                 inRecordingMode7 = False
              if inRecordingMode7:
                 myline_dist7.append(line)
              if line.startswith('DOCKED: USER    Run = 8'):
                 inRecordingMode8 = True
              elif line.startswith('DOCKED: ENDMDL'):
                 inRecordingMode8 = False
              if inRecordingMode8:
                 myline_dist8.append(line)
              if line.startswith('DOCKED: USER    Run = 9'):
                 inRecordingMode9 = True
              elif line.startswith('DOCKED: ENDMDL'):
                 inRecordingMode9 = False
              if inRecordingMode9:
                 myline_dist9.append(line)
              if line.startswith('DOCKED: USER    Run = 10'):
                 inRecordingMode10 = True
              elif line.startswith('DOCKED: ENDMDL'):
                 inRecordingMode10 = False
              if inRecordingMode10:
                 myline_dist10.append(line)
            list_dlg_results = [myline_dist1,myline_dist2,myline_dist3,myline_dist4,myline_dist5,myline_dist6,myline_dist7,myline_dist8,myline_dist9,myline_dist10]
            score_dict = {}
            score_list = []
            dlg_log.close()
        for results in list_dlg_results:
            print("this is results: ", len(results))
            for line in results:
                print("this is line: ", type(line))
                #dock_score = re.search("Estimated Free Energy of Binding    = (\S\d+.+)", line)
                #dock_score = re.search("Estimated Free Energy of Binding    = (\D+\d+\D+\S+\s+\w+\D\w+)", line)
                dock_score_only = re.search("Estimated Free Energy of Binding    = (\D+\d+\D+\S+)", line)
                if dock_score_only:
                    score = float(dock_score_only.group(1))
                    score_dict[list_dlg_results.index(results)] = score
                    if score not in score_list:
                        score_list.append(score)


                    #if score < -10:
                    #    print("this is score_dict: ", score_dict)
                    print("this is score ONLY: ", score)
                    #    print("this is the lowest item in score_dict: ", min(score_dict, key=score_dict.get))
        #print("this is the lowest item in score_dict: ", min(score_dict, key=score_dict.get))
        if not os.path.isfile("docking_summary_log.txt"):
            docking_file = pd.DataFrame(columns=['lowestEnergy_dlgfn','LE'])
            docking_file.to_csv("docking_summary_log.txt", index=False)
        else:
            pass
        docking_file = pd.read_csv("docking_summary_log.txt", delimiter=',', header=0, engine='python')
        #add_data = pd.DataFrame([get_dlg_conf.replace(".dlg",""),min(score_list)], columns=list('lowestEnergy_dlgfn''LE'))
        #docking_file.append(add_data)
        add_data={'lowestEnergy_dlgfn' : 'test', 'LE' : '-10'}
        #print("this is tyep get_dlg_conf.replace(\".dlg\",\"\"): ", type(get_dlg_conf.replace(".dlg","")))
        #docking_file['LE'] = min(score_list)
        #print("this is type min(score_list): ", min(score_list))
        docking_file = docking_file.append({'lowestEnergy_dlgfn' : get_dlg_conf.replace(".dlg",""), 'LE' : min(score_list) }, ignore_index=True, sort=False)
        print("this is add_data: ", type(add_data))
        print("this is docking_file: ", docking_file)
        docking_file.to_csv("docking_summary_log.txt",index=False)

        mol_to_extract = min(score_dict, key=score_dict.get)
        output = list_dlg_results[mol_to_extract]
        inRecordingModelowmol = False
        lowmol = []
        for line in output:
            if line.startswith('DOCKED: REMARK  Name'):
                inRecordingModelowmol = True
            elif line.startswith('DOCKED: TER'):
                inRecordingModelowmol = False
            if inRecordingModelowmol:
                lowmol.append(line)
        finallowmol=[]
        for line in lowmol:
            finallowmol.append(line.replace("DOCKED: ", "").replace("USER","REMARK"))

        finallowmol.append("TER")

        writefile = open(get_dlg_conf.replace(".dlg",".pdbqt"), "w")
        for i in finallowmol:
            writefile.write(i)
        writefile.close()


    def collect_and_log(self, dlg_files, list_of_dlg):
        print("AM I STARTING FROM PARSE")
        #print("this is list_of_dlg from parse_log: ", list_of_dlg)
        os.chdir(dlg_files)
        for get_dlg_conf in list_of_dlg:
            print("this is file from the first loop in parse_dlg: ", get_dlg_conf)
            self.find_info_in_file(get_dlg_conf)
        return
