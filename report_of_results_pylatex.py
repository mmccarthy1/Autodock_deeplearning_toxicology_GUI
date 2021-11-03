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

import numpy as np

from pylatex import Document, Section, Subsection, Tabular, Math, TikZ, Axis, \
    Plot, Figure, Matrix, Alignat, SubFigure, NoEscape, Command, LongTable, MultiColumn, LongTabu, NewPage
from pylatex.utils import italic
import os
import pandas as pd
import re
import shutil

if __name__ == '__main__':

    image_filename = os.path.join(os.path.dirname(__file__), 'predictions_heatmap.png')


def generate_text_report(doc):
    local_path = os.getcwd()
    list_of_figures = []
    for find_parsed_heatmap_file in os.listdir(local_path):
        if find_parsed_heatmap_file.endswith(".png"):
                print("this is the output of find_parsed_heatmap_file :", find_parsed_heatmap_file)
                list_of_figures.append(find_parsed_heatmap_file)


    list_of_proteins = []
    list_of_sankey = []
    list_of_logs = []
    for dir, subdirs, files in os.walk(protein_files):
        for name in files:
            protein_name = os.fsdecode(name)
            if protein_name.endswith(".pdbqt"):
                list_of_proteins.append(protein_name.replace(".pdbqt",""))

    for sankey_pics_location in os.listdir(dlg_files):
        if sankey_pics_location.endswith("sankey.png"):
            sankey_full_files = os.path.join(dlg_files, sankey_pics_location)
            list_of_sankey.append(sankey_pics_location)

    for mol_log_location in os.listdir(dlg_files):
        if mol_log_location.endswith("parsed_log.txt"):
            sankey_full_files = os.path.join(dlg_files, mol_log_location)
            list_of_logs.append(mol_log_location)

    with doc.create(Section('Rules that were used to define key residues and interactions on which to build models')):
        doc.append("Criteria for active and inactive molecules for each model.\
\n 1. Available binding interaction information from the literature was used to define key residues.\
\n Active molecules\
\n 2. Use residues that form a salt-bridge, pi interactions or hydrogen bonds etc. with an antagonist or agonist. \
\n 3. If other key residues are necessary for binding or activation then use those as well.\
\n 4. Use all molecules that have both interactions listed above if possible but prioritize the key salt-bridge and hydrogen bond. \
\nInactive molecules\
\n5. Use high docking score molecules that don't have the interactions above where possible. \
\n6. Then look for nearest neighbors using openmolecule datawarrior and used these molecules (matched by number to actives i.e. if actives is 900 then use the top 900 by neighbor count). \
\n7. First the set of inactive molecules with the most neighbors were used.\
\nAll molecules\
\n8. If possible I prioritized the molecules that match the receptor type binding information using the information in the binding database.\
\n9. Eliminated low measured affinity molecules (over 100nM) using the information in the binding database. \
\n10. If possible use molecule with docking score below -6.00kcal/mol\
\n \
\nNaming convention used for the models and proteins.\
\n  Proteins and binding sites are named by the following: PDBID_receptor-type_ActiveSiteMolecule\
\n  For example: 2FY4_CHAT_asCOA is PDB ID 2FY4, receptor-type is choline acetyltransferase(CHAT), and active site(as) molecule is coenzyme(COA)"
)

        with doc.create(Subsection('The model for human choline acetyltransferase (ChAT) bound to coenzyme A PDB: 2FY4 asCOA')):
            doc.append("Active molecule interactions for models\
\n1. CoA forms hydrogen bonds with GLN(541), SER(412), SER(440), ASP(414) and LYS(407), CoA also forms a hydrogen bond with GLN(144) of the p-loop which is described as an important interaction.\
\n2. Active molecule for the model were selected based on hydrogen bond with any of these residues and special consideration was given to molecule with paired hydrogen bonds with GLN(144) or SER(440) and any other of the 5 amino acids\
\n3. Phosphorylation of SER(440) can alter the activity of ChAT, so important residues to focus on for hydrogen bonding are SER(440) and GLN(144).\
\nInterpretation\
\n1. This structure is bound to CoA, therefore any molecules that score high in this conformation could compete with CoA especially considering the normal cellular concentration of CoA.\
\n2. It should be noted that ChAT can be phosphorylated on as many as six residues and this appears to alter its subcellular localization and activity. This, along with the fact that ChAT has many substrates, makes predicting its activity based on docking difficult"
          )
        with doc.create(Subsection('The model for human choline acetylcholine Alpha4 Beta2 nicotinic receptor (nAChR) PDB: 5KXI asNCT')):
            doc.append("Active molecule interactions for models.\
\n1. A key interaction in the brain nAChR is a cation-pi interaction between molecules and TRP(156).\
\n2. Active molecules also potentially had a hydrogen bond to TRP(156).\
\n3. To insure inactive molecules would occupy a similar position in the pocket, molecules with hydrophobic contacts to PHE(119), TRP(156), TYR(197), TYR(204) and VAL(111) were selected.\
\nInterpretation\
\n1. The crystal structure was bound to nicotine so the receptor is in an active conformation. Molecules label active should be considered agonists."
          )
        with doc.create(Subsection('The model for human neuronal Alpha2 homopentamer nicotinic acetylcholine receptor (nAChR) PDB: 5FJV asEPJ')):
            doc.append("Active molecule interactions for models.\
\n1. The interactions described as important for binding in general are interactions with the side chains of HIS(138) and HIS(146).\
\n2. This receptor was bound to agonist epibatidine (EPJ) which forms cation-pi and hydrogen bonds to TRP(178).\
\n3. All molecules used to build the model have a hydrogen bond, cation-pi interaction or salt-bridge with one of those three at least.\
\nInterpretation\
\n1. This receptor is agonist bound so molecule predicted to be active should be considered agonists."
          )
        with doc.create(Subsection('The model for human glycine receptor-alpha 3, PDB: 5CFB asSY9')):
            doc.append("Active molecule interactions for models\
\n1. For strychnine binding PHE(63 & 207) are described as essential and in this structure strychnine forms a hydrogen bond with PHE(159).\
\n2. Active molecules were filtered for hydrogen bonds with PHE(159) and pi interactions with PHE(63 & 207) .\
\n3. We also considered hydrogen bonds with ARG(65), which is important for glycine binding and many molecules formed hydrogen bonds with it in this pocket conformation.\
\nInterpretation\
\n1. This structure was bound to antagonist strychnine, molecules labeled active should be considered antagonists."
          )
        with doc.create(Subsection('The model for human GluA2o ligand binding domain, AMPA receptor PDB: 5YBF asGLU')):
            doc.append("Active molecule interactions for models\
\n1. In the glutamate or AMPA binding domain hydrogen bonds with ARG(485) and SER(654) or corresponding ARG(506) and SER(675) (in 5ybf) are most important for domain closure which leads to activation.\
\n2. Additionally, molecule hydrogen bonding with ARG(506), GLU(726), SER(675) and THR(676) were consider actives for the model.\
\nInterpretation\
\n1. This structure was bound to glutamate and an allosteric activator, therefore molecules labeled active should be considered allosteric activators."
          )
        with doc.create(Subsection('The model for human GluK2 ligand binding domain, Kinate receptor PDB: 3QXM asNDZ')):
            doc.append("Active molecule interactions for models\
\n1. In the Kinate binding domain hydrogen bonds with ARG(485) and SER(654) or corresponding ARG(492) and ALA(658) but THR(659) was used in 3qxm since the molecule crystalized in this structure formed a hydrogen bond with it.\
\n2. Additionally, molecule hydrogen bonding with GLU(707) and THR(710) were consider actives for the model.\
\nInterpretation\
\n1. This structure was bound to an agonist, therefore molecules labeled active should be considered agonists."
          )
        with doc.create(Subsection('The model for rat glutamate receptor, AMPA sub-type GluA2 PDB: 3KGC asZK1')):
            doc.append("Active molecule interactions for models\
\n1. In the glutamate or AMPA binding domain hydrogen bonds with ARG(485) and SER(654) or corresponding ARG(96) and SER(142) (in 3kgc) are most important.\
\n2. In this case the same residues interact with the antagonist, however binding site is stabilized in the apo or open conformation and without domain closure activation does not occur.\
\nInterpretation\
\n1. This structure was about to antagonist ZK200775, molecules labeled active should be considered antagonists."
          )
        with doc.create(Subsection('The model for human GluA2o allosteric binding domain, AMPA receptor PDB: 5YBG as8SO')):
            doc.append("Active molecule interactions for models\
\n1. The article published with this structure describes residues SER(518) PRO(515) and GLY(752) as import of activity for allosteric binders, however other articles focus on SER(750) instead of GLY(752).\
\n2. To collect actives, molecules that make hydrogen bonds to each residue were used.\
\nInterpretation\
\n1. This structure was bound to glutamate and an allosteric activator, therefore molecules labeled active should be considered agonists."
          )
        with doc.create(Subsection('The model for human GluN2A receptor, NMDA receptor PDB: 5I2N as67J')):
            doc.append("Active molecule interactions for models\
\n1. All actives include hydrogen bonds to THR759 or pi interaction with TYR535, as these residues are focused on for activity.\
\n2. Note the protein structure residues had to be renumbered and this was done by aligning the protein sequence with NCBI protein blast.\
\nInterpretation\
\n1. The structure was bound to a Positive Allosteric Modulator (PAM) so molecules labeled active should be considered PAMs."
          )
        with doc.create(Subsection('The model for human GABA(A) receptor PDB: 4COF asBEN')):
            doc.append("Active molecule interactions for models\
\n1. Activating interactions were hydrogen bonds with SER156 as well as hydrogen bonds and cation-pi interactions with TYR157 and hydrogen bonds with GLU155.\
\nInterpretation\
\n1. This structure was bound to an agonist, molecules labeled active should be considered agonists."
          )
        with doc.create(Subsection('The model for human GABA(B) receptor PDB: 4MS4 as2CO')):
            doc.append("Active molecule interactions for models\
\n1. There are a few key residues to form hydrogen bonds with for both agonist and antagonist. They are SER(130), SER(153), HIS(170) and GLU(349).\
\n2. Unique to agonist are hydrogen bonds with TYR(250) and a cation-pi interaction with TRP(278).\
\nInterpretation\
\n1. This crystal structure was bound to an agonist and therefore molecules labeled active should be considered agonists."
          )
        with doc.create(Subsection('The model for human D2 Dopamine receptor PDB: 6CM4 as8NU')):
            doc.append("Active molecule interactions for models\
\n1. The article describing this crystal structure points out that the following residues, when mutated, reduced binding of the inverse agonist used by more than tenfold. So all actives had a salt-bridge or hydrogen bond to least one of them, prior the applying the criteria above.\
\n2. ASP(114), THR(119), PHE(198), PHE(382), TRP(386), PHE(389), THR(412) and TYR(416)\
\nInterpretation\
\n1. Inverse agonist makes similar interactions as agonist but occupy only the deepest part of the orthosteric pocket\
\n2. This protein structure is an inactive conformation and crystalized with an inverse agonist, therefore docked molecules labeled active should be considered potential inverse agonists."
          )
        with doc.create(Subsection('The model for human D3 Dopamine receptor PDB: 3PBL asETQ')):
            doc.append("Active molecule interactions for models\
\n1. The article describing this crystal structure points out that a salt-bridge or hydrogen bonds with ASP110 are critical for high affinity binding.\
\n2. All actives for the model had either a salt-bridge or hydrogen bond with ASP(110)\
\nInterpretation\
\n1.  Inverse agonist make similar interactions as agonist but occupy only the deepest part of the orthosteric pocket\
\n2. This protein structure is an inactive conformation and crystalized with an inverse agonist, therefore docked molecules labeled active should be considered potential inverse agonists."
          )
        with doc.create(Subsection('The model for human D4 Dopamine receptor PDB: 5WIU asAQD')):
            doc.append("Active molecule interactions for models\
\n1. Molecules for actives were chosen based on hydrogen bonds with ASP(115) and SER(196) \
\n2. From those molecules I selected a subset that had hydrophobic interactions with LEU(187)\
\nInterpretation\
\n1. This protein structure is an inactive conformation and crystalized with an antagonist, therefore docked molecules labeled active should be considered a potential antagonists."
          )
        with doc.create(Subsection('The model for human serotonin 5-HT1B G protein-coupled receptor PDB: 4IAQ as2GM')):
            doc.append("Active molecule interactions for models\
\n1. Molecules for actives were chosen based on salt-bridge or hydrogen bonds with ASP(129)\
\n2. The original publication with the crystal structure describes mutations of ASP(129) to ALA abolish agonist activity\
\n3. The publication also mentions an aromatic cage around the ligand formed by TYR(359), TRP(327), PHE(331), and PHE(330). Changing any of the four to ALA abolishes activity\
\nInterpretation\
\n1. This protein structure is an active conformation and crystalized with an agonist, therefore docked molecules labeled active should be considered potential agonists."
          )
        with doc.create(Subsection('The model for human serotonin 5-HT1B G protein-coupled receptor PDB: 5V54 as89F')):
            doc.append("Active molecule interactions for models\
\n1. Molecules for actives were chosen based on salt-bridge or hydrogen bonds with ASP(129)\
\n2. The original publication with the crystal structure describes mutations of ASP(129) to ALA abolish agonist activity\
\n3. The publication also mentions an aromatic cage around the ligand formed by TYR(359), TRP(327), PHE(331), and PHE(330). Changing any of the four to ALA abolishes activity\
\nInterpretation\
\n1. This protein structure is an inactive conformation and crystalized with an inverse agonist, therefore docked molecules labeled active should be considered potential inverse agonists.\
\n2. Note, inverse agonist occupy only the deepest part of the othersteric pocket."
         )
        with doc.create(Subsection('The model for human purinoceptor 3 (P2X3) receptor PDB: 5SVK asATP')):
            doc.append("Active molecule interactions for models\
\n1. This study focused on the differences in binding of TNP-ATP (antagonist) and ATP (agonist).\
\n2. The main difference is the area of the pocked occupied by TNP-ATP and ATP, TNP-ATP is buried deeper into the binding site and limits the movement of residue PHE(174).\
\n3. Active interactions considered were hydrogen bonds with LYS(65), ASN(279), ARG(281), and LYS(299), while filtering out molecules that interact with PHE(174).\
\nInterpretation\
\n1. This pocket is large and bound to ATP, so molecules labeled active should be considered agonists.\
\n2. However, if molecules interact with PHE(174) they should be considered inactive."
         )
        with doc.create(Subsection('The model for human acetylcholinesterase PDB: 5HF5 asDEP')):
            doc.append("Active molecule interactions for models\
\n1. The most important residue for inhbition is SER(203). Phosphorylation of SER(203) is seen with inhibitor and if dealkylation occur inhbition is irreversible\
\n2. Actives molecules for the deep-learning model had hydrogen bonds with SER(203), HIS(447), GLY(121) and GLY(122).\
\nInterpretation\
\n1. Docking can not determine covalent bonding therefore only molecules with hydrogen bonds to SER(203) and the most residues mentioned should be conisdered.\
\n2. This receptor was bound to an inhbitor so molecule predict active should be considered inhibitors."
         )
        with doc.create(Subsection('The model for human Glycine Receptor alpha-3 receptor PDB: 5TIN 7C6')):
            doc.append("Active molecule interactions for models\
\n1. This study focused on an allosteric activatory bound near the glycine binding site.\
\n2. The activator binds partialy to loop B which contains TYR(161) (important for the activator) and GLU(157), SER(158) and PHE(159) (important for glycine binding).\
\n3. Considered important for activity were pi-interactions with ASP(84) and TYR(161) as well as hydrogen bonds with ARG(29).\
\nInterpretation\
\n1. This is an allosteric activation site bound to an activator, molecules predict active should be considered allosteric activators.\
\n"
         )
            doc.append(NewPage())
def add_heatmaps(doc, file_heatmap, file_heatmap_mass, file_heatpredmap, file_heatmap_scoreden):
    with doc.create(Section('The heatmap figues')):
        with doc.create(Subsection('heatmaps')):
                with doc.create(Figure(position='ht')) as thing:
                   #doc.append(Command('centering'))
                   with doc.create(SubFigure(
                        position='t', width=NoEscape(r'0.50\textwidth'))) as left_thing:
                        left_thing.add_image(file_heatmap, width=NoEscape(r'0.95\textwidth'))
                        left_thing.add_caption('Full heatmap of all scores')
                   with doc.create(SubFigure(
                        position='t', width=NoEscape(r'0.50\textwidth'))) as right_thing:
                        right_thing.add_image(file_heatmap_mass, width=NoEscape(r'0.95\textwidth'))
                        right_thing.add_caption('Filtered heatmap of scores per mass that equate to micromolar binding affinity')
                with doc.create(Figure(position='ht')) as thing:
                   doc.append(Command('centering'))
                   with doc.create(SubFigure(
                        position='t', width=NoEscape(r'0.50\linewidth'))) as top_thing:
                        top_thing.add_image(file_heatpredmap, width=NoEscape(r'0.95\linewidth'))
                        top_thing.add_caption('Heatmap of all predictions')
                   with doc.create(SubFigure(
                        position='t', width=NoEscape(r'0.50\linewidth'))) as bottom_thing:
                        bottom_thing.add_image(file_heatmap_scoreden, width=NoEscape(r'0.95\linewidth'))
                        bottom_thing.add_caption('Heatmap of predictions filtered based on score per mass and a predictions probablity over 50%')
                   #thing.add_caption("top stuff")
    doc.append(NewPage())
    return doc


def generate_full_table(doc):
        list_of_proteins = []
        list_of_sankey = []
        list_of_logs = []
        for dir, subdirs, files in os.walk(protein_files):
            for name in files:
                protein_name = os.fsdecode(name)
                if protein_name.endswith(".pdbqt"):
                    #print(protein_name.replace(".pdbqt",""))
                    list_of_proteins.append(protein_name.replace(".pdbqt",""))

        for sankey_pics_location in os.listdir(dlg_files):
            if sankey_pics_location.endswith("sankey.png"):
                sankey_full_files = os.path.join(dlg_files, sankey_pics_location)
                list_of_sankey.append(sankey_pics_location)

        for mol_log_location in os.listdir(dlg_files):
            if mol_log_location.endswith("parsed_log.txt"):
                sankey_full_files = os.path.join(dlg_files, mol_log_location)
                #print("this is list_of_proteins[2]",list_of_proteins[2] )
                list_of_logs.append(mol_log_location)
                #print("this is list of logs ", list_of_logs)
                    #j=j+1
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
        merge_col = []
        for j in list_of_proteins:
           k = k +1
           for i, elem in enumerate(list_of_logs):
            if list_of_proteins[k] in elem:
                #print("this is j: ", k)
                file_with_path = os.path.join(dlg_files, elem)
                log_contents = pd.read_csv(file_with_path, index_col=0, sep=',', engine='python')
                file = open(file_with_path, 'r')
                contents = file.read()
                protein_name = list_of_proteins[k]
                mol_name = elem.replace("_"+list_of_proteins[k]+"_parsed_log.txt", "")
                merge_mol_prot = mol_name+"_"+protein_name
                ori_log.append(elem)
                curr_protein.append(protein_name)
                curr_mol.append(mol_name)
                merge_col.append(merge_mol_prot)
                with open(file_with_path) as search:
                    for line in search:
                       #print(line)
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
                    "Molecule Name" : curr_mol,
                    #"Original log file name" : ori_log,
                    "Hydrogen bonds" : hbonds,
                    #"Hydrophobic interactions" : hydophob,
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
        result_input = pd.read_csv('docking_summary_and_prediction_key.txt', delimiter=',', header=0, engine='python', usecols=['lowestEnergy_dlgfn','protein_id','mol_name','LE','Score_per_Mass','probability','Predciton_heatmap_format','Score_per_Prediction'])[['lowestEnergy_dlgfn','protein_id','mol_name','LE','Score_per_Mass','probability','Predciton_heatmap_format','Score_per_Prediction']]

        with doc.create(Section('Table')):
            with doc.create(LongTabu("X[l] X[l] X[l] X[l] X[l] X[l]")) as table:
                 table.add_hline()
                 interact_df = pd.read_csv('table_of_interaction_info.csv', delimiter=',', header=0, engine='python')
                 print(interact_df.head())
                 new_input1 = pd.merge(left=interact_df, right=result_input, how='right', left_on='merge_column', right_on='lowestEnergy_dlgfn')
                 print("This is your merged dataframe: ", new_input1.head())
                 Class_prediction = []
                 for b in new_input1['Predciton_heatmap_format']:
                    #print(b)
                    i = round(b,3)
                    #print(i)
                    if 2.000 <= i <= 3.000:
                        #print("between 2 and 3")
                        Class_prediction.append('Active, high-interaction')
                    elif 1.000 <= i <= 2.000:
                        #print("between 1 and 2")
                        Class_prediction.append('Active, med-interaction')
                    elif 0.1 <= i <= 1.000:
                        #print("between 0 and 1")
                        Class_prediction.append('Active, low-interaction')
                    elif b == 0.0:
                        #print("right at 0.0")
                        Class_prediction.append('no prediction')
                    elif b >= 3.0:
                        #print("better than 4")
                        Class_prediction.append('Inactive')
                 new_input1['Class prediction'] = Class_prediction
                 new_input = new_input1.round(decimals=4)
                 result_col_names = ['mol_name','LE','Score_per_Mass','Class prediction','probability']
                 new_title_results = ['Molecule Name', 'AD4 dG', 'Score per Mass', 'Activity prediction', 'Probablity']
                 column_names = ['Molecule Name','Hydrogen bonds','pi-pi stacking interactions','T-stacking','cation-pi','Salt Bridges']
                 table.add_hline()
                 table.add_hline()
                 new_element = []
                 ####I have changed new_input to report_select in an effort to make it match the report printed... we'll see if this works
                 for element in new_input['Protein Name']:
                     if element not in new_element:
                         new_element.append(element)
                         table.add_hline()
                         table.add_hline()
                         table.add_row(new_title_results,strict=False)
                         table.add_hline()
                         table.add_row(column_names)
                         table.add_empty_row()
                         table.add_hline()
                         table.add_hline()
                         table.add_row((MultiColumn(6, align='c', data="Protein Name: "+element),))
                         table.add_hline()
                         table.add_hline()
                         newnew_input = new_input[new_input['Protein Name'] == element]
                         newresult_input = result_input[result_input['protein_id'] == element]
                         for row in newnew_input.index:
                             result_lookup = new_input.loc[row, 'Report Molecule Name']

                             table.add_row(list(new_input.loc[row, result_col_names]),strict=False)
                             table.add_row(list(new_input.loc[row, column_names]), strict=False)
                             table.add_hline()



def generate_custom_table(report_select, doc):
        list_of_proteins = []
        list_of_sankey = []
        list_of_logs = []
        for dir, subdirs, files in os.walk(protein_files):
            for name in files:
                protein_name = os.fsdecode(name)
                if protein_name.endswith(".pdbqt"):
                   list_of_proteins.append(protein_name.replace(".pdbqt",""))

        for sankey_pics_location in os.listdir(dlg_files):
            if sankey_pics_location.endswith("sankey.png"):
                sankey_full_files = os.path.join(dlg_files, sankey_pics_location)
                list_of_sankey.append(sankey_pics_location)

        for mol_log_location in os.listdir(dlg_files):
            if mol_log_location.endswith("parsed_log.txt"):
                sankey_full_files = os.path.join(dlg_files, mol_log_location)
                list_of_logs.append(mol_log_location)

        result_input = pd.read_csv('docking_summary_and_prediction_key.txt', delimiter=',', header=0, engine='python', usecols=['lowestEnergy_dlgfn','protein_id','mol_name','LE','Score_per_Mass','probability','Predciton_heatmap_format','Score_per_Prediction'])[['lowestEnergy_dlgfn','protein_id','mol_name','LE','Score_per_Mass','probability','Predciton_heatmap_format','Score_per_Prediction']]


        with doc.create(Section('Table')):
            with doc.create(LongTabu("X[l] X[l] X[l] X[l] X[l] X[l]")) as table:
                 table.add_hline()
                 interact_df = pd.read_csv('table_of_interaction_info.csv', delimiter=',', header=0, engine='python')
                 new_input1 = pd.merge(left=interact_df, right=result_input, how='right', left_on='merge_column', right_on='lowestEnergy_dlgfn')
                 Class_prediction = []
                 for b in new_input1['Predciton_heatmap_format']:
                    i = round(b,3)
                    if 2.000 <= i <= 3.000:
                        Class_prediction.append('Active, high-interaction')
                    elif 1.000 <= i <= 2.000:
                        Class_prediction.append('Active, med-interaction')
                    elif 0.1 <= i <= 1.000:
                        Class_prediction.append('Active, low-interaction')
                    elif b == 0.0:
                        Class_prediction.append('no prediction')
                    elif b >= 3.0:
                        Class_prediction.append('Inactive')
                 new_input1['Class prediction'] = Class_prediction
                 new_input = new_input1.round(decimals=4)
                 result_col_names = ['Molecule Name','LE','Score_per_Mass','Class prediction','probability']
                 new_title_results = ['Molecule Name', 'AD4 dG', 'Score per Mass', 'Activity prediction', 'Probablity']
                 column_names = ['Molecule Name','Hydrogen bonds','pi-pi stacking interactions','T-stacking','cation-pi','Salt Bridges']
                 table.add_hline()
                 table.add_hline()
                 all_elements = []
                 all_molecules = []
                 newnew_input = pd.merge(new_input, report_select, on=['Protein Name','Report Molecule Name'], how='inner')
                 for element in newnew_input['Protein Name']:
                     molecules = report_select['Report Molecule Name']
                     protein = newnew_input['Protein Name']
                     if element not in all_elements:
                         all_elements.append(element)
                         table.add_hline()
                         table.add_hline()
                         table.add_row(new_title_results,strict=False)
                         table.add_hline()
                         table.add_row(column_names)
                         table.add_empty_row()
                         table.add_hline()
                         table.add_hline()
                         table.add_row((MultiColumn(6, align='c', data="Protein Name: "+ element),))
                         table.add_hline()
                         table.add_hline()
                         table.add_hline()
                         print("this is newnew_input", newnew_input)
                         for row in newnew_input[newnew_input['Protein Name']== element].index:
                             table.add_row(list(newnew_input.loc[row, result_col_names]),strict=False)
                             table.add_row(list(newnew_input.loc[row, column_names]), strict=False)
                             table.add_empty_row()
                             table.add_hline()
        #doc.append(NewPage())

def add_sankey_figures(*fig, doc):
        print(len(fig))
        list_len = len(fig)
        if list_len >= 8:
            with doc.create(Section('The figues')):
              with doc.create(Subsection('Sankey diagrams')):
                with doc.create(Figure(position='ht')) as thing:
                   #doc.append(Command('centering'))
                   with doc.create(SubFigure(
                        position='t', width=NoEscape(r'0.50\textwidth'))) as left_thing:
                        left_thing.add_image(fig[0], width=NoEscape(r'0.95\textwidth'))
                        left_thing.add_caption(fig[4])
                   with doc.create(SubFigure(
                        position='t', width=NoEscape(r'0.50\textwidth'))) as right_thing:
                        right_thing.add_image(fig[1], width=NoEscape(r'0.95\textwidth'))
                        right_thing.add_caption(fig[5])
                with doc.create(Figure(position='ht')) as thing:
                   doc.append(Command('centering'))
                   with doc.create(SubFigure(
                        position='t', width=NoEscape(r'0.50\linewidth'))) as top_thing:
                        top_thing.add_image(fig[2], width=NoEscape(r'0.95\linewidth'))
                        top_thing.add_caption(fig[6])
                   with doc.create(SubFigure(
                        position='t', width=NoEscape(r'0.50\linewidth'))) as bottom_thing:
                        bottom_thing.add_image(fig[3], width=NoEscape(r'0.95\linewidth'))
                        bottom_thing.add_caption(fig[7])
            doc.append(NewPage())
        elif 6 <= list_len:
            with doc.create(Section('The figues')):
              with doc.create(Subsection('Sankey diagrams')):
                with doc.create(Figure(position='ht')) as thing:
                   with doc.create(SubFigure(
                        position='t', width=NoEscape(r'0.50\textwidth'))) as left_thing:
                        left_thing.add_image(fig[0], width=NoEscape(r'0.95\textwidth'))
                        left_thing.add_caption(fig[3])
                   with doc.create(SubFigure(
                        position='t', width=NoEscape(r'0.50\textwidth'))) as right_thing:
                        right_thing.add_image(fig[1], width=NoEscape(r'0.95\textwidth'))
                        right_thing.add_caption(fig[4])
                with doc.create(Figure(position='ht')) as thing:
                   doc.append(Command('centering'))
                   with doc.create(SubFigure(
                        position='t', width=NoEscape(r'0.50\linewidth'))) as top_thing:
                        top_thing.add_image(fig[2], width=NoEscape(r'0.95\linewidth'))
                        top_thing.add_caption(fig[5])
            doc.append(NewPage())
        elif 4 <= list_len:
            with doc.create(Section('The figues')):
              with doc.create(Subsection('Sankey diagrams')):
                with doc.create(Figure(position='ht')) as thing:
                   #doc.append(Command('centering'))
                   with doc.create(SubFigure(
                        position='t', width=NoEscape(r'0.50\textwidth'))) as left_thing:
                        left_thing.add_image(fig[0], width=NoEscape(r'0.95\textwidth'))
                        left_thing.add_caption(fig[2])
                   with doc.create(SubFigure(
                        position='t', width=NoEscape(r'0.50\textwidth'))) as right_thing:
                        right_thing.add_image(fig[1], width=NoEscape(r'0.95\textwidth'))
                        right_thing.add_caption(fig[3])
            doc.append(NewPage())
        elif 2 <= list_len:
            with doc.create(Section('The figues')):
              with doc.create(Subsection('Sankey diagrams')):
                with doc.create(Figure(position='ht')) as thing:
                   #doc.append(Command('centering'))
                   with doc.create(SubFigure(
                        position='t', width=NoEscape(r'0.50\textwidth'))) as left_thing:
                        left_thing.add_image(fig[0], width=NoEscape(r'0.95\textwidth'))
                        left_thing.add_caption(fig[1])
            doc.append(NewPage())
        else:
            pass

        x=0
        y=1
        return
def define_dlg_protein_log(dlg_temp, protein_temp):
    global dlg_files
    dlg_files = dlg_temp
    global protein_files
    protein_files = protein_temp
    return


def generate_report(report_select, list_len, file_heatmap, file_heatmap_mass, file_heatpredmap,file_heatmap_scoreden, dlg_temp, protein_temp, fileName, working_project):
        x=0
        y=1
        new_list = []
        fig_one = 0
        fig_two = 1
        fig_three = 2
        fig_four = 3
        define_dlg_protein_log(dlg_temp, protein_temp)
        geometry_options = {"tmargin": "1.5cm", "lmargin": "1.0cm", "rmargin" : "1.0cm"}
        doc = Document(geometry_options=geometry_options)
        generate_text_report(doc)
        add_heatmaps(doc, file_heatmap, file_heatmap_mass, file_heatpredmap,file_heatmap_scoreden)
        for fig in report_select.index:
            #figure_next = figure2 + 1
            if fig_one < list_len and fig_two < list_len and fig_three < list_len and fig_four < list_len:
                fig_1, cap_1 = report_select["File location"][fig_one],report_select["Caption"][fig_one]
                fig_2, cap_2 = report_select["File location"][fig_two], report_select["Caption"][fig_two]
                fig_3, cap_3 = report_select["File location"][fig_three], report_select["Caption"][fig_three]
                fig_4, cap_4 = report_select["File location"][fig_four], report_select["Caption"][fig_four]
                fig_one = fig_one + 4
                fig_two = fig_two + 4
                fig_three = fig_three + 4
                fig_four = fig_four + 4
                add_sankey_figures(fig_1, fig_2, fig_3, fig_4, cap_1,cap_2,cap_3,cap_4, doc=doc )
            elif fig_one < list_len and fig_two < list_len and fig_three < list_len and not fig_four < list_len:
                fig_1, cap_1 = report_select["File location"][fig_one],report_select["Caption"][fig_one]
                fig_2, cap_2 = report_select["File location"][fig_two], report_select["Caption"][fig_two]
                fig_3, cap_3 = report_select["File location"][fig_three], report_select["Caption"][fig_three]
                fig_one = fig_one + 4
                fig_two = fig_two + 4
                fig_three = fig_three + 4
                add_sankey_figures(fig_1, fig_2, fig_3, cap_1, cap_2, cap_3, doc=doc)
            elif fig_one < list_len and fig_two < list_len and not fig_three < list_len and not fig_four < list_len:
                fig_1, cap_1 = report_select["File location"][fig_one],report_select["Caption"][fig_one]
                fig_2, cap_2 = report_select["File location"][fig_two], report_select["Caption"][fig_two]
                fig_one = fig_one + 4
                fig_two = fig_two + 4
                add_sankey_figures(fig_1, fig_2, cap_1, cap_2, doc=doc)
            elif fig_one < list_len and not fig_two < list_len and not fig_three < list_len and not fig_four < list_len:
                fig_1, cap_1 = report_select["File location"][fig_one],report_select["Caption"][fig_one]
                fig_one = fig_one + 4
                add_sankey_figures(fig_1, cap_1, doc=doc)
            else:
                pass
        #add_heatmaps(doc, file_heatmap, file_heatmap_mass, file_heatpredmap,file_heatmap_scoreden)
        generate_custom_table(report_select, doc)
        ### I need to put this back in to generate the report
        doc.generate_pdf(fileName, clean_tex=False)
        if fileName not in os.listdir(working_project):
            print("the working directory is: ", working_project)
            print(fileName)
            try:
                shutil.move(fileName, working_project)
            except:
                pass
        return
