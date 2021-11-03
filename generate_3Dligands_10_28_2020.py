

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

import sys
import time
import os
from os import rename
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import PropertyMol
from rdkit.Chem import rdDistGeom
from rdkit.Chem import Descriptors
from rdkit.Chem import RDConfig
from rdkit.Chem import PandasTools
import pandas as pd
from concurrent import futures
from rdkit import RDConfig
from rdkit.Chem import Draw
from openbabel import openbabel
from openbabel import pybel


def generate2dimage(sdf_output_file, nam):
        print("this is sdf_output_file from gen2image: ", sdf_output_file)
        print("this is nam from gen2imgae; ", nam)
        TwoDm2 = Chem.SDMolSupplier(sdf_output_file)
        for m in TwoDm2:tmp=AllChem.Compute2DCoords(m)
        print("this is type(m): ", type(m))
        pic_file_name = nam+".png"
        Draw.MolToFile(m, pic_file_name)
        return


def generateconformations(smilesuppl, sdf_output_file, nam, writer):
        m2 = Chem.AddHs(smilesuppl)
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
        m2.SetProp('image', str(nam)+".png")
        m2.SetProp('_Name', str(nam))
        for id in ids:
            writer.write(m2, confId=id)
            generate2dimage(sdf_output_file, nam)
        #writer.close()
        return

def update_sdf_file(sent_sdf_file):
    suppl = Chem.SDMolSupplier(sent_sdf_file)
    sdf_output_file = sent_sdf_file.replace(".sdf","_temp.sdf")
    writer = Chem.SDWriter(sdf_output_file)
    nam_list = []
    for mol in suppl:
        print("this is prop: ", mol.GetProp('_Name'))
        name = mol.GetProp('_Name').replace("-","_")
        mwt = Descriptors.MolWt(mol)
        tpsa = Descriptors.TPSA(mol)
        logp = Descriptors.MolLogP(mol)
        rotbond = Descriptors.NumRotatableBonds(mol)
        hbdaccept = Descriptors.NumHAcceptors(mol)
        hbddon = Descriptors.NumHDonors(mol)
        fractCSP3 = Descriptors.FractionCSP3(mol)
        nam_list.append(name+".png")
        prop = Chem.PropertyMol.PropertyMol(mol)
        prop.SetProp('image', str(name)+".png")
        prop.SetProp('MolWt', str(mwt))
        prop.SetProp('tpsa', str(tpsa))
        prop.SetProp('logp', str(logp))
        prop.SetProp('rotatable_bonds', str(rotbond))
        prop.SetProp('hbond_accept', str(hbdaccept))
        prop.SetProp('hbond_donor', str(hbddon))
        prop.SetProp('fractionCSP3', str(fractCSP3))
        prop.SetProp('_Name', str(name))
        writer.write(prop)
        ##This method works better to make pics from an existing SDF file
        tmp=AllChem.Compute2DCoords(mol)
        pic_file_name = name+".png"
        Draw.MolToFile(mol, pic_file_name)
        #generate2dimage(sdf_output_file, name)


    writer.close()
    print(str(sdf_output_file)+" "+str(sent_sdf_file))
    os.rename(str(sdf_output_file), str(sent_sdf_file))
    sdf_output_file = sent_sdf_file
    generate_output_table(sdf_output_file)
    return

def generate_output_table(sdf_output_file):

    frame = PandasTools.LoadSDF(sdf_output_file,
                                     smilesName='SMILES',
                                     includeFingerprints=False)
    frame.set_index(['ID'], inplace=True)
    frame.to_csv('excel_molecules.csv', sep=':')
    return


def gen_rdkit_mols(smi_input_file):
    sdf_output_file = smi_input_file.replace("txt","sdf")
    writer = Chem.SDWriter(sdf_output_file)
    getmols = pd.read_csv(smi_input_file, header=None)
    mols_to_write = []
    i = 0
    for index, row in getmols.iterrows():
        i=i+1
        if Chem.MolFromSmiles(row[0]):
           smilesuppl = Chem.MolFromSmiles(row[0])
           nam = "mol_smi_"+str(i)
           generateconformations(smilesuppl, sdf_output_file, nam, writer)
        elif Chem.inchi.MolFromInchi(row[0]):
           inchisuppl = Chem.inchi.MolFromInchi(row[0])
           nam = "mol_inchi_"+str(i)
           generateconformations(inchisuppl, sdf_output_file, nam, writer)
    writer.close()
    generate_output_table(sdf_output_file)
    return


def gen_obabel_mols(fileName,ffield,steps):
    mol_smi = pybel.readfile("smi", fileName)
    mol_inchi = pybel.readfile("InChI", fileName)
    output_file = fileName.replace("txt","sdf")
    output = pybel.Outputfile("sdf", output_file, overwrite=True)

    #### This SMILES input from openbabel cannot handle commas in the text files.
    i = 0
    if mol_smi:
      for mymol in mol_smi:
         mymol_smi = mymol
         i=i+1
         mymol_smi.addh()
         mymol_smi.make3D(forcefield=ffield, steps=steps)
         mymol_smi.localopt(forcefield=ffield, steps=steps)
         #descvalues = mymol.calcdesc()
         descvalues = mymol.calcdesc(descnames=['MW', 'logP', 'TPSA', 'HBA1', 'HBA2', 'HBD', 'dbonds' ])
         mymol_smi.data.update(descvalues)
         mymol_smi.title = ("mol_smi_"+str(i))
         mymol_smi.data['title'] = ("mol_smi_"+str(i))
         mymol_smi.data['image'] = ("mol_smi_"+str(i)+".png")
         mymol_smi.draw(show=False, filename="mol_smi_"+str(i)+".png", update=False, usecoords=False)
         output.write(mymol_smi)

    if mol_inchi:
      for mymol in mol_inchi:
         mymol_inchi = mymol
         i=i+1
         mymol_inchi.addh()
         mymol_inchi.make3D(forcefield=ffield, steps=steps)
         mymol_inchi.localopt(forcefield=ffield, steps=steps)
         mymol_inchi.title = ("mol_inchi_"+str(i))
         #descvalues = mymol_inchi.calcdesc()
         descvalues = mymol.calcdesc(descnames=['MW', 'logP', 'TPSA', 'HBA1', 'HBA2', 'HBD', 'dbonds' ])
         mymol_inchi.data.update(descvalues)
         mymol_inchi.data['title'] = ("mol_inchi_"+str(i))
         mymol_inchi.data['image'] = ("mol_inchi_"+str(i)+".png")
         mymol_inchi.draw(show=False, filename="mol_inchi_"+str(i)+".png", update=False, usecoords=False)
         output.write(mymol_inchi)

    output.close()
    generate_output_table(output_file)
    return

    #frame = PandasTools.LoadSDF(sdf_output_file,
    #                                 smilesName='SMILES',
    #                                 includeFingerprints=False)
    #frame.to_csv('excel_molecules.csv', sep=':')
#generate2dimage(sdf_output_file)
