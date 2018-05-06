#! /usr/bin/env python

import os,sys
import PMPformat as PMP
from rdkit import Chem

### Function to extract data to pmp

def ProcessGDMAdata(fname, multipole=2):
    ''' Process GDMA output 
    Return AtomPropDict{pname:[values..]}, MolPropDict{pname: value}'''
    with open(fname) as f:
        lines = f.readlines()
        atomdatas={}
        moldatas = {}
        for i in range(multipole):
            if i == 0:
                atomdatas['Dipole']=[]
            elif i == 1:
                atomdatas['Quadrupole']=[]
            elif i == 2:
                atomdatas['Octopole']=[]           
            elif i == 3:
                atomdatas['Hexadecapole']=[] 
        for line in lines:
            if line[0] != "#" and len(line)>50:
                data=line.split()
                for i in range(multipole):
                    if i == 0:
                        atomdatas['Dipole'].append(float(data[6]))
                    elif i == 1:
                        atomdatas['Quadrupole'].append(float(data[7]))
                    elif i == 2:
                        atomdatas['Octopole'].append(float(data[8]))           
                    elif i == 3:
                        atomdatas['Hexadecapole'].append(float(data[9]))
            elif line[0] != "#":
                data=line.split()
                for i in range(multipole):
                    if i == 0:
                        moldatas['Dipole'] = float(data[1])
                    elif i == 1:
                        moldatas['Quadrupole'] = float(data[2])
                    elif i == 2:
                        moldatas['Octopole'] = float(data[3])           
                    elif i == 3:
                        moldatas['Hexadecapole'] = float(data[4])
        return atomdatas, moldatas


def ProcessAtomSolEng(fname):
    '''Read Atomic Solvation Energy data from MIBPB result.'''
    with open(fname) as f:
        lines = f.readlines()
        datas=map(float,[ line.strip() for line in lines])
    return datas

def ProcessPQRTA(fname, atomtype, charge="resp", radius="mbondi"):
    '''Read many data(Atomic data) and (Mol data) from PQRTA'''
    moldatas={}
    atomdatas={charge:[],radius:[],atomtype:[],"AtomArea":[]}
    with open(fname) as f:
        for line in f:
            if line[:6] == "REMARK":
                tmps=line.strip().split()
                if tmps[1] == "AREAS":
                    moldatas['MolArea']=float(tmps[2])
                elif tmps[1] == "VOLUMES":
                    moldatas['MolVolume']=float(tmps[2])
                elif tmps[1] == "AREA":
                    moldatas['Area_'+tmps[2]]=float(tmps[3])
            elif line[:6] == "ATOM  " or line[:6] == "HETATM":
                atomdatas[charge].append(float(line[54:62].strip()))
                atomdatas[radius].append(float(line[62:70].strip()))
                atomdatas[atomtype].append(line[80:88].strip())
                atomdatas["AtomArea"].append(float(line[88:100].strip()))
        return atomdatas, moldatas

def ProcessMol2AtomType(fname):
    '''Read SYBYL Atom Type from mol2 file'''
    datas=[]
    findatom=False
    with open(fname) as f:   
        for line in f:
            if not findatom and "@<TRIPOS>ATOM" in line:
                findatom=True
                continue
            if "@<TRIPOS>BOND" in line:
                break
            if findatom:
                datas.append(line.split()[5])
    return datas

def ProcessMIBPBresults(fname):
    '''Read Electrostatic solvation energy from MIBPB output'''
    with open(fname) as f:   
        for line in f:
            if "Electrostatics solvation engergy=:" in line:
                return float(line.strip().split(':')[1].strip())


## Using PQRTA file as input
fname = sys.argv[1]
bename = os.path.basename(fname)
bname = os.path.splitext(fname)
dname = os.path.dirname(fname)
if not dname:
    dname='.'
deleney=bname[0][:-5]
deleneynum=int(bename[8:-11])

pmpf = PMP.PMPFormator()
mol = pmpf.MolFromPDBFile(fname)

if len(sys.argv) == 2:
    mol = pmpf.MolMatchBondByMol2File(deleney+'.mol2')
    smiles = Chem.MolToSmiles(Chem.RemoveHs(mol,updateExplicitCount=True))
    ReallogS = -9999.999
elif len(sys.argv) > 2:
    import csv
    linecount = 0
    with open(sys.argv[2],'rb') as csvfile:
        csvread = csv.reader(csvfile)
        ReallogS = -9999.999
        for row in csvread:
            if deleneynum == linecount:
                smiles = row[9]
                ReallogS = float(row[8])
                break
            linecount += 1
        mol = pmpf.MolMatchBondBySmiles(smiles)
        
pmpf.SetMolProp(smiles, 'SMILES', ptype='s', plen=10)
if ReallogS != -9999.999:
    pmpf.SetMolProp(ReallogS, 'ExpLogS', ptype='f', plen=10, floatPoint=3)

atomtype="AT_gaff"
AtomDatas, MolDatas= ProcessPQRTA(bname[0]+'.pqrta', atomtype=atomtype,charge="resp", radius="mbondi")
for item,value in AtomDatas.items():
    if item != atomtype:
        pmpf.SetAtomsProp(value, item, ptype='f', plen=10,floatPoint=4)
    else:
        pmpf.SetAtomsProp(value, item, ptype='s', plen=6)            
for item,value in MolDatas.items():
    pmpf.SetMolProp(value, item, ptype='f', plen=10,floatPoint=4) 

sybyltype = ProcessMol2AtomType(deleney+'.mol2')
pmpf.SetAtomsProp(sybyltype, "AT_sybyl", ptype='s', plen=6)

mibpbOut = ProcessMIBPBresults(dname+'/mibpb5.log')
pmpf.SetMolProp(mibpbOut, 'ElecSolvEng', ptype='f', plen=10,floatPoint=6)

GDMAdata, MolGDMA = ProcessGDMAdata(deleney+'.gdma',multipole=2)
for item,value in GDMAdata.items():
    pmpf.SetAtomsProp(value, item, ptype='f', plen=10,floatPoint=6)
for item,value in MolGDMA.items():
    pmpf.SetMolProp(value, item, ptype='f', plen=10,floatPoint=6)

soldata = ProcessAtomSolEng(dname+'/AtomSoleng.txt')
pmpf.SetAtomsProp(soldata, "AtomSolEng", ptype='f', plen=10,floatPoint=6)

## Write the pmp file
pmpf.MolToPMPFile(deleney+'.pmp')

if False:
    print pmpf.GetPMP()
    print pmpf._mol.GetPropsAsDict()
    print pmpf._mol.GetAtomWithIdx(0).GetPropsAsDict()

    
    pmpf2 = PMP.PMPFormator(pmpfile=deleney+'.pmp')
    print Chem.MolToSmiles(pmpf2._mol)
    print pmpf2._mol.GetPropsAsDict()
    print pmpf2._mol.GetAtomWithIdx(0).GetPropsAsDict()
