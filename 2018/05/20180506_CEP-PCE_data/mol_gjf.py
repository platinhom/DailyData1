#! /usr/bin/env python

## Parameters
p_nproc=1
p_mem=1
p_basis="6-31G*"
p_OtherMethod="3-21G*"#"LanL2DZ"
p_method="HF/"+p_basis+" opt SCF=tight Test Pop=MK iop(6/33=2) iop(6/42=6) iop(6/50=1)"
excludeElement=['I']
noradiiElement=['Se','I']

import sys,os,StringIO

if not os.path.exists(sys.argv[1]):
	print "Error: No input file exist"
	sys.exit()

f=open(sys.argv[1])
fntmp=os.path.splitext(sys.argv[1])

fofile=open(fntmp[0]+".gjf",'w')

fo = StringIO.StringIO()
findatom=False
elements=set()
totCharge=0
for line in f:
	if 'BEGIN ATOM' in line:
		findatom=True
		continue
	if 'END ATOM' in line:
		break
	if findatom:
		tmp=line.split()
		element=tmp[3]
		elements.add(element)
		fo.write(' '+element+'        '+tmp[4]+'  '+tmp[5]+'  '+tmp[6]+'\n')
		if "CHG=" in line:
			totCharge+=int(line.split('CHG=')[1].split()[0])
fo.write('\n')

if (elements.intersection(set(excludeElement))):
	## Write normal elements
	useStr=""
	for elem in elements:
		if elem not in excludeElement:
			useStr+=elem+" "
	useStr+="0\n"+p_basis+"\n****\n"
	for elem in excludeElement:
		useStr+=elem+" "
	useStr+="0\n"+p_OtherMethod+"\n****\n"

	useStr+="\n"
	for elem in excludeElement:
		useStr+=elem+" 0\n"+p_OtherMethod+"\n"
	fo.write(useStr)
	fo.write('\n')



fofile.write('%nproc='+str(p_nproc)+'\n%mem='+str(p_mem)+'GB\n%nosave\n')
fofile.write('%chk='+fntmp[0]+'.chk\n')

if (elements.intersection(set(excludeElement))):
	fofile.write('# HF/gen opt pseudo=read SCF=tight Test Pop=(readradii,MK) iop(6/33=2) iop(6/42=6) iop(6/50=1)\n\n')
elif (elements.intersection(set(noradiiElement))):
	fofile.write('# HF/'+p_basis+' opt SCF=tight Test Pop=(readradii,MK) iop(6/33=2) iop(6/42=6) iop(6/50=1)\n\n')
else:
	fofile.write('# '+p_method+'\n\n')

fofile.write(os.path.basename(sys.argv[1])+'\n')
fofile.write('\n'+str(totCharge)+' 1\n')
fofile.write(fo.getvalue())

if (elements.intersection(set(noradiiElement))):
	if "I" in elements:
		fofile.write("I 2.0\n")
	if "Se" in elements:
		fofile.write("Se 1.9\n")
	fofile.write('\n')

fofile.write("ligand.gesp\n\n")

if (elements.intersection(set(noradiiElement))):
	if "I" in elements:
		fofile.write("I 2.0\n")
	if "Se" in elements:
		fofile.write("Se 1.9\n")
	fofile.write('\n')

fofile.write("ligand.gesp\n\n")

fofile.close()
fo.close()
f.close()
