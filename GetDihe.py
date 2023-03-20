from math import sqrt, cos, sin, acos, pi
import pandas as pd 
import os
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colors
from math import ceil, floor
from rama_config import RAMA_PREFERENCES
from rama_config import ROTA_PREFERENCES

pdb_columns = ['name', 'resn', 'resid', 'Chain', 'X', 'Y', 'Z', 'Element']
# Read in PDB file one line at a time, if the first four letter ar ATOM or HETA then it will parse the data into the 
# data frame, using the atome index int he PDB as the row index in the data frame. 
AAA_dict = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": 'V' }
A_dict = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS', 'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN','G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP','A': 'ALA', 'V':'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}
Methyl_groups = {"LEU" : ['CD1','CD2'], "ILE" : ['CD1'], "VAL" : ['CG1','CG2'], "THR" : ['CG2'], "MET" : ['CE'], "ALA" : ['CB']}
Methyls = ['ILE', 'LEU','VAL','MET','ALA', 'THR']
Aromatics = ['PHE', 'TYR']
Aromatic_groups = {"PHE" : ['CE1', 'CE2'], "TYR": ['CE1','CHE2']}
Pass_Atoms = ['N   ALA', 'H   ALA', 'N   ARG', 'H   ARG', 'N   ASN', 'H   ASN', 'N   ASP', 'H   ASP', 'N   CYS', 'H   CYS', 'N   GLU', 'H   GLU', 'N   GLN', 'H   GLN', 'N   GLY', 'H   GLY', 'N   HIS', 'H   HIS', 'N   ILE', 'H   ILE', 'N   LEU', 'H   LEU', 'N   LYS', 'H   LYS', 'N   MET', 'H   MET', 'N   PHE', 'H   PHE', 'N   SER', 'H   SER', 'N   THR', 'H   THR', 'N   TRP', 'H   TRP', 'N   TYR', 'H   TYR', 'N   VAL', 'H   VAL']


if len(sys.argv) == 1:
	print('Usage: GetDieh pdb chain')
	exit()

RAMA_PREF_VALUES = None
ROTA_PREF_VALUES = None

def _cache_RAMA_PREF_VALUES():
	## Values in data span -179 to 179 in incrments of 2, so values are translated to 0 to 360 in increments of 1
	f_path = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
	RAMA_PREF_VALUES = {}
	for key, val in RAMA_PREFERENCES.items():
		RAMA_PREF_VALUES[key] = np.full((360, 360), 0, dtype=np.float64)
		with open(os.path.join(f_path, val["file"])) as fn:
			for line in fn:
				if line.startswith("#"):
					continue
				else:
					x = int(float(line.split()[1]))
					y = int(float(line.split()[0]))
					RAMA_PREF_VALUES[key][x + 180][y + 180] \
						= RAMA_PREF_VALUES[key][x + 179][y + 179] \
						= RAMA_PREF_VALUES[key][x + 179][y + 180] \
						= RAMA_PREF_VALUES[key][x + 180][y + 179] \
						= float(line.split()[2]) 
	return RAMA_PREF_VALUES

def _cache_ROTA_PREF_VALUES():
	f_path = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
	ROTA_PREF_VALUES = {}
	for key, val in ROTA_PREFERENCES.items():
		ROTA_PREF_VALUES[key] = np.full((360, 360), 0, dtype=np.float64)
		with open(os.path.join(f_path, val["file"])) as fn:
			for line in fn:
				if line.startswith("#"):
					continue
				else:
					x = int(float(line.split()[1]))
					y = int(float(line.split()[0]))
					ROTA_PREF_VALUES[key][x][y] = float(line.split()[2])
	return ROTA_PREF_VALUES

def crossProduct(u,v):
	## Calculates the cross product of two 3d vectors (as 1-d arrays).
	prod = [0.,0.,0.]
	prod[0] = u[1]*v[2]-u[2]*v[1]
	prod[1] = u[2]*v[0]-u[0]*v[2]
	prod[2] = u[0]*v[1]-u[1]*v[0]
	return prod

def norm(v):
	mag = sqrt(v[0]**2 + v[1]**2 + v[2]**2)
	normv = [v[i]/mag for i in range(3)]
	return normv

def getOrthoNorm(v,ref):
	temp = [v[i] - ref[i]* dotProduct(v,ref) for i in range(3)]
	orth = norm(temp)
	return orth

def dotProduct(u,v):
	## Calculates the dot product between two vectors.
	return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]

def calcDihedrals(A,B,C,D):
	##Calculates dihedral angle side or main chain
	U2 = norm([C[i] - B[i] for i in range(3)])
	U1 = getOrthoNorm([A[i] - B[i] for i in range(3)], U2)
	U3 = getOrthoNorm([D[i] - C[i] for i in range(3)], U2)
	angle = 180/pi * acos(dotProduct(U1,U3))
	sign = dotProduct(crossProduct(U1,U3),U2)
	if sign < 0:
		angle = -angle
	return angle

# phi main chain torsion angle for atoms C-1,N,CA,C
# psi main chain torsion angle for atoms N,CA,C,N+1
# chi1 side chain torsion angle for atoms N,CA,CB,*G
# chi2 side chain torsion angle for atoms CA,CB,*G,*D
# chi3 side chain torsion angle for atoms CB,*G,*D,*E
# chi4 side chain torsion angle for atoms *G,*D,*E,*Z
# chi5 side chain torsion angle for atoms *D,*E,*Z, NH1
# return phi, psi
#MODEL        1
#
#SideDihe = {
# 'R':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD' ], ['chi3', 'CB', 'CG', 'CD', 'NE'], ['chi4', 'CG', 'CD', 'NE', 'CZ'], ['chi5', 'CD', 'NE', 'CZ' ,'NH1']],
# 'N':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'OD1']],
# 'D':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'OD1']],
# 'C':[['chi1', 'N', 'CA', 'CB','SG']],
# 'Q':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD'], ['chi3','CB', 'CG', 'CD', 'OE1']],
# 'E':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD'], ['chi3', 'CB', 'CG', 'CD', 'OE1']],
# 'H':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2','CA', 'CB', 'CG', 'ND1']],
# 'I':[['chi1', 'N', 'CA', 'CB','CG1'], ['chi2', 'CA', 'CB', 'CG1', 'CD1']],
# 'L':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD1']],
# 'K':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD'],['chi3', 'CB', 'CG', 'CD', 'CE'],['chi4', 'CG', 'CD', 'CE', 'NZ']],
# 'M':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'SD'],['chi3', 'CB', 'CG', 'SD', 'CE']],
# 'F':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD1']],
# 'P':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD']],
# 'S':[['chi1', 'N', 'CA', 'CB','OG']],
# 'T':[['chi1', 'N', 'CA', 'CB','OG1']],
# 'W':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD1']],
# 'Y':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD1']],
# 'V':[['chi1', 'N', 'CA', 'CB','CG1']]}
SideDihe = {
 'I':[['chi1', 'N', 'CA', 'CB','CG1'], ['chi2', 'CA', 'CB', 'CG1', 'CD1']],
 'L':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD1']],
 'M':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'SD']],
 'F':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD1']],
 'Y':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD1']],
 'W':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD1']]}


def plot_phi_psi_ramachandran(pdf, PhiDF, PsiDF):
	print('Plotting Phi-Psi')

	global RAMA_PREF_VALUES

	if RAMA_PREF_VALUES is None:
		RAMA_PREF_VALUES = _cache_RAMA_PREF_VALUES()

	for res in PhiDF.index.to_list():
		normals = {"phi":[],"psi":[]}
		outliers = {"phi":[],"psi":[]}
		aa_type = PhiDF.loc[res,'type']
		for mnum in range(1,21,1):
			if RAMA_PREF_VALUES[aa_type][int(PsiDF.loc[res,mnum])+ 180][int(PhiDF.loc[res,mnum]) + 180] < RAMA_PREFERENCES[aa_type]["bounds"][1]:
				outliers["phi"].append(PhiDF.loc[res,mnum])
				outliers["psi"].append(PsiDF.loc[res,mnum])
			else:
				normals["phi"].append(PhiDF.loc[res,mnum])
				normals["psi"].append(PsiDF.loc[res,mnum])

		fig, ax =plt.subplots(figsize=(3,3))
		plt.imshow(RAMA_PREF_VALUES[aa_type], cmap=RAMA_PREFERENCES[aa_type]["cmap"],
				norm=colors.BoundaryNorm(RAMA_PREFERENCES[aa_type]["bounds"], RAMA_PREFERENCES[aa_type]["cmap"].N),
				extent=(-180, 180, 180, -180))
		ax.scatter(normals["phi"], normals["psi"],marker='o',s= 30,facecolors='blue', edgecolors= 'none', linewidth=1.0)
		ax.scatter(outliers["phi"], outliers["psi"],marker='o',s= 30,facecolors='red', edgecolors= 'none', linewidth=1.0)
		ax.set_xlabel('Phi')
		ax.set_ylabel('Psi')
		ax.set_title(res)
		ax.set_xlim([-180,180])
		ax.set_xticks([-180,-120,-60,0,60,120,180])
		ax.set_yticks([-180,-120,-60,0,60,120,180])
		ax.set_ylim([-180,180])
		ax.grid(visible=True, which='major', axis='both',linestyle='--')
		plt.tight_layout()
		pdf.savefig()
		plt.close()

def plot_chi1_chi2_ramachandran(pdf, chi1DF, chi2DF):
	print('Plotting Chi1-Chi2')
	global ROTA_PREF_VALUES

	if ROTA_PREF_VALUES is None:
		ROTA_PREF_VALUES = _cache_ROTA_PREF_VALUES()

	for res in chi1DF.index.to_list():
		normals = {"chi1":[],"chi2":[]}
		outliers = {"chi1":[],"chi2":[]}
		aa_type = res[0]
		for mnum in range(1,2,1):
			if ROTA_PREF_VALUES[aa_type][int(chi2DF.loc[res,mnum])][int(chi1DF.loc[res,mnum])] < ROTA_PREFERENCES[aa_type]["bounds"][1]:
				outliers["chi1"].append(chi1DF.loc[res,mnum])
				outliers["chi2"].append(chi2DF.loc[res,mnum])
			else:
				normals["chi1"].append(chi1DF.loc[res,mnum])
				normals["chi2"].append(chi2DF.loc[res,mnum])

		fig, ax =plt.subplots(figsize=(3,3))
		plt.imshow(ROTA_PREF_VALUES[aa_type], cmap=ROTA_PREFERENCES[aa_type]["cmap"],
				norm=colors.BoundaryNorm(ROTA_PREFERENCES[aa_type]["bounds"], ROTA_PREFERENCES[aa_type]["cmap"].N),
				extent=(0, 360, 360, 0))
		ax.scatter(normals["chi1"], normals["chi2"],marker='o',s= 30,facecolors='blue', edgecolors= 'none', linewidth=1.0)
		ax.scatter(outliers["chi1"], outliers["chi2"],marker='o',s= 30,facecolors='red', edgecolors= 'none', linewidth=1.0)
		ax.set_xlabel('chi1')
		ax.set_ylabel('chi2')
		ax.set_title(res)
		ax.set_xlim([0,360])
		ax.set_xticks([0,60,120,180,240,300,360])
		ax.set_yticks([0,60,120,180,240,300,360])
		ax.set_ylim([0,360])
		ax.grid(visible=True, which='major', axis='both',linestyle='--')
		plt.tight_layout()
		pdf.savefig()
		plt.close()




in_pdb = sys.argv[1]
###----------------------------------------------------------------------------###
##	Extract the starting and ending lines for the 20 NMR models in the file.	##
##	Build a separate dictionary for each set of coordinates with the format:	##
##								A#-atom:[x,y,z] 								##
##	where A# is the residue single letter code and index 						##
###----------------------------------------------------------------------------###
def extract(in_pdb, Sequence, outdir):
	outname = in_pdb.split('/')[-1]
	Starts, Ends = [], []
	for mnum in range(1,21,1):
		start = open(in_pdb).readlines().index('MODEL' + '%9s\n' %str(mnum))
		exec('Coor' + str(mnum) + ' = {}')
		Starts.append(start)
		Ends.append(start-1)
	Ends.append(len(open(in_pdb).readlines()))
	Ends = sorted(Ends)[1:]
	pdb = open(in_pdb).readlines()
	n = 0
	Starts = (sorted(Starts))
	for (start,end) in zip(Starts,Ends):
		n+=1
		# print('Reading coordinates for model %d' %n)
		Coor = eval('Coor' + str(n))
		for x in range(start,end,1):
			line = pdb[x]
			if line[0:4] == "ATOM" or line[0:4] == 'HETA':
				index = '{:}{:}-{:}'.format(AAA_dict[line[17:20].strip()],line[22:26].strip(),line[12:16].strip())
				Coor[index] = [float(line[30:38]),float(line[38:46]),float(line[46:54])]
	PhiDF =  pd.DataFrame(columns=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,'mean','stdv','type'])
	PsiDF =  pd.DataFrame(columns=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,'mean','stdv'])
	chi1DF = pd.DataFrame(columns=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,'mean','stdv'])
	chi2DF = pd.DataFrame(columns=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,'mean','stdv'])
	#chi3DF = pd.DataFrame(columns=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,'mean','stdv'])
	#chi4DF = pd.DataFrame(columns=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,'mean','stdv'])
	#chi5DF = pd.DataFrame(columns=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,'mean','stdv'])

	for mnum in range(1,21,1):
		Coords = eval('Coor' + str(mnum))
		# print('Coor' + str(mnum))
		# print(Coords)
		for i in range(1,len(Sequence)-1,1):
			phi = calcDihedrals(Coords[Sequence[i-1]+ '-C'],Coords[Sequence[i]+ '-N'],Coords[Sequence[i]+ '-CA'],Coords[Sequence[i]+ '-C'])
			PhiDF.loc[Sequence[i],mnum] = np.round(phi,1)
			psi = calcDihedrals(Coords[Sequence[i]+ '-N'],Coords[Sequence[i]+ '-CA'],Coords[Sequence[i]+ '-C'],Coords[Sequence[i+1]+ '-N'])
			PsiDF.loc[Sequence[i],mnum] = np.round(psi,1)
		for res in Sequence: 
			if res[0] in SideDihe.keys():
				for dihe in SideDihe[res[0]]:
					diheDF = eval(dihe[0] + 'DF')
					if res+ '-' + dihe[1] in Coords.keys() and res+ '-' + dihe[2] in Coords.keys() and res+ '-' + dihe[3] in Coords.keys() and res+ '-' + dihe[4] in Coords.keys():
						ang = calcDihedrals(Coords[res+ '-' + dihe[1]],Coords[res+ '-' + dihe[2]],Coords[res+ '-' + dihe[3]],Coords[res+ '-' + dihe[4]])
						if ang < 0: ang = ang + 360.0
						diheDF.loc[res,mnum] = np.round(ang,1)


	PhiDF['mean'] = np.round(PhiDF.mean(axis=1),2)
	PhiDF['stdv'] = np.round(PhiDF.std(axis=1),2)
	PsiDF['mean'] = np.round(PsiDF.mean(axis=1),2)
	PsiDF['stdv'] = np.round(PsiDF.std(axis=1),2)
	chi1DF['mean'] = np.round(chi1DF.mean(axis=1),2)
	chi1DF['stdv'] = np.round(chi1DF.std(axis=1),2)
	chi2DF['mean'] = np.round(chi2DF.mean(axis=1),2)
	chi2DF['stdv'] = np.round(chi2DF.std(axis=1),2)

	for i in range(1,len(Sequence)-1,1):
		if Sequence[i+1][0] == 'P': PhiDF.loc[Sequence[i],'type'] = "PRE-PRO"
		if Sequence[i][0] == 'P':PhiDF.loc[Sequence[i],'type'] = "PRO"
		if Sequence[i][0] == 'G':PhiDF.loc[Sequence[i],'type'] = "GLY"
		else: PhiDF.loc[Sequence[i],'type'] = "General"
	for res in Sequence: 
		if res[0] in SideDihe.keys():
			chi1DF.loc[res,'type'] = res[0]

	PhiDF.to_csv(outdir + outname.replace('.pdb', '_Phi.csv'))
	PsiDF.to_csv(outdir + outname.replace('.pdb', '_Psi.csv'))
	chi1DF.to_csv(outdir + outname.replace('.pdb', '_Chi1.csv'))
	chi2DF.to_csv(outdir + outname.replace('.pdb', '_Chi2.csv'))
	#chi3DF.to_csv(outdir + in_pdb.replace('.pdb', '_chi3.csv'))
	#chi4DF.to_csv(outdir + in_pdb.replace('.pdb', '_chi4.csv'))

	pdf = PdfPages(outdir + outname.replace('.pdb', '_Phi-Psi.pdf'))
	plot_phi_psi_ramachandran(pdf,PhiDF, PsiDF)
	pdf.close()

	pdf = PdfPages(outdir + outname.replace('.pdb', '_Chi1-Chi2.pdf'))
	plot_chi1_chi2_ramachandran(pdf,chi1DF, chi2DF)
	pdf.close()

	# for res in chi1DF.index.to_list():
	# 	fig, ax = plt.subplots(figsize=(3,3))
	# 	ax.scatter(chi1DF.loc[res][:-2],chi2DF.loc[res][:-2],marker='o',s= 30,facecolors='none', edgecolors= 'blue', linewidth=1.0)
	# 	ax.set_xlabel('Chi1')
	# 	ax.set_ylabel('Chi2')
	# 	ax.set_title(res)
	# 	ax.set_xlim([0,360])
	# 	ax.set_xticks([0,60,120,180,240,300,360])
	# 	ax.set_yticks([0,60,120,180,240,300,360])
	# 	ax.set_ylim([0,360])
	# 	ax.grid(visible=True, which='major', axis='both',linestyle='--')
	# 	plt.tight_layout()
	# 	pdf.savefig()
	# 	plt.close()
	# pdf.close()
