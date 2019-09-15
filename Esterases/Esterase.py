import sys, os, glob, math
import pandas as pd
#from ResidueTypes import NON_POLAR, POLAR, AROMATIC, POSITIVE_CHARGE, NEGATIVE_CHARGE

NON_POLAR = ["gly", "ala", "val", "leu", "ile", "pro","met"]
AROMATIC = ["phe", "tyr", "trp"]
POLAR = ["ser", "thr", "cys", "asn", "gln"]
POSITIVE_CHARGE = ["lys", "arg", "his"]
NEGATIVE_CHARGE = ["asp", "glu"]

#GENERATE A PARENT CLASS TO COPE WITH PDB OBJECTS

__all__ = ["Esterase", "ActiveSiteDescriptors"]
__author__="Ruben Canadas Rodriguez"
__mail__="ruben.canadas@bsc.es"
__maintainer__="Ruben Canadas Rodriguez"
__version__=1.0


class Error(Exception):

	pass

class ActiveSiteParserError(Error):

	pass


class SiteMapError(Error):

	pass



class Esterase(object):

	"""

	This class have the methods to find the canonical catalytic triad for an esterase given a pdb file.
	It basically tries to find Ser-His-Asp triad using some distance thresholds. Maybe some esterase have differences
	like a glutamic acid instead of an aspartate
	TODO: Implement glutamic option!!
	
	"""


	def __init__(self, pdb, verbose=True):

		self.__pdb = pdb
		self.__protein = self._PDBParser()
		self.__verbose = verbose


	def __len__(self):

		"""
		Returns the length of the file
		"""

		return len(self.__protein)

	def __str__(self):

		return "esterase: {}".format(self.__pdb)

	

	def _PDBParser(self):
		try:
			with open(self.__pdb ,"r") as infile:
				return infile.readlines()
		except IOError:
			exit("pdb file does not exist!")


	def _ResidueDetector(self, residue_name):


		resid = {}
		if residue_name == "SER": atom_name = "HG"
		elif residue_name == "HIS_1": #We have his_1 and his_2 since we are taking into account two atoms for histidine since it interacts with serine and aspartate
			residue_name = "HIS"
			atom_name = "NE2"
		elif residue_name == "HIS_2":
			residue_name = "HIS"
			atom_name = "HD1"
		elif residue_name == "ASP": atom_name = "OD2"
		elif residue_name == "GLU": atom_name = "OE2"
		elif residue_name == "LYS": atom_name = "NZ"
		elif residue_name == "TYR": atom_name = "OH"	
		for line in self.__protein:
			if (line[0:6].strip() == "ATOM") and (line[17:20].strip() == residue_name) and (line[12:16].strip() == atom_name):
				x = float(line[30:38].strip())
				y = float(line[38:46].strip())
				z = float(line[46:54].strip())
				residue_number = int(line[22:26])
				resid["{}_{}".format(residue_name, residue_number)] = [x,y,z]
		return resid


	def _DetectActiveSite(self):

		"""
		This method finds the catalytric triad of a given 
		esterase (in a pdb) according to pre-fixed 
		thresholds
		"""

		def ComputeDistance(atom1, atom2):

			r = [atom2[0] - atom1[0], atom2[1] - atom1[1], atom2[2] - atom1[2]]
			return math.sqrt(r[0] ** 2 + r[1] ** 2 + r[2]**2)

		residues = ["SER", "HIS_1", "HIS_2", "ASP", "GLU", "LYS", "TYR"]
		result = list(map(self._ResidueDetector, residues)) #In python 3 it returns an iterator, thus list(map())
		catalytic_residues = None
		active_site_type = None
		ser_his = []
		his_asp = []
		for key1, value1 in result[0].items():
			for key2, value2 in result[1].items():
				if (ComputeDistance(value1, value2) <= 4.5):
					ser_his.append("{}-{}".format(key1,key2))

		for key3, value3 in result[2].items():
			for key4, value4 in result[3].items():
				if(ComputeDistance(value3, value4) <= 4.5): 
					his_asp.append("{}-{}".format(key3, key4))

		#We take those residues which have the same histidine, meaning that  his is interacting with a serine and an aspartate and the Cat. triad is being formed!
		for elem1 in ser_his:
			serine = elem1.split("-")[0]
			histidine1 = elem1.split("-")[-1]
			for elem2 in his_asp:
				histidine2 = elem2.split("-")[0]
				aspartic = elem2.split("-")[-1]
				if histidine1==histidine2: catalytic_residues = "{}-{}-{}".format(serine, histidine1, aspartic)

		his_glu = []
		if catalytic_residues is None: #If the aspartate is not found is because maybe there is a glutamate instead
			for key3, value3 in result[2].items():
				for key5, value5 in result[4].items():
					if(ComputeDistance(value3, value5) <= 4.5): 
						his_glu.append("{}-{}".format(key3, key5))


		for elem1 in ser_his:
			serine = elem1.split("-")[0]
			histidine1 = elem1.split("-")[-1]
			for elem2 in his_glu:
				histidine2 = elem2.split("-")[0]
				aspartic = elem2.split("-")[-1]
				if histidine1==histidine2: catalytic_residues = "{}-{}-{}".format(serine, histidine1, aspartic)

		ser_lys = []
		lys_tyr = []
		if catalytic_residues is None:
			for key6, value6 in result[0].items():
				for key7, value7 in result[5].items():
					if (ComputeDistance(value6, value7) <= 4.0):
						ser_lys.append("{}-{}".format(key6, key7))

		if catalytic_residues is None:
			for key8, value8 in result[5].items():
				for key9, value9 in result[6].items():
					if (ComputeDistance(value8, value9) <= 4.0):
						lys_tyr.append("{}-{}".format(key8, key9))

		for elem1 in ser_lys:
			ser = elem1.split("-")[0]
			lys1 = elem1.split("-")[-1]
			for elem2 in lys_tyr:
				lys2 = elem2.split("-")[0]
				tyr = elem2.split("-")[-1]
				if lys1 == lys2:
					catalytic_residues = "{}-{}-{}".format(ser, lys1, tyr)

		# Classify the type of active site
		if "ASP" in catalytic_residues:
			active_site_type = 0

		elif "GLU" in catalytic_residues:
			active_site_type = 1

		elif "LYS" in catalytic_residues:
			active_site_type = 2

		if self.__verbose:
			print("catalytic triad of {} is {}: ".format(self.__pdb, catalytic_residues))
		#We need to return the ["SER", num_residue] form to use it in ProteinTopology module
		return catalytic_residues, ["SER", int(catalytic_residues.split("-")[0].split("_")[-1])],active_site_type


##################################################################################################################################################################


class ActiveSiteDescriptors(object):


	def __init__(self, pdb, resid, active_type, radi=10):

		"""
		Constructor of the class
		"""

		self.__resid = resid #list indicating type and num. ex: ["Glu", 88]
		#self.__selection_type = "CA"
		self.__radius = radi
		self.__pdb = pdb
		self.__lines = self.__PDBParser()
		self.__coords = self._ActiveSiteParser()
		self.__neighboring_residues = self._GetNeighboringAtoms()
		self.__active_type = active_type


	def __str__(self):
		
		return "{}".format(self.__pdb)

	def __PDBParser(self):

		try:
			infile = open(self.__pdb, "r")
		except IOError:
			exit("PDB {} could not be opened".format(self.__pdb))

		lines = infile.readlines(); infile.close()
		return lines

	def __CreateResultsPath(self, name=""):

		if not os.path.exists(os.path.join(os.getcwd(), name)):
			os.mkdir(os.path.join(os.getcwd(), name))
		return os.path.join(os.getcwd(), name)

	@property
	def residue(self):
		return self.__resid

	@property
	def radius(self):
		return self.__radius

	@radius.setter
	def radius(self, value):

		self.__radius = value

	def _ActiveSiteParser(self):

		try:
			coords = []
			for line in self.__lines:
				if (line[17:20].strip() == self.__resid[0]) and (int(line[22:26].strip()) == self.__resid[1]): 
				#and (line[12:16].strip() == self.__selection_type):
					x = float(line[30:38].strip())
					y = float(line[38:46].strip())
					z = float(line[46:54].strip())
					coords.append([x,y,z])
			return coords
		except:
			raise ActiveSiteParserError("Active site could not be parsed")


	def _GetNeighboringAtoms(self):

		def ComputeDistance(atom1, atom2):


			r = [abs(atom2[0] - atom1[0]), abs(atom2[1] - atom1[1]), abs(atom2[2] - atom1[2])]
			return math.sqrt(r[0] ** 2 + r[1] ** 2 + r[2]**2)

		path = self.__CreateResultsPath(name="results_topology_active_site")
		neighbors_file = open(os.path.join(path, "neighbors_{}.xyz".format(self.__pdb.split("/")[-1])), "w")
		neighboring_residues = []
		for line in self.__lines:
			if ("ATOM" in line): # and (line[12:16].strip() == "CA"):
				x = float(line[30:38].strip())
				y = float(line[38:46].strip())
				z = float(line[46:54].strip())
				resid_name = line[17:20]
				resid_num = line[22:26]
				atom_coords = [x,y,z]
				for x,y,z in self.__coords:
					if ComputeDistance(atom_coords, [x,y,z]) <= self.__radius:
						residue = [str(resid_name), int(resid_num)]
						if residue not in neighboring_residues:
							neighboring_residues.append(residue)
						neighbors_file.write("{},{},{}\n".format(x,y,z))
		neighbors_file.close()
		return neighboring_residues


	def _ActiveSiteResidueTypes(self):

		"""
		This method computes the percentatge of types of residues
		Types are based on physicochemical properties of the aminoacids
		"""

		resids = [elem[0] for elem in self.__neighboring_residues]
		total_resids = len(resids)
		types = {"polar": 0, "non_polar": 0, "aromatic": 0, "positive_charge": 0, "negative_charge": 0}
		for resid in resids:
			if resid.lower() in POLAR:
				types["polar"] += 1
			elif resid.lower() in NON_POLAR:
				types["non_polar"] += 1
			elif resid.lower() in AROMATIC:
				types["aromatic"] += 1
			elif resid.lower() in POSITIVE_CHARGE:
				types["positive_charge"] += 1
			elif resid.lower() in NEGATIVE_CHARGE:
				types["negative_charge"] += 1
			else:
				raise ValueError("Non-existing residue type!")

		return {key:round(value/float(total_resids),3) for key, value in types.items()}, self.__neighboring_residues, self.__active_type #Dividing the number of residues for the total number of residues in active site

#################################################################################################################################################################

class SiteMapDescriptors(object):


	def __init__(self, pdb, path):

		self.__pdb = pdb
		self.__path = path
		self.__schrodinger_utilities_path = "/opt/schrodinger2018-4/utilities"
		self.__schrodinger_path = "/opt/schrodinger2018-4"
		self.__cat, self.__resid = Esterase(os.path.join(self.__path, self.__pdb))._DetectActiveSite()
		self.__jobname = self.__pdb.split("/")[-1].split(".")[0]

	def ConvertPDBToMAE(self):

		os.system("{}/pdbconvert -ipdb {} -omae {}".format(self.__schrodinger_utilities_path, os.path.join(self.__path, self.__pdb), os.path.join(self.__path, self.__pdb.split(".")[0]+".mae")))

	def SiteMap(self):

		"""
		verbosity 3 in order to obtain the eval file with 
		the desired properties
		"""
		try:
			os.system("{}/sitemap -siteasl \"res.num {}\" -sitebox 5 -maxsites 1 -maxvdw  0.8 -maxdist 12 -enclosure 0.3 -LOCAL -verbosity 3 -j {} -prot {}".format(self.__schrodinger_path, self.__resid[1], self.__jobname, os.path.join(self.__path,self.__pdb.split(".")[0]+".mae")))
		except Exception as e:
			raise SiteMapError("SiteMap could not be used: {}".format(e))
		
	def ParseSiteMapFiles(self):

		file_name = self.__jobname+"_site_1_eval.log"
		infile = open(file_name, "r")
		lines = infile.readlines(); infile.close()
		nums_list = None
		#properties = ["SiteScore", "size", "Dscore", "volume", "exposure", "enclosure", "contact", "hydrophobic", "hydrophilic", "balance", "don/acc"]
		for idx,line in enumerate(lines):
			if "SiteScore size" in line:
				scores = lines[idx]
				nums = lines[idx+1]
				nums_list = [nums.split("   ")[i].strip() for i in range(11)]
		return nums_list


if __name__=="__main__":

	ComputeProps = True
	CreateSiteMap = False
	#PATH = "/home/rubencr/Desktop/esterase_project/GraphCalculations/prova"
	#PATH ="/home/rubencr/Sergi/esterase_project/Protein_structures/Site_Map"
	pdbs = [files for files in os.listdir(PATH) if files.endswith("pdb")]
	properties = ["SiteScore", "size", "Dscore", "volume", "exposure", "enclosure", "contact", "hydrophobic", "hydrophilic", "balance", "don/acc"]
	names = []
	df = pd.DataFrame(columns=properties)
	if CreateSiteMap:
		for pdb in pdbs:
			try:
				site = SiteMapDescriptors(pdb, PATH)
				site.ConvertPDBToMAE()
				site.SiteMap()

			except Exception as e:
				print("error {}".format(e))
				continue
	if ComputeProps:
		for idx,pdb in enumerate(pdbs):
			try:
				site = SiteMapDescriptors(pdb, PATH)
				propert = site.ParseSiteMapFiles()
				print("pro ", properties)
				names.append(pdb)
				for i,prop in enumerate(properties):
					print("prop1 {} prop2 {}".format(prop, properties[i]))
					df.loc[idx, prop] = propert[i]
			except:
				continue

		df["esterase"] = names
		df.set_index("esterase", inplace=True); df.sort_index(inplace=True)
		df.to_csv(os.path.join(PATH, "SiteMapDescriptors.csv"))
		print("datafra ", df)
			#SiteMapDescriptors().ConvertPDBToMAE()
			#SiteMapDescriptors().SiteMap()
	#SiteMapDescriptors().ParseSiteMapFiles()
