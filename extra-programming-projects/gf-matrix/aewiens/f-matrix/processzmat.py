import re,regex
from molecule import Molecule

class InputFile(object):

	def __init__(self,header,footer,body_template):
		
		self.Template = header + body_template + footer

	def makeInput(self,coordinates):

		return self.Template.format(*coordinates)


class ProcessZmat(object):

	def __init__(self, inputString):
		
		inputString  = inputString.replace('{', '{{').replace('}', '}}')
		
		match = None
		for match in re.finditer(regex.Zmatrix(), inputString, re.MULTILINE):
			pass
		if not match:
			raise Exception('Cannot find start of internal coordinate block')

		start, end       = match.start(), match.end()
		zmatBlock        = inputString[start:end]
		header,footer    = inputString[:start],inputString[end:]
		bodyTemplate     = re.sub('\s+(-?\d+\.?\d+)', ' {:> 17.12f}', zmatBlock)
		self.template    = header + bodyTemplate + footer
		self.inputObject = InputFile(header,footer,bodyTemplate)


if __name__ == '__main__':
	f    = open("template.dat","r").read()
	#mol  = Molecule(f)
	zmat = ProcessZmat(f)
	print( zmat.inputObject.makeInput(mol.coords) )
