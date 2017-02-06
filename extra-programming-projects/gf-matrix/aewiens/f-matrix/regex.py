import re
from masses import charge

def capture(string):
  return r'({:s})'.format(string)

def zero_or_more(string):
  return r'(?:{:s})*'.format(string)

## Andreas
## def zero_or_one(string):

def one_or_more(string):
  return r'(?:{:s})+'.format(string)

def two_or_more(string):
  return r'(?:{:s}){{2,}}'.format(string)

def isAtom(string):

    atomLabels = [c for c in charge]
    if string.upper() in atomLabels: return True
    else: return False


# fundamental components
character        = r'[a-zA-Z]'
unsigned_integer = r'\d+'
signed_double    = r'-?\d+\.\d+'
space            = r'\s'
endline          = r'\n'
atom             = r'[A-Z][a-z]?'
Label            = r'\w' + zero_or_more( r'\w')
coordLabel       = capture(Label) + r'\s+\=\s+\-?\d+\.\d+\n'
coordinate       = Label + r'\s+\=\s+(\-?\d+\.\d+)\n'


def Zmatrix():
	
	zmatLine1  = atom + zero_or_more(space) + r'\n'
	zmatLine2  = atom + r'\s+\d+\s+\w\w?' + zero_or_more(space) + r'\n'
	zmatLine3  = atom + r'\s+\d+\s+\w\w?\s+\d+\s+\w\w?' + zero_or_more(space) + r'\n'
	zmatLine4  = atom + r'\s+\d+\s+\w\w?\s+\d+\s+\w\w?\s+\d+\s+\w\w?' + zero_or_more(space) + r'\n'

	coordBlock = one_or_more(Label + r'\s+\=\s+\-?\d+\.\d+' + zero_or_more(space) + r'\n')
	zmatBlock  = zmatLine1 + zmatLine2 + zero_or_more(zmatLine3) + zero_or_more(zmatLine4) +  coordBlock

	return zmatBlock
 
if __name__ == '__main__':
	f = open("template.dat","r")
	zmat = f.read()
