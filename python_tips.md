# Python Tips and Best Practices

Now that you have learned the basics of Python, you may have realized that there
are many different ways to write code. A big problem with large codebases is
that you can't always tell what each line of someone else's code does, or worse,
you can't remember what each line of your own code does.

## Pass by Assignment

```python
>>> a = 1
>>> b = a
>>> a = 2
>>> print(a, b)
2 1
>>> a = list(range(10))
>>> b = a
>>> a[0] = 5
>> print(b)
[5,1,2,3,4,5,6,7,8,9]
```

What happened here? In the first block, when we changed `a`, `b` stayed the
same. In the second block, we changed `a` and `b` changed as well. Python passes
values by assignment. This means that instead of copying each value in the lists
to a new list, it sets each variable as a reference to the same list. This
reduces the amount of copying necessary, but can lead to some headaches if you
are not careful.

```python
>>> spam = 1
>>> def change_me(var):
	  print(var)
	  var += 1
	  print(var)
>>> change_me(spam)
1
2
>>> spam
1
```

In this case we set spam to a value, and then when we called `change_me()`, we
incremented `var` by one. This did  not change spam, as it is not in the same scope.
Instead, when the function was called, the value of `var` was set to that of `spam,
1. Then var was incremented. However, if we were to have passed a list and
incremented the value of an item in the list, the underlying list would then change.

```python
>>> spam = [0,1,2,3]
>>> def change_me(var):
	  print(var)
	  var[0] = 9
	  print(var)

>>> change_me(spam)
[0, 1, 2, 3]
>>> [9, 1, 2, 3]
>>> spam
[9, 1, 2, 3]
```

In this case, spam is in the global scope. When we pass it to change_me(),
'var' now references the same list as 'spam'. When we change the list that var
points to, it changes the list that 'spam' points to. This is very important to
be aware of as we should be careful when we write functions that mutate the
input variables, as it has the potential to lead to some confusing bugs if we
don't realize that is happening.

Here is another example of Python not acting exactly as you would expect it to.

```python
>>> a = [0]*4
>>> print(a)
[0, 0, 0, 0]
>>> a[0] = 1
>>> print(a)
>>> [1, 0, 0, 0]
>>> b = [[0]*4]*4
>>> print(b)
[[0, 0, 0, 0],
 [0, 0, 0, 0],
 [0, 0, 0, 0],
 [0, 0, 0, 0]]
>>> b[0][0] = 4
>>> print(b)
[[4, 0, 0, 0],
 [4, 0, 0, 0],
 [4, 0, 0, 0],
 [4, 0, 0, 0]]
```

Situations like this are often best served by one of the many libraries in Python,
In this case, matrix class in NumPy might be useful.


## Scoping

```python
>>> i = 1
>>> spam = 1
>>> for i in range(10):
	  spam += i
>>> print(spam, i)
46 9
>>> for i in range(10):
	  eggs = 1
	  for j in range(10):
		  eggs += j
>>> print(eggs, i, j)
46 9 9
```

The previous chunk of code can be very dangerous. While this code will work in
Python, it won't work in C++ as `eggs`, `i`, and `j` would be considered out of
scope. If the code breaks before `i` reaches the end of the range, it can lead
to unexpected results, but can also be useful if you want to determine when the
break happened.

There are two types of scopes in Python

* local - variables defined inside of the current def block
* global - variables defined for the entire file

```python
>>> eggs = 1
>>> def spam():
		print(eggs)
  	spam()
1
>>> def spam():
		print(eggs)
		eggs += 1
UnboundLocalError local variable 'eggs' referenced before assignment
```

This is a classic example of why you need to be careful of scoping. In the
first block, the compiler recognizes 'eggs' as a variable from a higher scope.
When it is asked to print, it knows exactly what to do. However, in the second
block, a needs to be a local variable in order to set its value. The compiler
recognizes this, and prevents us from printing the value of eggs, as it is now
a local variable and its value has not been defined.

Since each file that you make can be imported as a module, all global variables
in a module can be accessed by any file that imports the module.


## Modules
When importing from a module, it is often easiest to import all the functions
from a module as follows.

```python
	from spam import *
```

While this works well for small modules that you wrote, it can be dangerous for
larger modules where you don't know everything that is in it. What if your module
had a function called range that iterated over letters. If you imported it and
overwrote the built-in range function, you would probably end up with an error
later on in your code. It is much safer and organized to instead import only
individual functions from a module or to import the

```python
>>> from spam import my_function
>>> my_function()

>>> import spam
>>> spam.my_function()

>>> import really_long_name_module as mod
>>> mod.my_function()
```

Doing this will help prevent you from cluttering your namespace with
uneccessary functions that you might end up writing over, as well as
conceptually organizing everything. A good example of this is with the math
module.

```python
from math import sqrt
import math
```

The first option limits what you import, and is useful if you only need the
square root function. The second is useful if you are going to be doing a
lot of math, as you then you don't have a whole bunch of functions like
ciel(), fabs(), factorial(), floor(), etc. that you might end up accidently
writing over.

Later in my lecture, we will build our own module and import it into another
Python script.


## Default Arguments

Sometimes we want to write functions with default values for our arguments.

```python
>>> def spam(eggs=1):
		print(eggs)
>>> spam()
1
>>> spam(9)
9
```

Here, if we call spam() it will default to printing 1, but we could also change
the value to whatever we want. We can also have mandatory arguments; they must
be in the beginning of the function declaration.

```python
>>> def spam(a, b, c=3):
		print(a, b, c)
>>> spam(1, 2)
1 2 3
>>> spam(1, 2, 'hey')
1 2 'hey'
```


## Iterators

When we have a `for` loop, we typically iterate over a list of values. We can
make our own iterators if we would like to step through our list on our own.

```python
>>> nums = [1,2,3]
>>> it = iter(nums)
>>> it.next()
1
>>> next(it)
2
>>> it.next()
3
>>> it.next()
StopIteration:
```

How is this useful? We can open a file and iterate over its lines as follows

```python
>>> f = open('geom.xyz')
>>> f is f.__iter__()
>>> f.next()
```

Generators are also the basis for list, set, and dictionary comprehensions.

```python
>>> [ i**2 for i in range(10) ]
>>> { i**2 for i in range(10) }
>>> { i:i**2 for i in range(10) }
```
In Python 2.x, the range function returns a list--potentially a very large one.
In Python 3.x, range returns an iterator like object (technically it is a lazy
evaluated sequence). The advantage of this is that we don't have to store a
large list of numbers that we will discard anyways.

```python
>>> range(10000000)
>>> xrange(10000000)
```

## Documenting Your Code

Documenting your code is always a good idea. While Python is exceptionally
readable, when you come back to your code a month later, you will often wonder
what it does. When documenting a function, it is important to say what the
function is designed to do, what the parameters should be, and what the
function will return. Hopefully your documentation will be good enough that
someone who is looking at your code can know exactly what the function will do,
without needing to look at how you implemented your code.

```python
  # Too little documentation
	def dot(vec1, vec2):
		if len(vec1) != vec2: return None
		return sum([i*j for i,j in zip(a,b)])

  # Excessive documentation
	def dot(vec1, vec2):
		'''
		Calculates the dot product of the two vectors
		
		In mathematics, the dot product, or scalar product (or sometimes inner
		product in the context of Euclidean space), is an algebraic operation
		that takes two equal-length sequences of numbers (usually coordinate
		vectors) and returns a single number. This operation can be defined either
		algebraically or geometrically. Algebraically, it is the sum of the
		products of the corresponding entries of the two sequences of numbers.
		Geometrically, it is the product of the Euclidean magnitudes of the two
		vectors and the cosine of the angle between them. The name "dot product"
		is derived from the centered dot " · " that is often used to designate this
		operation; the alternative name "scalar product" emphasizes the scalar
		(rather than vectorial) nature of the result.
		
		In three-dimensional space, the dot product contrasts with the cross product
		of two vectors, which produces a pseudovector as the result. The dot product
		is directly related to the cosine of the angle between two vectors in
		Euclidean space of any number of dimensions.
    
		:param vec1: the first vector, any iterable containing numbers
		:param vec2: the second vector, any iterable containing numbers
		:return: the dot product
		'''
		
		if len(vec1) != len(vec2):  # if the vectors are not equal
			return None # return nothing
		dot = 0 # Initialize the dot product to 0
		for i in range(len(vec1)):  # for each element in vec1
			dot += vec1[i]*vec2[i]  # add the product of the ith elements to the dot product

		return dot  # return the dot product
```

## Object Oriented Coding

In Python, object-oriented programming is the standard paradigm. What follows
is a very brief overview of some of the basics of object-oriented programming,
and hopefully you will read further into other types of programming (especially
functional programming) and weigh the various merits of the different forms.

```python
class animal(object):
  def __init__( self, name ):
		self.name = name

	def talk(self):
		pass
	
	def whoAmI(self):
		return self.name

class dog(animal):
	def talk(self):
		print( 'woof' )

class cat(animal):
	def talk(self):
		print( 'meow' )

if __name__ == '__main__':
	a = animal('bob')
	c = cat('Evil')
	d = dog('Spot')
	print(a.name)
	print(a.whoAmI())
	a.talk()
	c.talk()
	print(c.whoAmI())
	d.talk()
```

While we are talking about object oriented programming in Python, here are a
few of the hidden things you can do in Python that will make your lives much
easier.

```python
class page(object):
	def __init__(self, words):
		self.words = words
	
	def __len__(self):
		return len(self.words)
	
	def print_page(self):
		print(self.words)

class book(object):
	
	def __init__(self, pages):
		self.pages = []
		for p in pages:
			self.pages.append( page(p) )
  
  def __len__(self):
			return len(self.pages)
	
	def __getitem__(self, pagenum):
		return self.pages[pagenum]
	
	def __setitem__(self, pagenum, value):
		self.pages[pagenum] = value
	
	def print_book(self, begin=0, end=None):
		for book_page in self.pages[begin:end]:
			book_page.print_page()

if __name__ == "__main__":
	p = page('Hello')
	print('Printing page:')
	p.print_page()

	b = book(['I am Sam.', 'Sam I am.', 'Do you like green eggs and ham?', 'The End' ])
	print('\nPrinting book:')
	b.print_book()
	print()
	b[2].print_page()
```

## Learning Numpy
```python
>>> import numpy as np
>>> a = np.array([1,2,3])
>>> b = np.array([3,4,5])
>>> a*b
array([3,8,15])
>>> a.dot(b)
26
>>> np.cross(a,b)
array([-2, 4, -2])
```

Here I created two vectors and was able to perform basic linear algebra
functions on them.

```python
>>> import numpy as np
>>> a = np.matrix( [[1, 2],
                    [-1, 4]] )
>>> a.shape
(2, 2)
>>> print(a.T)
matrix(	[[1, -4],
         [ 2, 4]])
>>> a.I
matrix(	[[0.66666667, -0.3333333],
         [0.16666667, 0.16666667]])
```

The NumPy module is very powerful. Even more advanced features are available in
SciPy and I would highly recommend looking into it. Here are a few

```python	
>>> import numpy as np
>>> from scipy import linalg
>>> a = np.array( [[1, 2],
                   [-1, 4]] )
>>> a.shape
(2, 2)
>>> a.T
array(	[[1, -4],
         [ 2, 4]])
>>> linalg.inv(a)
array(	[[0.66666667, -0.3333333],
         [0.16666667, 0.16666667]])
>>> linalg.eig(a)
(array([ 2.,  3.]),
 array([[-0.89442719, -0.70710678],
        [-0.4472136 , -0.70710678]]))
```

## Code Formatting

In general, all code should be formatted as follows.

* Tab = 4 spaces (in your .vimrc type `set expandtab`)
* Wrap all line at 79 characters
* Place a space on each side of the following: `=`, `+`, `*`, `/`
* Use a space after every comma
* Two newlines before and after every function (unless in a class, then use only one)

If you get bored I would highly recommend reading through [PEP
8](https://www.python.org/dev/peps/pep-0008/), the Python style guideline.
If nothing else it will motivate you to get back to work (or have a nice nap).
