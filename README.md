Summer Program
==============

Repository for all summer program related programming. If you have any questions,
remember that Google and [Stack Overflow](http://stackoverflow.com/) are your best
friends.

Projects
--------

<ol start="0">
  <li><a href="https://github.com/CCQC/summer-program/blob/master/0/instructions.pdf">Molecule Class</a></li>
  <li><a href="https://github.com/CCQC/summer-program/blob/master/1/instructions.pdf">Vibrational Frequencies</a></li>
  <li><a href="https://github.com/CCQC/summer-program/blob/master/2/instructions.pdf">Finite-difference Hessians</a></li>
  <li><a href="https://github.com/CCQC/summer-program/blob/master/3/instructions.pdf">Closed-shell RHF</a></li>
  <li><a href="https://github.com/CCQC/summer-program/blob/master/4/instructions.pdf">Closed-shell MP2</a></li>
  <li><a href="https://github.com/CCQC/summer-program/blob/master/5/instructions.pdf">UHF</a></li>
  <li><a href="https://github.com/CCQC/summer-program/blob/master/6/instructions.pdf">MP2</a></li>
</ol>


Python
------
Download the [Anaconda package with Python3.4](http://continuum.io/downloads#34)
(you'll have to click on 'I want Python3.4'). To learn how to code Python, use
[Learn Python the Hard Way](http://learnpythonthehardway.org/book/). It will take
time, but this is the best resource I have found, as it skips all the hand-holding
that other resources do. One note, all print lines will need parantheses, a major
change from Python2. Also, starting on excersise 6, see if you can learn how to use
the [string format](https://docs.python.org/3.5/library/string.html#string-formatting)
instead of the '%' to format strings.

Git
---
To learn how to use the basics of Git, follow the Github [Git
tutorial](https://try.github.io/). For a more advanced understanding of how Git
works and for an overview of what version control is, read the [git-scm
tutorial](http://git-scm.com/book/en/v2/Getting-Started-About-Version-Control).

A few helpful hints
* The usual workflow is
    - `git pull`
    - do some work
    - `git add my_files.ext`
    - `git commit -m "My message"`
    - `git push`
* Always remember to use `git pull` before doing a `git push`, otherwise it 
may give you an error message if your local repository is not up to date with
the one on github.
* When commiting code, always provide at least a one-line summary of the
changes you made. If you need to write more, follow the [following
guidelines](http://chris.beams.io/posts/git-commit/)


Organization
------------

 ```
summer-program
│   README.md
│
├───0
│   │   instructions.pdf
│   │
│   ├───username1
│   │   code.py
│   │
│   ├───username2
│   │   code.py
│   ...
│
├───1
...
```

