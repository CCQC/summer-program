Summer Program
==============

Repository for all summer program related programming. If you have any questions,
remember that Google and [Stack Overflow](http://stackoverflow.com/) are your best
friends.


Projects
--------

<ol start="0">
  <li><a href="https://github.com/CCQC/summer-program/blob/master/0/instructions.pdf">Geometry</a>
     (<a href="https://github.com/CCQC/summer-program/blob/master/0/advanced.pdf">Advanced</a>)</li>
  <li><a href="https://github.com/CCQC/summer-program/blob/master/1/instructions.pdf">Frequencies</a>
     (<a href="https://github.com/CCQC/summer-program/blob/master/1/theory.pdf">Theory</a>)</li>
  <li><a href="https://github.com/CCQC/summer-program/blob/master/2/instructions.pdf">Hessian</a>
     (<a href="https://github.com/CCQC/summer-program/blob/master/2/theory.pdf">Theory</a>)</li>
  <li><a href="https://github.com/CCQC/summer-program/blob/master/3/instructions.pdf">HF</a>
     (<a href="https://github.com/CCQC/summer-program/blob/master/3/theory-part-1.pdf">Theory 1</a>,
      <a href="https://github.com/CCQC/summer-program/blob/master/3/theory-part-2.pdf">Theory 2</a>)</li>
  <li><a href="https://github.com/CCQC/summer-program/blob/master/4/instructions.pdf">MP2</a>
     (<a href="https://github.com/CCQC/summer-program/blob/master/4/exercises.pdf">Exercises</a>)</li>
  <li><a href="https://github.com/CCQC/summer-program/blob/master/5/instructions.pdf">Spin-Orbital UHF</a>
     (<a href="https://github.com/CCQC/summer-program/blob/master/5/exercises.pdf">Exercises</a>)</li>
  <li><a href="https://github.com/CCQC/summer-program/blob/master/6/instructions.pdf">Spin-Orbital MP2</a>
     (<a href="https://github.com/CCQC/summer-program/blob/master/6/exercises.pdf">Exercises</a>)</li>
  <li><a href="https://github.com/CCQC/summer-program/blob/master/7/instructions.pdf">CIS</a>
     (<a href="https://github.com/CCQC/summer-program/blob/master/7/exercises.pdf">Exercises</a>,
      <a href="https://github.com/CCQC/summer-program/blob/master/7/theory.pdf">Theory</a>)</li>
  <li><a href="https://github.com/CCQC/summer-program/blob/master/8/instructions.pdf">CCD</a>
     (<a href="https://github.com/CCQC/summer-program/blob/master/8/theory.pdf">Theory</a>)</li>
</ol>


Extra Programming Projects
--------------------------
 - [Internal Coordinate Vibrational Analysis -- The GF Matrix Method](https://github.com/CCQC/summer-program/tree/master/extra-programming-projects/gf-matrix/)
 - [CEPA0 Energy](https://github.com/CCQC/summer-program/tree/master/extra-programming-projects/cepa0/)
 - [Numerical Solution to the Linear 1D Schrödinger Equation -- The Numerov Method](https://github.com/CCQC/summer-program/tree/master/extra-programming-projects/numerov/)
 - [Spin-Integrated UHF](https://github.com/CCQC/summer-program/tree/master/extra-programming-projects/spin-integrated-uhf/)
 - [UHF Natural Orbitals](https://github.com/CCQC/summer-program/tree/master/extra-programming-projects/uhf-natural-orbitals/)
   ([Theory](https://github.com/CCQC/summer-program/blob/master/extra-programming-projects/uhf-natural-orbitals/theory.pdf))


Extra Applications Projects
---------------------------
 - coming soon...


Algebraic methods
-----------------
 - [Tensors](https://github.com/CCQC/summer-program/blob/master/algebraic-methods/tensors.pdf)
 - [Motivating second quantization](https://github.com/CCQC/summer-program/blob/master/algebraic-methods/motivating-second-quantization.pdf)
 - [Normal ordering and Wick's theorem](https://github.com/CCQC/summer-program/blob/master/algebraic-methods/normal-ordering.pdf)


Python
------
Download the [Anaconda package with Python3.4](http://continuum.io/downloads#34)
(you'll have to click on 'I want Python3.4'). To learn how to code Python, use
[Learn Python the Hard Way](http://learnpythonthehardway.org/book/). It will
take time, but this is the best resource I have found, as it skips all the
hand-holding that other resources do. One note, all print lines will need
parantheses, a major change from Python2. Also, starting on excersise 6, see if
you can learn how to use the [string
format](https://docs.python.org/3.5/library/string.html#string-formatting)
instead of the '%' to format strings. For a general overview of best practices
in Python, read [python_tips.md](python_tips.md).


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
* If you forgot to `git pull` before you committed your work and you did not
edit the same files as anyone else (e.g. you only worked in your directory),
run `git rebase` and it will reorganize commits to remove the unecessary merge.
* When commiting code, always provide at least a one-line summary of the
changes you made. If you need to write more, follow the [following
guidelines](http://chris.beams.io/posts/git-commit/)


Organization
------------

 ```
summer-program/
|   .gitignore
|   python_tips.md
│   README.md
│
├───extra-files/
│
├───0/
│   │   instructions.pdf
│   │
│   ├───username1/
│   │   code.py
│   │
│   ├───username2/
│   │   code.py
│   ...
│
├───1/
...
```

