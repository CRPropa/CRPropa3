# Contributing to CRPropa
First of all thank you for considering contributing to CRPropa. We are happy to
include contributions into the CRPropa repository that are beneficial to the
community. However, to keep the code working and maintainable, we ask any
contributors to keep the coding style of their contributions consistent and
also add documentation and tests along with their code. Consider that your code
is more often read than written.

If you start working on a future contributions, please create a new issue to
notify other users and developers on your planned work early on. Thus others
can collaborate with you early on and double work is avoided.

Ideally new contributors should:

  * fork the CRPropa repository
  * create a feature branch within their repository to work on their project
  * create a pull request when ready

The code will be reviewed before including in CRPropa. To ease the review
process, please adhere to the following coding conventions / requirements.
There should always be a good reason for an exception.

If the code cannot be added to the main repository, please consider releasing
it as a plugin maintained in a separate repository.

## Coding Conventions / Requirements
Contributions to CRPropa should follow two design principles:

  * **KISS**: Keep it simple, not stupid. Any solution / contribution should be
	  as simple as possible, but not simpler. Software (and physics) is complex
		enough without adding unnecessary complexity.
  * **YAGNI**:  You ain't gonna need it. This is basically a consequence of the
    KISS principle. Often designing a new feature one immediately has ideas about
    additional potential uses that can be very easily added in the future, if the
    feature is only made a little bit more general (complex). Adding this
    complexity now should be avoided in favour of later refactoring of the code.
    If the feature is not needed or immediately used, there should be no
    preparatory code added as it requires maintenance. 

In addition to those general ideas, we use the following conventions:

### Equations and Units
  * SI units are used throughout the code in all interfaces. Definitions in
    units.h should be used to make this explicit. Write `foo(1*meter)` instead of
    `foo(1)`
  * The origin of all equations should be referenced inline and properly
    commented

### Code Behavior
  * Code should not output on the terminal, but use the logging facilities
  * Code should not use preprocessor conditionals, e.g. to enable a debugging
    mode.
  * All code should be written to allow parallelization, in particular all
    modules need to be stateless
  * All methods should adhere `const` correctness
  * Dependencies should be trivial to install, i.e. shipped with virtually all
    Linux distributions or be shipped with CRPropa in a specific version.  Users
    should not need to compile dependencies.
  * Local variables should be initialized in same line as declaration and
    preferably be declared in scope of `if`, `while`, `for` statements unless
    this has severe performance implications due to calling constructors.

### Documentation + Tests
  * All definitions of classes etc. should be documented using doxygen
  * All new features should have an usage example in the wiki
  * All code and new functionality should be covered by unit tests, coverage
    should be [checked](https://github.com/CRPropa/CRPropa3/wiki/Code-Coverage)

### C++ Code
  * CamelCase is used for naming, i.e. `FooBar::doSomething` instead of e.g.
    `foo_bar::do_something`.
  * Class names should be capitalized, methods and variables not
  * Names should be rather verbose, e.g. `c_light` vs. `c`, `setDescription` vs.
    `setDesc`, etc.
  * Tabs are used for indentation, spaces for alignment.
  * No spaces should be used in parentheses, but between statement and
    parentheses `if (condition) { ... }` instead of `if( condition){ ... }`
	* Spaces hwould be used around operators, i.e. `x = a * b + 2` instead of
		`x=a*b+2`
  * Brace is on the same line as the statement (K&R style), i.e.
    `while(x == y){`

    `}`
  * Includes should be relative to the base path of CRPropa and not use UNIX
    relative paths such as ., ..
  * Includes should be in order of
    other CRPropa headers, other libraries, C++ standard libs, separated by blank
    lines
  * If possible headers should be included in the implementation (cpp file) and not the definition (.h) file.

### Python Code
  * All code needs to run on python 3 (and 2.7)
  * Code should follow [PEP8](https://www.python.org/dev/peps/pep-0008/)

Additional suggestions for good styled C++-code can be found e.g. in [Google's
C++ style guide](https://google.github.io/styleguide/cppguide.html).
