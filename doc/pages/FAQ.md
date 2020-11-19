## FAQ

### I have an installation issue...

In case of an installation issue make sure that you did not miss something in
the [Installation](Installation.md)
section of this documentation. Then take a look at the list of already reported [issues
tagged wth
"Installation"](https://github.com/CRPropa/CRPropa3/issues?utf8=%E2%9C%93&q=is%3Aissue+label%3Ainstallation+).
If your issue is not listed, open a new ticket.

Please be aware that we cannot support general issues regarding an operating
system setup, a third party package installation, an unsupported third party
package (libraries, compilers) version or similar, therefore do not report them
to us. In case of doubt whether something is an CRPropa issue or not, report
the problem to the issue tracker.

### How to specify the seed for the random number generator?
The random number seed of the global random number singleton can be set with
```python
Random_seedThreads(seed)  # seed from 0 - 2^32-1
```
For example you get
```python
In [17]: crpropa.Random_seedThreads(42)

In [18]: crpropa.Random_instance().rand()
Out[18]: 0.37454011439684315

In [19]: crpropa.Random_seedThreads(42)

In [20]: crpropa.Random_instance().rand()
Out[20]: 0.37454011439684315
```


### How to define source positions from a matter density grid?

```python
grid = ScalarGrid(Vector3d(0, 0, 0), 256, 256, 256, 100 * kpc)
loadGridFromTxt(grid, 'some_density_grid.txt')
#loadGrid(grid, 'some_density_grid.raw')
density = crpropa.SourceDensityGrid(grid)
```

For continuous source positions
```python
source.add(density)
```

For n discrete source positions
```python
positions = SourceMultiplePositions()
p = ParticleState()
for i in range(n):
    density.prepare(p)
    pos.add(p.getPosition())
source.add(positions)
```

### How to use version tracking of a python steering card?

Good practice to secure the reproducibility of results obtained with CRPropa is to track which version of CRPropa is used for particular steering card (python code). CRPropa provides the following helper function:
``crpropa.declare_version("3.1-135-g9ec850f")``
When a new crpropa version is installed a warning message will be displayed to remind the user of potential differences:
```
2017-12-07 09:32:00 [WARNING] Version mismatch! To clear this warning
review the python code for potential incompatibilities and update
its version declaration or install the declared version of CRPropa.
- CRPropa version: 3.1-136-g75fcd3b
- Python code version: 3.1-135-g9ec850f
Use git diff to inspect the differences:
  git diff 3.1-135-g9ec850f 3.1-136-g75fcd3b
```

A more casual check of versions is also allowed ``crpropa.declare_version("3.1")`` where the warning will be displayed only when the tag is different, e.g., 3.1 != 3.2.

### I did not find an answer to my question...

Before opening [a new ticket](https://github.com/CRPropa/CRPropa3/issues/new), one should first try to:
1. Look at this documentation.
3. Check [previously asked questions](https://github.com/CRPropa/CRPropa3/issues?utf8=%E2%9C%93&q=label%3Ausage-question+).
4. Use search in [the issue tracker](https://github.com/CRPropa/CRPropa3/issues?utf8=%E2%9C%93&q=is%3Aissue).
