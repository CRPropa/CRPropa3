#### Seed the random number generator 
The random number seed can be set with
```python
Random_seedThreads(seed)  # seed from 0 - 2^32-1
```


#### Source positions from a matter density grid

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

#### Version tracking of a python steering card

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