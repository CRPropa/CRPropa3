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
