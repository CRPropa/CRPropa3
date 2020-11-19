## Code Coverage
`lcov` is needed for this feature.

* compile crpropa with coverage support
```
cmake -DENABLE_COVERAGE=TRUE ..
```

*  Execute tests
```
make test
```

*  Create coverage report
```
make coverage
```
The final report is in ```coverageReport/index.html```
