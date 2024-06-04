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
The final report is in ```~/build/coverageReport/index.html```

### Coverage report of the current master branch
The coverage report of the current master branch can be accesed `by clicking here <coverageReport/index.html>`_ now.
