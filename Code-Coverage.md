## Code Coverage

* compile crpropa with coverage spport
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
