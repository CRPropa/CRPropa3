The Doxygen documentation is tracked in the branch gh-pages and is hosted on https://crpropa.github.io/CRPropa3.
It is updated automatically using travis. 

The documentation can be updated and committed manually using the following procedure. Note that the cahanges are overwritten by the next commit by travis!

1) go to CRPropa3/doc
```
cd doc
```

2) clone the branch gh-pages to html
```
git clone https://github.com/CRPropa/CRPropa3.git --branch gh-pages --single-branch html
```

3) update the documentation
```
doxygen Doxyfile
```

4) commit the changes
```
cd html
git add -u
git commit
git push
```