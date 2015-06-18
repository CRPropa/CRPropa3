The Doxygen documentation is tracked in the branch gh-pages and is hosted on https://crpropa.github.io/CRPropa3.
The documentation has to be updated and committed manually.
To update use the following procedure:

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