## GNU Debugger (GDB)

* compile crpropa with debugging symbols:
```
cmake -DCMAKE_BUILD_TYPE:STRING=Debug ..
```

* run a crpropa script within gdb:
```
gdb --args python python_script.py
```
(typing `run` inside gdb will run the script).