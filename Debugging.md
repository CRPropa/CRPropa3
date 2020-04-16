## GNU Debugger (GDB)

Compile CRPropa with debugging symbols:
```
cmake -DCMAKE_BUILD_TYPE:STRING=Debug ..
```

Run a crpropa script within gdb:
```
gdb --args python python_script.py
```
(typing `run` inside gdb will run the script).

## Profiling

In programming and especially in numerical simulations which are almost by definition computing intensive, it is critical to discover bottlenecks in a code design or, in other words, to find out where does the program spend most of the time. There are many ways to do this analysis but a common one is to do the code profile.

There are different types of profilers available, but here we will focus on the [Valgrind](http://valgrind.org/) framework and a tool called callgrind. Callgrind can generate call-graphs from which one can easily see where are the slowest parts of a program. Valgrind provides an instrumenting profiling which is accurate, support multi-threaded programs. Unfortunately, it is really slow (5-10 times slower than the normal program execution).

### Simple example

Here is a simple example with a few arbitrary addition just to prevent compiler optimisations. The primary goal is to inspect four functions (`main`, `sleepy1`, `sleepy2` and `sleepy3`) and to see how much time does the program spend at each of them. 

```cpp
/* main.cpp */
#include <iostream>

using namespace std;

void sleepy1( int );
void sleepy2( int );
void sleepy3( int );

void sleepy1( int i ){
    int sum = 0;
    for( int j = 0; j<100; ++j ){
      sum += j+i;
      sleepy3( sum );
    }
}

void sleepy2( int counter ){
    int sum = 0;
    for( int j = 0; j<100000; ++j ){
      sum += counter;
    }
}

void sleepy3( int counter ){
    int sum = 0;
    for( int j = 0; j<1000; ++j ){
      sum += counter;
    }
}

int main( void ){
 
  for( int i = 0; i < 20; ++i ){
    sleepy1( i );
    sleepy2( i );
  }
  
  return 0;
}
```

Compile it with debugging symbols (the "-g" flag) which introduce "markers" in the executable, so that human-readable function names are displayed in the output:
```
g++ -g main.cpp -o main
```
and then run it with Valgrind:
```
valgrind --tool=callgrind ./main
```
The result is saved in callgrind.out.* which can be easily read by, for example, KDE's [kcachegrind](http://kcachegrind.sourceforge.net):
```
kcachegrind callgrind.out.*
```

    Calls list

    Graphical representation

From the resultsm one can observe how many times some function is called, how much time in the relative scale does it cost and so on. There are two types of costs involved: '''self''' and '''incl.''' To understand the difference, we cite KCachegrind's manual:

> What is the difference between Incl. and Self?
>
> These are cost attributes for functions regarding some event type. As functions can call each other, it makes sense to distinguish the cost of the function itself (“Self Cost”) and the cost including all called functions (“Inclusive Cost”). “Self” is sometimes also referred to as “Exclusive” costs.

### Profile CRPropa

To compile CRPropa with debug symbols one should use the Debug directive in CMake:
```
cmake -DCMAKE_BUILD_TYPE:STRING=Debug .. 
```
Valgrind is run as follows:
```
valgrind --tool=callgrind python steering_card.py
```
### References
* [KCachegrind Manual](https://docs.kde.org/stable4/en/kdesdk/kcachegrind/kcachegrind.pdf)
* [Callgrind manual](http://valgrind.org/docs/manual/cl-manual.html)
* [GCC debugging options](https://gcc.gnu.org/onlinedocs/gcc-3.4.5/gcc/Debugging-Options.html)
* [Python profiling](http://www.blog.pythonlibrary.org/2014/03/20/python-102-how-to-profile-your-code/)
