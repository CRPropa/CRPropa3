## Using CRPropa from C++
Although is highly recommended to use and [to extend](Extending CRPropa) CRPropa with the default Python binding, there is also a possibility to use it within a C++ code when Python or SWIG are not available or desirable.

To include CRPropa classes and definitions it is sufficient to include ``crpropa/CRPropa.h``, example:

```c++
#include "crpropa/CRPropa.h"

using namespace crpropa;

int main(void) {
    ModuleList sim;
    sim.add(new SimplePropagation(1*kpc, 10*Mpc));
    sim.add(new Redshift());
    sim.add(new PhotoPionProduction(new CMB()));
    sim.add(new PhotoPionProduction(new IRB_Kneiske04()));
    sim.add(new PhotoDisintegration(new CMB()));
    sim.add(new PhotoDisintegration(new IRB_Kneiske04()));
    sim.add(new NuclearDecay());
    sim.add(new ElectronPairProduction(new CMB()));
    sim.add(new ElectronPairProduction(new IRB_Kneiske04()));
    sim.add(new MinimumEnergy(1*EeV));

    ref_ptr<Observer> obs = new Observer();
    obs->add(new Observer1D());
    obs->onDetection(new TextOutput("events.txt", Output::Event1D));
    obs->onDetection(new TextOutput());
    sim.add(obs);

    ref_ptr<Source> source = new Source();
    source->add(new SourceUniform1D(1*Mpc, 1000*Mpc));
    source->add(new SourceRedshift1D());

    ref_ptr<SourceComposition> composition =
        new SourceComposition(1*EeV, 100*EeV, -1);
    composition->add(1,  1,  1);
    composition->add(4,  2,  1);
    composition->add(14, 7,  1);
    composition->add(56, 26, 1);
    source->add(composition);

    sim.setShowProgress(true);
    sim.run(source.get(), 2000, true);
}
```

Compiler such as ``gcc`` should have an access to the header and to CRPropa's shared library (``libcrpropa.so``). If one uses the conda package installation as described [here](Installation#Over-Conda-Package-(recommended)), gcc line would look like:
```sh
g++ example.cpp -o run -I$CONDA_PREFIX/include/ -L$CONDA_PREFIX/lib/ -lcrpropa
```
However, a Makefile should be employed in a general case.
