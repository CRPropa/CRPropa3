## Using CRPropa from C++
Although is highly recommended to use and [to extend](Extending CRPropa) CRPropa with the default Python binding, there is also a possibility to use it within a C++ code when Python or SWIG are not available or desirable.

To include CRPropa classes and definitions it is sufficient to include ``crpropa/CRPropa.h``, example:

```c++
#include "crpropa/CRPropa.h"

using namespace crpropa;

int main(void){

        ModuleList sim;
        sim.add(new SimplePropagation(1*kpc, 10*Mpc));
        sim.add(new Redshift());
        sim.add(new PhotoPionProduction(CMB));
        sim.add(new PhotoPionProduction(IRB));
        sim.add(new PhotoDisintegration(CMB));
        sim.add(new PhotoDisintegration(IRB));
        sim.add(new NuclearDecay());
        sim.add(new ElectronPairProduction(CMB));
        sim.add(new ElectronPairProduction(IRB));
        sim.add(new MinimumEnergy(1*EeV));

        ref_ptr<Observer> obs = new Observer();
        obs->add(new ObserverPoint());
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
        sim.run(source, 2000, true);

        return 0;
}
```

Compiler such as ``gcc`` should have an access to the header and to CPRropa's shared library (``libcrpropa.so``). If one used paths from [here](Installation), gcc line would look like:
```
g++ example.cpp -o run -I$HOME/.local/include/ -L$HOME/.local/lib/ -lcrpropa
```
however, Makefile should be employed in a general case.
