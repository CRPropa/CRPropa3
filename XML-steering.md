The XML steering uses 'steering cards', files with a number of commands in XML syntax, to setup and control the simulation. The XML steering is intended to provide backwards compatibility to CRPropa2. Due to the fixed setup the XML steering cannot provide access to all functionalities of CRPropa3.

To run a simulation execute
```sh
 crpropa-xmlrun mysteeringcard.xml
```

For a complete list of xml settings, refer to the CRPropa 2 manual.

#### Changes to CRPropa2 cards
Since CRPropa 2 does not fully enforce the XML-1.0 standard, a few modifications to existing XML cards can be necessary. 
 * XML-cards can have only one root node. The root node is
   ```
   <CRPropa>
   ....
   </CRPropa>
   ```
 * All nodes including the optional header node need to be closed, e.g.
   ```<Option1> ... </Option1>``` or ```<Option2/>```
 * Values need to be embraced by quotes, e.g.
   ```<TrajNumber value = "1000"/>``` or ```<TrajNumber value = '1000'/>```

#### Example
The following XML syntax is an example of a 3D simulation without cosmology or galactic magnetic field deflections.

Each XML card should be opened with the root node
```xml
<CRPropa>
```
The order of the following nodes does not matter, as long as they are enclosed by the root node.

The number of particles to inject at the source, the lowest energy (in EeV) down to which the UHECRs will be followed, and the maximum time (in Mpc) that the particles will be propagated, are set with
```xml
 <TrajNumber    value = "1000" />
 <MinEnergy_EeV value = "1"    />
 <MaxTime_Mpc   value = "2000" />
```

The output type, file type and file name are specified by
```xml
 <Output type = "Events">
    <File type = "ROOT" option = "force"> outputfilename.root </File>
 </Output>
```
Other possibilities instead of "Events" for the output type are "None" and "Full Trajectories". Another possibility for the file type is ROOT. So instead of the previous block another example would be
```xml
 <Output type = "Full Trajectories">
    <File type = "ASCII" option = "force"> outputfilename.txt </File>
 </Output>
```
in which case, instead of only the final state of the particle at an observer, the particle state after each time-step will be written to the output file.

In a 3D simulation the environment in which the simulation is done should be specified like this (for a (20 Mpc)^3 simulation box):
```xml
 <Environment type = "LSS">
    <Xmin_Mpc value = "-10" />
    <Xmax_Mpc value = "10"  />
    <Ymin_Mpc value = "-10" />
    <Ymax_Mpc value = "10"  />
    <Zmin_Mpc value = "-10" />
    <Zmax_Mpc value = "10"  />
 </Environment>
```

For a 1D simulation the environment type should be set to "One Dimension" like this:
```xml
 <Environment type = "One Dimension" />
```

When, in 3D, magnetic fields are included, the details of the environment are specified in the magnetic field box so that, for the environment, this line is enough:
```xml
  <Environment type = "LSS" />
```

An example for a Kolmogoroff-type magnetic field would look like this:
```xml
 <MagneticField type = "Kolmogoroff">
    <Nx value = "128" />
    <Ny value = "128" />
    <Nz value = "128" />
    <Step_Mpc value = "0.5" />
    <Origin>
       <X_Mpc value = "0" />
       <Y_Mpc value = "0" />
       <Z_Mpc value = "0" />
    </Origin>
    <SpectralIndex value = "0.66" />
    <RMS_muG value = "0.001" />
    <Kmin value = "0.025" />
    <Kmax value = "0.5" />
 </MagneticField>
```

However, in the case without magnetic fields the environment details are specified in the "Environment" block and the magnetic field should be set to "Null" like this:
```xml
 <MagneticField type = "Null" />
```

Other available options for the magnetic field type are: "1D", "Uniform" and "LSS-Grid". For detailed examples of these options see the XML files in the "test/" folder after installing CRPropa 3.

Next the interactions used in the simulation should be specified. This should be done in the "Interactions" block, where available types are "None" (switching off all interactions) and "Sophia", where by default all interactions (photodisintegration, pair production, pion production and decay) are switched on. Within this block the maximum step size for one interaction should be specified. So a simple example of this block would be:
```xml
 <Interactions type = "Sophia">
    <MaxStep_Mpc value = "1" />
 </Interactions>
```

Within this block you can switch off photodisintegration, pair production, pion production and decay with the flags "NoPhotodisintegration", "NoPairProd", "NoPionProd" and "NoDecay" respectively. Furthermore it is possible to switch off interactions with the infrared background by adding the flag "NoIRO". In 1D it is possible to switch off redshift evolution by adding the flag "NoRedshift" to this block, in 3D redshift evolution is available as an a posteriori correction or by doing a 4D simulation, however by default redshift evolution is switched off. If secondary photons and/or neutrinos should also be propagated, the flags "SecondaryPhotons" and/or SecondaryNeutrinos" should be added to the "Interactions" block. An example without photodisintegration but including secondary neutrinos would be:
```xml
 <Interactions type = "Sophia">
    <MaxStep_Mpc value = "1" />
    <NoPhotodisintegration   />
    <SecondaryNeutrinos      />
 </Interactions>
```

Furthermore the integrator should be specified by adding the block "Integrator". At the moment the only possible type is "Cash-Karp RK". This integrator requires values for "Epsilon" and the minimal step size 
"MinStep_Mpc" like this:
```xml
 <Integrator type = "Cash-Karp RK">
    <Epsilon value = "1e-5"    />
    <MinStep_Mpc value = "1e-4" />
 </Integrator>
```

Next the source details should be specified in the block "Sources". Within this block the number of sources, their position, the spectrum at the sources and the species of nuclei injected at the sources should be specified. This is an example for one discrete source injecting iron, magnesium and proton (twice as many iron nuclei as the other species) with a power law spectrum of dN/dE∝E^-2.2 up to an energy of 300 EeV for each species:
```xml
 <Sources type = "Discrete">
    <Number value = "1" />
    <PointSource>
       <CoordX_Mpc value = "3" />
       <CoordY_Mpc value = "3" />
       <CoordZ_Mpc value = "3" />
    </PointSource>
    <Spectrum type = "Power Law">
       <Ecut_EeV value = "300" />
       <Alpha    value = "2.2" />
    </Spectrum>
    <Particles type = "Nuclei">
       <Number_Of_Species value = "3" />
       <Species MassNumber = "56" ChargeNumber = "26" Abundance = "20" />
       <Species MassNumber = "24" ChargeNumber = "12" Abundance = "10" />
       <Species MassNumber = "1"  ChargeNumber = "1"  Abundance = "10" />
    </Particles>
 </Sources>
```

Another possibility for the "Sources" type is "Continuous", in which case a "Density" type should be specified as "Uniform" (for an example see the XML files in the "test/" folder) or "Grid", in which case the density should be specified in an ASCII file. Another possibility for the "Spectrum" type is "Monochromatic", in which case all particles will be injected with the same energy or rigidity specified with the flags "Energy_EeV value = " and "Rigidity_EeV value = " respectively. In the case of a power law spectrum, instead of the "Ecut_EeV value = " flag, the flag "Rigidity_EeV value = " can be specified, in which case the maximum energy at injection is the supplied rigidity value times the charge of the injected nucleus. 

And finally, for a 3D event output, the observer needs to be specified. In this example the observer is a sphere with radius of 0.5 Mpc:
```xml
 <Observers type = "Spheres around Observers">
    <Number     value = "1"   />
    <Radius_Mpc value = "0.5" />
    <SphereObserver>
       <CoordX_Mpc value = "0" />
       <CoordY_Mpc value = "0" />
       <CoordZ_Mpc value = "0" />
    </SphereObserver>
 </Observers>
``` 

Another option for the "Observers" type is a "Spheres around Source", in which case only a single source should be specified and the particles are detected when they reach a certain distance from that source (for an example see the XML files in the "test/" folder).

And as a final line the XML-file should be closed with
```xml
 </CRPropa>
```
