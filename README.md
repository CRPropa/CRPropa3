CRPropa
========

![stable release](https://img.shields.io/badge/stable\_release-3.2.1-darkblue)
[![Build status](https://github.com/crpropa/crpropa3/actions/workflows/testing.yml/badge.svg)](https://github.com/crpropa/crpropa3/actions/)
[![Average time to resolve an issue](https://isitmaintained.com/badge/resolution/CRPropa/CRPropa3.svg)](https://isitmaintained.com/project/CRPropa/CRPropa3)
[![Percentage of issues still open](https://isitmaintained.com/badge/open/CRPropa/CRPropa3.svg)](https://isitmaintained.com/project/CRPropa/CRPropa3)

[![DOI:10.1088/1475-7516/2022/09/035](http://img.shields.io/badge/DOI-10.1088/1475-7516/2022/09/035.svg)](<https://doi.org/10.1088/1475-7516/2022/09/035>)
[![arXiv](https://img.shields.io/badge/arXiv-2208.00107-b31b1b.svg)](https://arxiv.org/abs/2208.00107)
[![ascl:2208.016](https://img.shields.io/badge/ascl-2208.016-blue.svg?colorB=262255)](https://ascl.net/2208.016)

CRPropa is a publicly available simulation framework to study the propagation
of ultra-high-energy nuclei up to iron on their voyage through an
(extra)galactic environment. It takes into account: pion production,
photodisintegration and energy losses by pair production of all relevant
isotopes in the ambient low-energy photon fields as well as nuclear decay.
CRPropa can model the deflection in (inter)galactic magnetic fields, the
propagation of secondary electromagnetic cascades and neutrinos for a multitude
of scenarios for different source distributions and magnetic environments. It
enables the user to predict the spectra of UHECR (and of their secondaries),
their composition and arrival direction distribution. Additionally, the
low-energy Galactic propagation can be simulated by solving the transport
equation using stochastic differential equations. CRPropa features a very
flexible simulation setup with python steering and shared-memory
parallelization.


## Interactive Online Demo
You can try out CRPropa online at [vispa.physik.rwth-aachen.de](https://vispa.physik.rwth-aachen.de/).
Use the guest login and go to the CRPropa example via "VISPA Cluster" --> "Open Examples".


## Installation and Documentation
To install CRPropa, download and unzip either the

* [latest release](https://github.com/CRPropa/CRPropa3/releases/latest) (recommended),
* or [current development version](https://github.com/CRPropa/CRPropa3).

Installation instructions, usage examples  and API documentation can be found on the [documentation web site of
CRPropa](https://crpropa.github.io/CRPropa3/).


## Support
Please use the [ticket system](https://github.com/CRPropa/CRPropa3/issues) for
support and in case of general questions. Please browse also the
[documentation](https://crpropa.github.io/CRPropa3/) and previous support requests on
[installation](https://github.com/CRPropa/CRPropa3/issues?utf8=%E2%9C%93&q=is%3Aissue+label%3Ainstallation+)
and
[usage](https://github.com/CRPropa/CRPropa3/issues?utf8=%E2%9C%93&q=label%3Ausage-question+)
of CRPropa before opening a new ticket.

To receive announcements etc., please subscribe to our mailing list by sending
a mail with subject: subscribe crpropa-user to sympa@desy.de from the address
you wish to subscribe.

## How to cite CRPropa
If you use CRPropa 3.2 for your research, please cite

**JCAP (2022) no. 09, 035; [arXiv:2208.00107](https://arxiv.org/abs/2208.00107)**

as well as [additional publications](https://crpropa.github.io/CRPropa3/pages/howto_cite_crpropa.html) dependent on the components you are using.


## Publications based on CRPropa
An extensive list of publications using CRPropa can be found via
[inSPIRE](https://inspirehep.net/search?ln=en&ln=en&p=refersto%3Arecid%3A1322902+or+refersto%3Arecid%3A1432676+or+refersto%3Arecid%3A1242078&of=hb&action_search=Search&sf=earliestdate&so=d&rm=&rg=25&sc=0).


## Plugins
Plugins are extensions of the core CRPropa framework, but they are not maintained by the CRPropa developer team. Instructions to install plugins can be found in the [documentation](https://crpropa.github.io/CRPropa3/pages/example_notebooks/extending-CRPropa/extending-CRPropa.html#Plugins:-Integrate-Custom-C++-Code-to-CRPropa%E2%80%99s-Python-Steering).

Make sure to correctly cite the plugins when using them.

| Name | Purpose | Link |
| ---- | ------- | ---- |
| FieldlineIntegrator | Magnetic Field Analysis | <https://github.com/lukasmerten/CRPropa_FieldLineIntegrator> |
| grplinst | Plasma Instabilities | <https://github.com/rafaelab/grplinst> |
| monopole | Magnetic Monopole Studies | https://github.com/chchristie/monopole/tree/main |
| ROOTOutputPlugin | Output into root file format | https://github.com/CRPropa/ROOTOutputPlugin |
