Interrupting simulations on runtime
------------------------------------------------

CRPropa simulations can be interrupted on runtime with the `SIGTERM` or `SIGINT` signals. 
If the user defines an output for the interruption (called `InterruptAction`) all candidates which are currently in the simulation will be passed to this output. 
In the error stream the user will see a message denoting the number of candidates which have not been started yet. 
If the simulation was run with a `candidateVector` as source, the indices of the candidates which have not been started yet will be printed or written to the file.
For a simulation with a source interface, a restart with the missing number of candidates will be sufficient to continue the simulation.

.. toctree:: 
   :caption: Using a candidateVector as source
   :maxdepth: 1
   
   example_notebooks/interrupting_simulations/interrupt_candidateVector.ipynb

.. toctree:: 
   :caption: Using a source interface
   :maxdepth: 1

   example_notebooks/interrupting_simulations/interrupt_source.ipynb


