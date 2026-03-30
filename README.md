# PhD_Code
Code developed during my PhD at the LDQM - LS-Efetov groups at ICFO and then at LMU. 
All the code posted here are the python functions which I then used in Jupyter Notebooks to do the in depth data analysis. All of these except exportDatatoText, ICElog_extraction and Qcodesdb_to_txt, are meant to be function files, not code which you run direcly. 

- Capacitor_simulation is the code for the electrostatic simulations solving the Poisson equation
- exportDatatoText_v8 converts the hdf5 data given by Labber to a txt.
- Fraunhofer_simulation has very basic way to plot Fraunhofers of JJs and SQUIDs
- homemade_functions has the basic architecture to make nice plots easier. Used more than basic plotting.
- Ic_extraction extracts the critical current of the Fraunhofer patterns using different criteria
- IcB_SimulationFrom_Jcx give different current density distribution profiles
- ICElog_extraction is used to convert the data from the log files of ICE Oxford into usable txt files
- JJ_current_distribution contains the functions to calculate the current denstiy distribution from a Fraunhofer pattern using the Fourier transform
- Qcodesdb_to_txt exports the files of the Qcodes database into txt files as the same as the ones extracted from Labber.
- LDQM_dataplot has all the functions to plot directly the data of the txt files of data to plot in the jupyter notebook.
- Superconductor_analysis functions for extracting critical current from minimums and maximums, make derivatives, plot derivative, and some of the simulations of Fraunhofer patterns with different magnetic signals.
