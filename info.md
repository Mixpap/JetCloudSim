### How to use PLUTO 
#### for cloud jet interaction and cooling

1. Copy the PLUTO folder in your home directory (or anywhere you want)
2. In your .bashrc file (usually its in the home directory) add the following:
    * ``export PLUTO_DIR='[your directory of PLUTO]'``
    
    * ``alias ipluto='python $PLUTO_DIR/setup.py'``

    for example:
    
    ```
    export PLUTO_DIR='/home/mixalis/PLUTO'
    alias ipluto='python $PLUTO_DIR/setup.py'
    ```

3. Copy the PLUTO template files in a new folder*
4. Open a terminal in this folder and run the command: 
    * ipluto`
5. In the Setup problem you can change the configuration of the problem 
    * Physics(HD, MHD, RMHD =relativistic MHD etc)
    * DIMENSION and COMPONENTS
    * BODY FORCE
    * COOLING
    * ...
    * USER_DEF_PARAMETERS

6. In the PHYSICS Menu (after we presed ENTER) we keep as it is
7. In the Parameters Menu we see the problem parameters (later) -nothing to do here
8. In the Makefile we choose Linux.gcc.defs for a single core execution, Linux.mpicc.defs for multithreding (you need mpi lib for this)
9. Hit ENTER and press Q to quit to terminal
10. Run the command: ```make```
11. Run the simulation by giving the command 
    * ```./pluto``` (for single thread) 
    * ```mpirun -n 4 pluto``` (for 4 threads for example)
    
### definitions.h
In this file we have the configurations from the `ipluto` script. The important configuration options here are the `USER_DEF_PARAMETERS`. In this option we put how many parameters we are going to give. The paremeter names are given in the `user-defined paremeters (label)` section

### pluto.ini
In this file we have the space-time grid and output configuration. 
* In the `[Grid]` section we give the box of our simulation and the grid density (power of 2 for multithreading)
* In the `[Time]` section we give the Timescale of the whole simulation
* In the `[Boundary]` section we tell to PLUTO we are going to use `useredef` boundary conditions for the `X2-beg` boundary.
* In the `[Static Frid Output]` section we tell PLUTO to gives us a dbl (for binary format output) and a ppm (a simple picture) every 300 timescales
* In the `[Parameters]` section we give our problem parameters (later)

### init.c
In this file we provide the initial conditions