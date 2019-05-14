### Compile the unit tests 

To compile the code, you can either do one of the two things:

#### Manual
- Change to the lammps5Nov16/src/STUBS directory, then build dummy MPI lib:

   $ make

- Change to the lammps5Nov16/src directory, then build Lammps code (donot forget to install nufeb package):

   $ make serial_nufeb

   $ make serial_nufeb mode=lib

- Change to the unittests/GTest directory, then build google test lib:

   $ make 

- In the main directory, build all tests:

   $ make 

#### Automated
- In the main directory, build everything above by running:

   $ ./BuildLibrarys

### Running the unit tests

- You can either test each single function by running:

   $ ./test_name

- Or run all tests:

   $ ./run.sh
