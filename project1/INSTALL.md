How to start using the code
===========================
*Nicole Potenza, Davi Augusto Mattoso, Miguel Angel Díaz Gutiérrez*

## Before compiling the code

The first step is to clone the directory in your computer, use the following command to do it:

```
jjdjdjdjdj
```

As the code is written in C, a C compiler is necessary (we recommend to use GCC because it was used to compile the code).

The code uses TREXIO an opensource library of quantum chemical data, because of that is necessary, as first step, to install the library HDF5. This last one refers
to Hierarchical Data Format, a software and file format for manage, process and store of heterogenous data (in the following link the documentation is provided: https://support.hdfgroup.org/documentation/).
So, the first step, if Ubuntu is used, is to type in the terminal:
  
  `sudo apt install libhdf5-dev`

  or, if MacOS is used:

  `brew install hdf5`

  The next step is to dowload TREXIO, this can be done by using the following link: https://github.com/TREX-CoE/trexio/releases/download/v2.5.0/trexio-2.5.0.tar.gz . To install TREXIO use:

```
tar -zxvf trexio-2.5.0.tar.gz
cd trexio-2.5.0
./configure
make
sudo make install
```

This will install the library in `/usr/local/lib`.

## Compiling

The following command will link TREXIO in the compilation of the program:

```
gcc -I/usr/local/include -L/usr/local/lib -ltrexio main.c -o hf_mp2_calculation
```

This will generate the program `hf_mp2_calculation`. To run it use the following command:

```
./hf_mp2_calculation
```

By
