This is the codespace for RTR beam experiment.
## Prerequisites

Make sure you have the following installed:

- GCC (g++)
- ROOT framework

## Compilation
Before running the code, You should compile the code.  
Of cource, you can use the root interpriter for executing the code, but I don't check how it works.  
#### To compile the code ```A.c```, run the following command:

```sh
g++ A.c -o A `root-config --cflags --libs`
```

## Data Aquation and Derivation procedure

### 1. Acquire the validation data of MPPC by pulse laser.
### 2. Check the data by ```readout/val_v3.c```.
Usage is 
```sh
./val_v3 "path of validation data"
``` 
Please ensure all of the fit goes well and check how many channel of MPPC is not working by checking pdf or root.
### 3. Acquire the dark-noise validation data of MPPC by clock generator. 
### 4. Check the data by ```readout/dark.c``` 
Usage is 
```sh
./dark "path of dark-noise validation data"
``` 
Please ensure all of the fit goes well and check how many channel of MPPC is not working by checking pdf or root.
### 5. Acquire the run-data.
### 6. After the run, repeat 1~4.
### 7. Write config file.
config file format is like ```readout/config```
### 8. Make TTree of our run data by ```readout/rawdata_tree.C```
Usage is 
```sh
./rawdata_tree "path of config file"
``` 
root file is ommitted to the path where inputted easiroc raw data exists.

