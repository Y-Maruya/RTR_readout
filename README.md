This is the codespace for RTR beam experiment.
## Prerequisites

Make sure you have the following installed:

- GCC (g++)
- ROOT framework

## Compilation
Before running the source code, You should compile the source code.  
Of cource, you can use the root interpriter for executing the source code, but I don't check how it works.  
#### To compile the source code ```A.c```, run the following command:

```sh
g++ A.c -o A `root-config --cflags --libs`
```

## Data Aquation and Derivation procedure

### 1. Acquire the validation data of MPPC by pulse laser.
Turn off the ready bottum, Connect the lemo cable to pulse laser inside the beam area.  
And then execute ```./easiroc``` inside ```software2```directly, then ```./udp``` inside ```software2/UDPcontrol2```.  
Set the voltage of MPPC by ```1``` 10 V->20 V-> 30 V-> 40 V-> 50 V step by step, then check the actual voltage ```3```. 
:::note alert
Then write the voltage and to the time to experimental note. 
::: 
:::note alert
絶対に、一回で50V に電圧を上げないでください。!!!!
:::
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

