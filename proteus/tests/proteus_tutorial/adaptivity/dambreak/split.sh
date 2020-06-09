#!/bin/bash 

#mpirun -np 4 zsplit Reconstructed.dmg Reconstructed.smb 4-Proc/ 4
aprun -n 4 zsplit Reconstructed.dmg Reconstructed.smb 4-Proc/ 4

