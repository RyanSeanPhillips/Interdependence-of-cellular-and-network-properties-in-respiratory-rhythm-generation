#!/bin/bash

################
### Compile code
################
g++ -std=c++0x -O2 spk.cc -o spk


##Spike height examples

#####################
##Example Bursters
#####################
#####Increase in gspk
./spk -DT 0.1 -d .285 0.285 -s 1 -T 40 -w 0.0 -prob 0.0 -kbath 8.5 -gSPK 00 00 -gAHP 0 0 -nap 3.33 3.33 -gleak 3.5 3.5  data/2_0_0.sp -o  >data/2_0_0.hst 2>data/2_0_0.ca&
./spk -DT 0.1 -d .285 0.285 -s 1 -T 40 -w 0.0 -prob 0.0 -kbath 8.5 -gSPK 02 02 -gAHP 0 0 -nap 3.33 3.33 -gleak 3.5 3.5  data/2_0_1.sp -o  >data/2_0_1.hst 2>data/2_0_1.ca&
./spk -DT 0.1 -d .285 0.285 -s 1 -T 40 -w 0.0 -prob 0.0 -kbath 8.5 -gSPK 04 04 -gAHP 0 0 -nap 3.33 3.33 -gleak 3.5 3.5  data/2_0_2.sp -o  >data/2_0_2.hst 2>data/2_0_2.ca&
./spk -DT 0.1 -d .285 0.285 -s 1 -T 40 -w 0.0 -prob 0.0 -kbath 8.5 -gSPK 06 06 -gAHP 0 0 -nap 3.33 3.33 -gleak 3.5 3.5  data/2_0_3.sp -o  >data/2_0_3.hst 2>data/2_0_3.ca

#####Increase in gahp
./spk -DT 0.1 -d .285 0.285 -s 1 -T 40 -w 0.0 -prob 0.0 -kbath 8.5 -gSPK 00 00 -gAHP 00 00 -nap 3.33 3.33 -gleak 3.5 3.5 data/2_1_0.sp -o  >data/2_1_0.hst 2>data/2_1_0.ca&
./spk -DT 0.1 -d .285 0.285 -s 1 -T 40 -w 0.0 -prob 0.0 -kbath 8.5 -gSPK 00 00 -gAHP 6  6  -nap 3.33 3.33 -gleak 3.5 3.5 data/2_1_1.sp -o  >data/2_1_1.hst 2>data/2_1_1.ca&
./spk -DT 0.1 -d .285 0.285 -s 1 -T 40 -w 0.0 -prob 0.0 -kbath 8.5 -gSPK 00 00 -gAHP 12 12 -nap 3.33 3.33 -gleak 3.5 3.5 data/2_1_2.sp -o  >data/2_1_2.hst 2>data/2_1_2.ca&
./spk -DT 0.1 -d .285 0.285 -s 1 -T 40 -w 0.0 -prob 0.0 -kbath 8.5 -gSPK 00 00 -gAHP 18 18 -nap 3.33 3.33 -gleak 3.5 3.5 data/2_1_3.sp -o  >data/2_1_3.hst 2>data/2_1_3.ca



#############
####Example Network Rhythms
#############
###uncoupled
#./spk -dt .025 -DT 20 -d .28 .28 -s 100 -T 80 -w 0.2 -prob 0.13 -kbath 8.5 -gSPK 00 00  -rastonly 0 -gAHP 00 00 data/3_0_0.sp -o  >data/3_0_0.hst 2>data/3_0_0.ca&
#./spk -dt .025 -DT 20 -d .28 .28 -s 100 -T 80 -w 0.2 -prob 0.13 -kbath 8.5 -gSPK 15 15  -rastonly 0 -gAHP 00 00 data/3_1_0.sp -o  >data/3_1_0.hst 2>data/3_1_0.ca&
#./spk -dt .025 -DT 20 -d .28 .28 -s 100 -T 80 -w 0.2 -prob 0.13 -kbath 8.5 -gSPK 00 00  -rastonly 0 -gAHP 30 30 data/3_2_0.sp -o  >data/3_2_0.hst 2>data/3_2_0.ca
#/ip3 -dt .025 -DT 20 -d .28 .28 -s 100 -fc 1.0 1.0 -T 150 -w 0.2 -fw 1 1 -prob 0.13 -kbath 8.5  -fsca 0.0 0.0 -gSPK 00 00 -gAHP 00 00 -rastonly 0 data/fig2/v2/sup/0_4_1.sp -o  >data/fig2/v2/sup/0_4_1.hst 2>data/fig2/v2/sup/0_4_1.ca&

