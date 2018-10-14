#!/bin/sh

#for cpu in $(seq  1 1 8); do 
#    for seed in $(seq  0 1 20); do 
#	./G4dataAnalysis ../Pola0degSeeds/G4data$cpu-$seed.root | grep 'MF' >> 0degSeed.txt
#    done
#done

#for ang in $(seq  0 1 21); do 
    #cp NP0-90deg/M2hitsNP$ang.root M2hitsNP.root
    #./G4dataAnalysis Pola0-90deg/G4data$ang.root | grep 'MF' >> 0-90deg.txt
 #   ./G4dataAnalysis NP0-90deg/G4dataNP$ang.root  > log.txt
    #mv M2hitsNP.root NP0-90deg/M2hitsNP$ang.root
  #  mv D1hitNP.root NP0-90deg/D1hitNP$ang.root
#done

#while read pola alpha 
#do
#    ./G4dataAnalysis ChristofferNBOnaxis/G4dataAlphaPFNB0deg$alpha-$pola.root | grep 'MF' >> AlphaPola.txt
#done < PF-Alpha_sample100.txt

for alpha10 in $(seq  -30 1 10); do 
    alpha=$(echo "($alpha10/10.0)" | bc -l)
    ./G4dataAnalysis Alpha-3.0-1.0/G4dataAlpha$alpha.root | grep 'MF' >> Alpha-3.0-1.0_2.txt
done

exit 0
