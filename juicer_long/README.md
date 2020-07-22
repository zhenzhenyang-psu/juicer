## Juicer directory, contains scripts/, references/
#juicer long pipeline used CPU version of the juicer scripts
#common folder and juicer_long folder should be placed under scripts directory under ${juiceDir} 
cd scripts
ln -s CPU/common .
ln -s juicer_long . 
