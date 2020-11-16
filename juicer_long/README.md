## Juicer directory, contains scripts/, references/  
juicer long pipeline used CPU version of the juicer scripts  
Directories "common" and "juicer_long" should be placed under scripts directory under ${juiceDir}  
cd scripts    
ln -s CPU/common .  
##juicerDir contains scripts directory and juicer_long directory \

## Run commands
t=20  
bash juicer_long_0.1.3.sh\. 
-g hg19 -D  /public/home/yangzhzh/projects/0_aiden_lab/home_juicer -t "$t" -s none -i 100 -r 1 -S early
