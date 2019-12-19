```shell
#!/bin/bash

#-------------------------------------------------------------------------------

#This script swap between gaps and straight lines "|" (easy to modify and to use other characters)

#-------------------------------------------------------------------------------                     

DIRINPUTS=**_[enter_your_path]_**/scripts/test   ####Set the path to your input files working directory

#####Set working directory to the input files directory
cd $DIRINPUTS

#####Create path to store the output files
mkdir swap_outputs
DIROUTPUTS=swap_outputs

#####list the file sin the input directory
ls *.fastq > reads-files.txt

#####Swap gaps-straight lines within the files
while read FILE; do
    
  sed 's/ /|/g' $FILE > $DIROUTPUTS/swaped_$FILE
  if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi
  
done < reads-files.txt

#-------------------------------------------------------------------------------

# End
echo '*********************END*****************************'
exit 0

#-------------------------------------------------------------------------------
```
