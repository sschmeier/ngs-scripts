#!/bin/bash
# USAGE: /bin/bash rsync_dir.sh dir-to-sync dir-to-sync-to
#
# Convenience script for doing a rsync call of one dir to another dir, 
# with some error logging.
# Use absolute paths when calling the script.
# 
# dir-to-sync   : dir or user@remote:dir 
# dir-to-sync-to: dir or user@remote:dir 
# 


# date_time
now=$(date +"%Y-%m-%d_%H%M%S");

# Log-files
FILE2=rsync_${now}.log;
FILE3=rsync_${now}.err;

echo '----------------------------' > $FILE2;
echo `date` >> $FILE2;
echo '----------------------------' >> $FILE2;
rsync -av --no-links --delete --exclude "*~" --exclude ".git" $1 $2 >> $FILE2 2> $FILE3;
echo '----------------------------' >> $FILE2;
