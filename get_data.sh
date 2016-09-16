#!/bin/bash
# DOWNLOAD files from URLs in a file "urls.txt"
# One url to file-to-download per line.
# will download to folder ./data/{FILENAMEBASE}

while read URL; do
  echo Processing: ${URL}

  # Build string for creating directory
  FILEBASE=$(basename ${URL})
  DIR=./data/${FILEBASE%%.*} # remove all after .

  # Create dir
  echo Creating dir: ${DIR}
  mkdir -p ${DIR}

  # download file to dir
  cd ${DIR} && { curl -O ${URL} ; cd -; }
done <urls.txt
