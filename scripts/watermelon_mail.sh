#!/bin/bash

while getopts "c:f:m:s:t:" opt; do
  case $opt in
    c) 
      CC_FLAG_ADDRESS=" -c ${OPTARG}"
      ;;
    f)
      mail_from=$OPTARG
      ;;
    m)
      mail_message=$OPTARG
      ;;
    s)
      mail_subject=$OPTARG
      ;;
    t)
      mail_to=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
    :)
        echo "Option -$OPTARG requires an argument" >&2
        exit 1
        ;;
  esac
done

export EMAIL=${mail_from}; echo -e "${mail_message}" | mutt -s "${mail_subject}" ${CC_FLAG_ADDRESS} ${mail_to}
