#!/bin/bash

while getopts "f:m:s:t:" opt; do
  case $opt in
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

if [[ ${WATERMELON_EMAIL_ENABLED+x} ]]; then
    EMAIL=${mail_from};
    echo -e "${mail_message}" | mutt -s \"${mail_subject}\" ${mail_to}
fi