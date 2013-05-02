#!/bin/bash

# check if the executable wants to dynamically import any of the malloc/free
# related symbols.

for i in $@; do
   
	nm $i | egrep -e "U malloc$" -e "U memalign$" -e "U posix_memalign$" -e "U free$" -e "U realloc$" -e "U calloc$"   

    if [ $? -eq "0" ]; then
	echo "ERROR: $i: memory allocation without rax_ prefix detected."	
    fi
done

exit 0;


