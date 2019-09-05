#!/bin/sh

#-------------------------------------------------------------------------------
# Author: Thomas Uchdorf, t.uchdorf@mpie.de
#
# Description: This little helper script generates files that serve for tests 
#              of the file system functionality of SPHInX.
#
# Usage: ./createTestResults.sh <command> <filename for the output>
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# init ()
#-------------------------------------------------------------------------------
init ()
{
   param1=$1;
   param2=$2; 
   # Using relative targets for the links
   target="${param2##*/}Target" 

   # Doing nothing in case of "Nothing"
   if [[ ${param1} = "File" ]]
   then
      echo "I was born in ${param2}">${param2}
   elif [[ ${param1} = "Dir" ]]
   then
      mkdir ${param2};
   elif [[ ${param1} = "FileSymLink" ]]
   then
      echo "I was born in ${target}">"${param2}Target"
      ln -sf ${target} ${param2}
   elif [[ ${param1} = "DirSymLink" ]]
   then
      mkdir "${param2}Target"
      ln -sf ${target} ${param2}
   elif [[ ${param1} = "EmptySymLink" ]]
   then
      # Willingly skipping the target creation
      ln -sf ${target} ${param2}
   fi
}

#-------------------------------------------------------------------------------
# determinePreAndPost ()
#-------------------------------------------------------------------------------
determinePreAndPost ()
{
   param1=$1;
   if [[ -L ${param1} ]]
   then
      if [[ -f ${param1} ]]
      then
         return 3
      elif [[ -d ${param1} ]]
      then
         return 4
      else
         return 5
      fi
   elif [[ -f ${param1} ]]
   then
      return 0
   elif [[ -d ${param1} ]]
   then
      return 1
   else
      return 2
   fi

   echo "Error: Could not stat ${param1}.">&2
   exit 4
}

#-------------------------------------------------------------------------------
# main ()
#-------------------------------------------------------------------------------

operation="cp_a"
outfile="default.outp"

if [[ ! -z ${1} ]]
then
   operation=${1}
else
   echo "Warning: As no operation was specified as first parameter, ${operation} is used as default one."
fi

if [[ ! -z ${2} ]]
then
   outfile=${2}
else
   echo "Warning: Since no output file was specified as second parameter, ${outfile} is used as default one."
fi

# Writing a header to the output file
echo "operation srcType destType result postSrcType postDestType">${outfile}

srcDir="testDir"
destDir="testDir"

# Checking whether the precondition for a successful procedure are fullfilled
if [[ -d ${srcDir} ]]
then
echo "Error: ${srcDir} already exists">&2
exit 1
fi
if [[ -d ${destDir} && ${srcDir} != ${destDir} ]]
then
echo "Error: ${destDir} already exists">&2
exit 2
fi

# Creating the directories for the operation simulations
mkdir ${srcDir}
if [[ ${srcDir} != ${destDir} ]]
then
   mkdir ${destDir}
fi

ARR=("File" "Dir" "Nothing" "FileSymLink" "DirSymLink" "EmptySymLink" "Self")
         
# Initializing the number that will be used to create unique identifiers
i=0
# Iterating over all types
for src in ${ARR[*]}
do
   # Excluding the special case "Self" as source type
   if [[ ${src} != "Self" ]]
   then
      # Traversing all types
      for dest in ${ARR[*]}
      do
         #echo "${i}: ${src} to ${dest}";
         # Creating unique identifiers
         let srcId=i
         let destId=i+7*6
         # Performing a special treatment in case the type is "Self"
         if [[ ${dest} = "Self" ]]
         then
            curSrc="${srcDir}/${src}${srcId}"
            curDest="${srcDir}/${src}${srcId}"
            init ${src} ${curSrc}
         else
            curSrc="${srcDir}/${src}${srcId}"
            curDest="${destDir}/${dest}${destId}"
            init ${src} ${curSrc}
            init ${dest} ${curDest}
         fi

         # Printing the current test to the standard output 
         echo "${i}: ${curSrc} to ${curDest}";
         # Checking which operation is to be tested
         if [[ ${operation} = "cp_a" ]]
         then
            # Executing the actual test
            `cp -a ${curSrc} ${curDest} 2>/dev/null`
         elif [[ ${operation} = "mv" ]]
         then
            # Executing the actual test
            `mv ${curSrc} ${curDest} 2>/dev/null`
         else
            echo "Error: An unknown operation was specified as second parameter.">&2
            exit 3;
         fi
         # Storing the actual result
         actResult=$?

         determinePreAndPost ${curSrc}
         postSrcType=$?
         determinePreAndPost ${curDest}
         postDestType=$?

         # Writing all relevant test parameters to a logfile
         echo "${curSrc} ${curDest} ${actResult} ${postSrcType} ${postDestType}">>${outfile};

         # Incrementing the number used to create unique identifiers
         let i=i+1
      done
   fi
done

# Cleaning up
#rm -r ${srcDir} 

# Returning that the program has successfully finished
exit 0
