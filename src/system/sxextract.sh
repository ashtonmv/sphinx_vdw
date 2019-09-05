#!/bin/sh
#----------------------------------------------------------------------------
# Extract input files directly out of *.cpp source files. The following
# tags are used:
#
# //= filename
# //+ line
# 
# Example: extract the following lines and save them as 'myfile.sx'
#
# //= myfile.sx
# //+ include <parameters.sx>;
# //+
# //+ structure  {
# //+    ...
# //+ }
#
#----------------------------------------------------------------------------
#
# 04/07/03   Sixten Boeck
#            boeck@sfhingx.de
#
#----------------------------------------------------------------------------

# --- print only a single line of a file
#     usage:
#        getLine /etc/passwd  20
#
#     prints the 20th line of the file /etc/passwd
getLine()
{
   file=$1
   line=$2
   res=`head -$line $file | tail -1`
   echo "$res"
}



# --- if invoked with -clean remove all files specified with the //= token
cleaning='no'
argc=$#
if [ $argc = 1 ]; then
   if [ $1 = "-clean" ]; then
      cleaning='yes'
   fi
fi

# --- loop over all *.cpp files
infiles=`ls *.cpp`
for infile in $infiles; do

   if [ ! -e $infile ]; then
      echo "$infile not found"
      exit 1
   fi

   outfiles=`cat $infile  | grep    //= | sed s/\\\/\\\/=//g | sed s/\ //g`
   for outfile in $outfiles; do

      echo "parsing $infile..."

      if [ $cleaning = 'no' ]; then


         # --- remove file if it exists already
         test -e $outfile && rm $outfile

         # --- do we have to create a subfolder first?
         folder=`echo $outfile \
               | awk -F '/' '{for (i=1;i<NF;i++) a=a $i "/"; print a; }'`
         if [ ! -z $folder ]; then
            echo "   creating folder $folder"
            mkdir -p $folder
         fi

         # --- create output file and write some header information
         echo "// This file has been extracted automatically from file" \
              >> $outfile
         echo "// $infile" >> $outfile
         echo >> $outfile

     
         # --- if is possible to define more than 1 one file
         #     extract output between
         #        //= the_matching_filename
         #        //+ ...
         #        //= 
         #     or
         #        //= the_matching_filename
         #        //+ ...
         #        EOF
         echo "   generating file $outfile"
         linenr=`cat $infile  | grep -n //= | grep $outfile | awk -F ':' '{print $1}'`
         nlines=`wc -l $infile | awk '{print $1}'`
         iline=$linenr
         line1=$iline
         while [ $iline -le $nlines ]; do
            iline=`expr $iline + 1`
            contents=`head -$iline $infile | tail -1` 
# extract the iline-th line
            newtoken=`echo $contents | grep -v //+`
            if [ ! -z "$newtoken" ]; then
               break;
            fi
            echo $contents | grep //+ | sed s/\\/\\/+//g | sed s/^\ //g >> $outfile
         done
         line2=$iline
         echo "   extracted lines $line1-$line2"

      else

         echo "   removing file $outfile"
         rm $outfile

      fi

   done

done
