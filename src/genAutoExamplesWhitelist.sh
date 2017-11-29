#!/bin/sh

##########################################################################
# Shellscript:	read pair styles from liggghts in files
# Author     :	Josef Kerbl <josef.kerbl@dcs-computing.com>
# Date       :	2017-05-17
# Category   :	LIGGGHTS Utility
##########################################################################
# Description
#    searches a path for liggghts files and reads the pair style
##########################################################################

PN=$(basename "$0") # Program name
VER='0.1'

Usage () {
    echo >&2 "$PN - searches paths recursively for LIGGGHTS example input scripts (in.*), and extracts the pair styles, version $VER
usage: $PN [-o outputfile] [-s LIGGGHTS src directory path] [-i] [-v] [-f] [-a] [paths|files]
    -o: output filename (default: ./style_contact_model_autoExamples.whitelist)
    -s: specify LIGGGHTS src directory path (default: ./)
    -i: do not search recursively, look in specified files
    -a: search recursive in all files, not only in.*
    -f: force outfile to be overwritten
    -v: verbose

Example:
    $PN ../examples/LIGGGHTS"
    exit 1
}

Msg () { echo >&2 "$PN: $*"; }
Fatal () { Msg "$@"; exit 1; }

ProcessFile ()
{
      file="$1"
      if [ "$verboseFlag" = "true" ]; then echo "processing $file"; fi
      style=$(grep 'pair_style' "$file")

      #read line by line
      printf '%s\n' "$style" | while IFS= read line ; do

        ### remove comments from line
        line=$(echo "$line" | grep -o '^[^#]*')
        line=$(echo "$line" | sed 's/"//g')
        wStyle="GRAN_MODEL("
        #do checking on style
        if [ -n "$line" ]; then
          if [ "$verboseFlag" = "true" ]; then echo "found: $line"; fi
          #check normal model
          i=0
          tmp=""
          for model in $normal_models; do

            if echo "$line" | grep -q "\b$( eval "echo \"\$normal_ident_$i\"")\b" ; then
              tmp=${model}
              break
            fi
            i=$((i + 1))
          done
          if [ -n "$tmp" ]; then
            wStyle="${wStyle}${tmp}, "
          else
            wStyle="${wStyle}HERTZ, "
          fi

          #check tangential model
          i=0
          tmp=""
          for model in $tangential_models; do
            if echo "$line" | grep -q "\b$( eval "echo \"\$tangential_ident_$i\"")\b" ; then
              tmp=${model}
              break
            fi
            i=$((i + 1))
          done
          if [ -n "$tmp" ]; then
            wStyle="${wStyle}${tmp}, "
          else
            wStyle="${wStyle}TANGENTIAL_OFF, "
          fi

          #check cohesion model
          i=0
          tmp=""
          for model in $cohesion_models; do
            if echo "$line" | grep -q "\b$( eval "echo \"\$cohesion_ident_$i\"")\b" ; then
              tmp=${model}
              break
            fi
            i=$((i + 1))
          done
          if [ -n "$tmp" ]; then
            wStyle="${wStyle}${tmp}, "
          else
            wStyle="${wStyle}COHESION_OFF, "
          fi

          #check rolling model
          i=0
          tmp=""
          for model in $rolling_models; do
            if echo "$line" | grep -q "\b$( eval "echo \"\$rolling_ident_$i\"")\b" ; then
              tmp=${model}
              break
            fi
            i=$((i + 1))
          done
          if [ -n "$tmp" ]; then
            wStyle="${wStyle}${tmp}, "
          else
            wStyle="${wStyle}ROLLING_OFF, "
          fi

          #check surface model
          i=0
          tmp=""
          for model in $surface_models; do
            if echo "$line" | grep -q "\b$( eval "echo \" \$surface_ident_$i\"")\b" ; then
              tmp=${model}
              break
            fi
            i=$((i + 1))
          done
          if [ -n "$tmp" ]; then
            wStyle="${wStyle}${tmp})"
          else
            wStyle="${wStyle}SURFACE_DEFAULT)"
          fi


          ###### check if model exists in whitelist
          if [ -z "$(grep "${wStyle}" "$outfile")" ] ; then #"$result2" ] ; then
            echo "$wStyle" >> "$outfile"
            if [ "$verboseFlag" = "true" ]; then echo "extracted unique style $wStyle"; fi
          fi
        fi
      done
}

#defaults
recflag="true"
outfile="style_contact_model_autoExamples.whitelist"
ligSrcPath="./"
forceFlag="false"
verboseFlag="false"
allFlag="false"

while getopts ahifvqo:s:: opt
do
    case "$opt" in
      o)  outfile="$OPTARG";;
      s)  ligSrcPath="$OPTARG";;
      i)  recflag="false";;
      f)  forceFlag="true";;
      a)  allFlag="true";;
      v)  verboseFlag="true";;
      h)  Usage;;
      \?) Usage;;
    esac
done
shift $(( OPTIND - 1 ))

## if no argument, show help
if [ $# -lt 1 ] ; then
  Usage
fi

##### DO STUFF
tmpfile=filelist.tmp

if [ -e "$outfile" ]; then
  if [ $forceFlag = "false" ]; then
    Fatal "outfile exists: $outfile. Run with -f to overwrite"
  else
    rm "$outfile"
  fi
fi
touch "$outfile"

if [ -e "$tmpfile" ]; then
  if [ $forceFlag = "false" ]; then
    Fatal "tmpfile exists: $tmpfile. Run with -f to overwrite"
  else
    rm "$tmpfile"
  fi
fi
touch "$tmpfile"

sed_ex="sed -E" # BSD sed
sed --version | grep -i gnu > /dev/null 2>&1
[ $? -eq 0 ] && sed_ex="sed -r" # GNU sed

surface_models=$( grep -s -E '^SURFACE_MODEL' "$ligSrcPath"/surface_model_*.h | $sed_ex 's/.*SURFACE_MODEL\((.+),\s*(.+),\s*(.+)\)/\1/' )
normal_models=$( grep -s -E '^NORMAL_MODEL' "$ligSrcPath"/normal_model_*.h | $sed_ex 's/.*NORMAL_MODEL\((.+),\s*(.+),\s*(.+)\)/\1/' )
tangential_models=$(grep -s -E '^TANGENTIAL_MODEL' "$ligSrcPath"/tangential_model_*.h | $sed_ex 's/.*TANGENTIAL_MODEL\((.+),\s*(.+),\s*(.+)\)/\1/')
cohesion_models=$(grep -s -E '^COHESION_MODEL' "$ligSrcPath"/cohesion_model_*.h | $sed_ex 's/.*COHESION_MODEL\((.+),\s*(.+),\s*(.+)\)/\1/')
rolling_models=$(grep -s -E '^ROLLING_MODEL' "$ligSrcPath"/rolling_model_*.h | $sed_ex 's/.*ROLLING_MODEL\((.+),\s*(.+),\s*(.+)\)/\1/')

surface_ident=$(grep -s -E '^SURFACE_MODEL' "$ligSrcPath"/surface_model_*.h | $sed_ex 's/.*SURFACE_MODEL\((.+),\s*(.+),\s*(.+)\)/\2/')
normal_ident=$(grep -s -E '^NORMAL_MODEL' "$ligSrcPath"/normal_model_*.h | $sed_ex 's/.*NORMAL_MODEL\((.+),\s*(.+),\s*(.+)\)/\2/')
tangential_ident=$(grep -s -E '^TANGENTIAL_MODEL' "$ligSrcPath"/tangential_model_*.h | $sed_ex 's/.*TANGENTIAL_MODEL\((.+),\s*(.+),\s*(.+)\)/\2/')
cohesion_ident=$(grep -s -E '^COHESION_MODEL' "$ligSrcPath"/cohesion_model_*.h | $sed_ex 's/.*COHESION_MODEL\((.+),\s*(.+),\s*(.+)\)/\2/')
rolling_ident=$(grep -s -E '^ROLLING_MODEL' "$ligSrcPath"/rolling_model_*.h | $sed_ex 's/.*ROLLING_MODEL\((.+),\s*(.+),\s*(.+)\)/\2/')

## create pseudo arrays for access
i=0
for ident in $surface_ident; do
  eval surface_ident_$i="$ident"
  i=$((i + 1))
done

i=0
for ident in $normal_ident; do
  eval normal_ident_$i="$ident"
  i=$((i + 1))
done

i=0
for ident in $tangential_ident; do
  eval tangential_ident_$i="$ident"
  i=$((i + 1))
done

i=0
for ident in $cohesion_ident; do
  eval cohesion_ident_$i="$ident"
  i=$((i + 1))
done

i=0
for ident in $rolling_ident; do
  eval rolling_ident_$i="$ident"
  i=$((i + 1))
done

if [ "$recflag" = "true" ]
then
  for path in "$@"; do
    if [ "$allFlag" = "false" ]
    then
      find "$path" -name 'in.*' > "$tmpfile"
    else
      find "$path" -type f > "$tmpfile"
    fi

    while read line; do
      ProcessFile "$line"
    done < "$tmpfile"
  done
else
  for file in "$@"; do
    ProcessFile "$file"
  done
fi

rm "$tmpfile"
####
# final output
echo "parsing of input scripts completed"
echo "they are added to the whitelist: $outfile , only \"style_contact_model_autoExamples.whitelist\" is automatically invoked"
