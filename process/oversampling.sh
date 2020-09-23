#!/bin/csh -f
# @(#) dindex

#-----------------------------------------------------------------
# Set path
#-----------------------------------------------------------------
# set Input_Dir   = "/n/seasasfs02/hnesser/TROPOMI/oversampling_input_csvs/"
set Input_Dir   = "/n/holyscratch01/jacob_lab/hnesser/TROPOMI/oversampling/input_csvs"
set Output_Dir  = "/n/holyscratch01/jacob_lab/hnesser/TROPOMI/oversampling/output_csvs"
# set Output_Dir  = "/n/seasasfs02/hnesser/TROPOMI/oversampling_output_csvs/"

# Make the savedir, if necessary
if ( ! -d $Output_Dir ) mkdir -p $Output_Dir

#-----------------------------------------------------------------
# Output resolution you want
#-----------------------------------------------------------------
set Res = 0.01

# set filenames = ("201909_latlim.csv")
set filenames = (`find $Input_Dir -name "*latlim.csv"`)

foreach Input_Filename ($filenames)
#-----------------------------------------------------------------
# Set file names
#-----------------------------------------------------------------
set Input_Filename = $Input_Filename
set addition = "_oversampled"
set Output_Filename = "$Input_Filename$addition"
set inin = $Input_Dir$Input_Filename
echo $inin
echo $Output_Filename

#-----------------------------------------------------------------
# Call RegridPixel.x, and pass user inputs
#-----------------------------------------------------------------
./RegridPixels.x<<EOF
$Output_Dir
$inin
$Output_Filename
$Res
EOF

quit:
end

exit
