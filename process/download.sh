#!/bin/bash
#################### DOWNLOAD TROPOMI ####################
# This is a script that calls the SRON provided dhusget.sh
# script to download TROPOMI data.
# Reference: https://scihub.copernicus.eu/userguide/BatchScripting
#
# Download options:
# -d : URL of the Data Hub Service to be polled
# -u : username
# -p : password
# -m : mission name <EXCLUDED>
# -i : instrument name <EXCLUDED>
# -t : search for products ingested in the last x hours
#      <EXCLUDED>
# -s : ingestion date from <EXCLUDED>
# -e : ingestion date to <EXCLUDED>
# -S : sensing date from
# -E : sensing date to
# -f : search for products ingested after date and time
#      provided in file <EXCLUDED>
#      ** CHANGE TO THIS OPTION
# -c : coordinates (format: lon1,lat1:lon2,lat2)
# -T : product type
# -l : maximum number of results per page (max: 100)
# -P : page number
# -C : write the list of products in a specified CSV
# -o : download the 'manifest', 'product', or 'all'
# -O : save the product ZIP file (not actually a zip) in
#      the specified folder with the specified filename
##########################################################

# Define the directories
#download_dir="/n/holyscratch01/jacob_lab/hnesser/TROPOMI"
#download_subdir="${download_dir}/downloads"
#script_dir="/n/home04/hnesser/TROPOMI"
#last_download_name="lastdownload"
#last_download_dir="/n/seasasfs02/hnesser/TROPOMI"
download_dir=$1
download_subdir=$2
download_subdir="${download_dir}/${download_subdir}"
script_dir=$3
last_download_dir=$4
last_download_name=$5

# dhusget.sh options
d="https://s5phub.copernicus.eu/dhus"
u="s5pguest"
p="s5pguest"
f="${last_download_dir}/${last_download_name}"
T="L2__CH4___"
l=100
q="OSquery-result.xml"

# cd into directory where data will be downloaded
cd ${download_dir}

# Copy the needed scripts
cp ${script_dir}/process/dhusget.sh ./
cp ${script_dir}/python/download_page_count.py ./

# Save the original start date as a variable.
# We will use this start date for the remainder of
# the downloads
start_date=$(head -1 ${f})

# Set the end date to the first day of last month
# (giving time for lag in retrievals)
end_date=$(date -d "-1 month" +"%Y-%m-01T00:00:00.000Z")

# Cancel the job if start_date == end_date
if start_date == end_date; then
    echo "Start date is the same as end date."
    exit 1
fi

# Define a function to download a page of data
# any time after that start date
download_page()
{
    i=$1
    ./dhusget.sh -d "${d}" -u "${u}" -p "${p}" -S "${start_date}" -E "${end_date}" -T "${T}" -l "${l}" -P "${i}" -C "./products_list_$i.csv" -o "product" -O "${download_subdir}" -R "failed_MD5_$i.csv" > log.log
}

# We start by downloading the first page
echo "Downloading files after start date ${start_date} and before ${end_date}"
echo "Downloading page 1"
download_page 1

# This gives us a file query_result.xml that we can
# parse in order to find the total number of results
source activate invpy
page_count=`python download_page_count.py ${q} ${l}`
conda deactivate

# Now create an array from 2 to page_count
page_numbers=($(seq 30 1 ${page_count}))
#page_numbers = ($(seq 4 -1 2))

# Then download the results from each of those pages.
for page_number in "${page_numbers[@]}"
do
    echo "Downloading page ${page_number}/${page_count}"
    download_page ${page_number}
done

# Finally, update the last download date
echo "${end_date}" > ${f}
