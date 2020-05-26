page_numbers=($(seq 0 1 5))

for page_number in "${page_numbers[@]}"
do
./dhusget.sh -d "https://s5phub.copernicus.eu/dhus" -u "s5pguest" -p "s5pguest" -S 2019-04-07T00:00:00.000Z -T "L2__CH4___" -l 100 -P "$page_number" -C "./products_list_$page_number.csv" -o "product" -O "/n/seasasfs02/hnesser/TROPOMI/downloads_20190318/"
done
