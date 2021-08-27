# this script is used to check the progress of calc_bkgd jobs based on the size of the 1MB "." stdout to log

# chromosome lengths file
lengths=/ifs/data/c2b2/gs_lab/dam2214/linked_selection/data/ch_features/chr_len_all.txt

# check for flags or pattern for scanning the logs
pat='AA_Map'
flag=''

# get kwargs from command line
for var in "${@:1}"; do 
    eval $var
done

PATTERN="${PWD}/*${pat}*"
FLAG=$flag

# if [[ $pat ]] ; then
#     PATTERN="/ifs/data/c2b2/gs_lab/dam2214/run/*"${pat}"*"
#     echo "SUBSET MATCHING: "$PATTERN    
# fi

# process the log files
# grep -H '^\.' $PATTERN returns "FILENAME:....." for all files matching $PATTERN  and containing a "..." progress bar
# updated c returns percentiles instead of MB so progress bars can be summed directly to get the % completion
grep -H '^\.' $PATTERN | awk/parse_pct_bar.awk -v flag=$FLAG | sort -k2n
# grep -H '^\.' $PATTERN | awk/parse_log.awk -v flag=$FLAG $lengths - | sort -k2n

# print usage message
# echo 'OPTIONS: check c shows completed files; check f shows incomplete files; check <glob> shows glob matching files'
echo 'usage: check [flag=FLAG (c=complete f=in_progress)] [pat=GLOB_PATTERN]'
