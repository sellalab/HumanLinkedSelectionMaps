#!/usr/bin/awk -f

substr($5, 3) != $4 && substr($6, 3) != $4{
    print
}