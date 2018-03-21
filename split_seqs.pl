#!/usr/bin/perl
#Created by Nicolas Brodu
#Modified by David Chen
#Purpose: takes a large number of unsorted image files generated from many "$cycles" or repetitions of "$scans" containing several images, "$imperscan,"
#from a directory/directories and makes symbolic links to them sorted in separated directories

use POSIX;
use File::Path "mkpath";

# experiment parameters:
# – 20 cycles of 60 scans => 1200 scans for this sequence
# – each scan comprises 360 images => 432000 images for this sequence

$cycles = 2;
$scans = 10;
$seqscans = $cycles * $scans;
$imperscan = 360;
$images = $imperscan * $seqscans;

for $d (glob("GLOBAL_FILEPATH/*")) {
    print "Processing directory $d\n";
    for $f (glob("IMAGE_FILEPATH/*")) {
        if ($f =~ /.*_(\d+).tiff/) {
            # use 0-based arithmetic, much simpler for modulo ops
            $idx = $1 - 1;#change this - number depending on the first picture in the data
            $seq = floor($idx/$images);
            $idx_within_seq = $idx - $seq * $images;
            $scan = floor($idx_within_seq / $imperscan);
            $slice = $idx_within_seq - $scan * $imperscan;
            $idx += 1; $scan += 1; $slice += 1;
            $seq += 1; # specific to this experiment
            #print "$idx $seq $scan $slice => $f\n";
            mkpath("OUTPUT_FILEPATH/scan_$scan");
            symlink("../../$f","OUTPUT_FILEPATH/scan_$scan/slice_$slice.tiff");
        }
    }
}
