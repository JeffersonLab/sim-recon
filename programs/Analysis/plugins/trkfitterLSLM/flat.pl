#!/usr/bin/env perl
$dist = 800;
$delta = 2000;
for ($pos = -1000.0; $pos < 1000.0; $pos += 10.0) {
    $doca = abs($pos);
    $resid = $doca - $dist;
    $resid2 = $resid*$resid + 2*$delta;
    if ($doca > $dist) {
	$resid3 = $resid2 - $delta;
    } else {
	$resid3 = $delta;
    }
    print "$pos $resid $resid2 $resid3\n";
}
exit;
