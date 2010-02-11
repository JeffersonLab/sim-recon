#!/usr/bin/env perl

$infile = $ARGV[0];
$outfile = "solenoid_1500_poisson";

# Read in inout file
open(INPUT, $infile);
@lines = <INPUT>;
close(INPUT);

# Open output file for writing
open(OUTPUT, ">$outfile");
print OUTPUT "#\n";
print OUTPUT "# This file generated from the file \"$infile\" via the script poisson2calibDB.pl\n";
print OUTPUT "# Lines begining with \"#--\" below are comments copied from the original file.\n";
print OUTPUT "#\n";

# loop over lines in input file
$my_header_printed = 0;
foreach $line (@lines){
	@tokens = split(/\s+/, $line);
	$ntokens = @tokens;
	$is_comment = 0;
	if($ntokens != 10){
		$is_comment=1;
	}else{
		foreach (@tokens){if(!/[0-9\.\+\-E]/ && length($_)!=0) {$is_comment=1;}} # non-digits indicate comments
	}
	
	if($is_comment){
		print OUTPUT "#-- ".$line;
	}else{
		if($my_header_printed==0){
			print OUTPUT "#\n";
			print OUTPUT "# NOTE: The values from the original file were changed from cm/Gauss to\n";
			print OUTPUT "#       inches/Tesla below (despite any comments above).\n";
			print OUTPUT "#\n";
			print OUTPUT "#     x            y            z            Bx          By           Bz\n";			
			$my_header_printed=1;
		}
	
		$r  = $tokens[1]/2.54;		# convert to inches
		$z  = $tokens[2]/2.54;		# convert to inches
		$Br = $tokens[3]/10000.0;	# convert to Tesla
		$Bz = $tokens[4]/10000.0;	# convert to Tesla
		if($Br==0){$Br = "0.00000000";} # make columns line up better for visual inspection of ASCII file
		
		print OUTPUT "$r \t 0.0000\t $z\t $Br\t 0.0000\t $Bz\n";
	}
}

# Close output file
close(OUTPUT);
