#!/usr/bin/perl

for my $i (1..64){
	printf "#define T$i 0x%08x\n",(int (4294967296*abs(sin($i))))."\n";
}
