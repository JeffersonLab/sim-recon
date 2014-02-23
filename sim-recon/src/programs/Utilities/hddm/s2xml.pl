#!/usr/bin/perl -w

# purpose: use the perl script schema-hddm to generate from the schema
#          file event.xsd the two hddm data format file hits.xml and
#          reaction.xml which are needed by file event.xml
#
# Author: B. Zihlmann Sep. 2007
#


use strict;
use FileHandle;

sub Usage(){

    printf "\n Usage: ./s2xml.pl\n";

    printf("\n reads schema file event.xsd and generates the two hddm files\n");
    printf(" reaction.xml and hits.xml that are used as include files in event.xml,\n"); 
    printf(" which in turn can be used by hddm-c to generate the files hddm_s.c and\n");
    printf(" and hddm_s.h\n");

}

if (($ARGV[0] eq "--help")||
    ($ARGV[0] eq "--h")||
    ($ARGV[0] eq "-help")||
    ($ARGV[0] eq "-h")){

    Usage();
    exit;
}

my $sf = "event.xsd";
my $of1 = "event.out";

my $cmd = "/export/home/zihlmann/HallD/src/libraries/HDDM/schema-hddm $sf > $of1";

system($cmd);

if(! (-e $of1)){
    printf "Error generating temporary file $of1 from schema file $sf\n";
    exit;
}

my $INF = new FileHandle("<$of1") or
    die "Error open temporary file $of1\n";

my $headerline = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
my $ofr = "reaction.xml";
my $ofh = "hits.xml";

my $OUTF = new FileHandle();
my $FH = 0;
while(<$INF>){

    if (/\<reaction/){
	$OUTF = new FileHandle(">$ofr") or
	    die "Error open file $ofr\n";
	print $OUTF $headerline;
	$FH=1;
    }
    if (/\<hitView/){
	$OUTF = new FileHandle(">$ofh") or
	    die "Error open file $ofh\n";	
	print $OUTF $headerline;
	$FH=1;
    }
    if ($FH){
	print $OUTF $_;
    }
    if (/\<\/reaction/){
	undef $OUTF;
	$FH=0;
    }
     if (/\<\/hitView/){
	undef $OUTF;
	$FH=0;
    }
   
}
undef $INF;
system("rm -f $of1");
