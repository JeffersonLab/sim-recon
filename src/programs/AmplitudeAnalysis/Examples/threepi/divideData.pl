#!/usr/bin/perl

$lowMass = 0.7;
$highMass = 2.0;
$nBins = 65;

# put a limit on the number of data events to process
# gen MC and acc MC smaples are not limited
$maxEvts = 100000;

$bin = "/Users/mashephe/iu_cvs/GlueXAmpExe";
$workingDir = "/Users/mashephe/amptools/examples/threepi_example";

$dataFile = "$workingDir/threepi_res_0pol.root";
$accMCFile = "$workingDir/threepi_flat.root";
$genMCFile = "$workingDir/threepi_flat.root";
$cfgTempl = "$workingDir/threepi_unpol_TEMPLATE.cfg";

$fitDir = "$workingDir/threepi_0pol_fit/";

if( ! -d $fitDir ){

  mkdir $fitDir;
}

chdir $fitDir;

for( $i = 0; $i < $nBins; ++$i ){

  mkdir "bin_$i" unless -d "bin_$i";
}

@dataParts = split /\//, $dataFile;
$dataTag = pop @dataParts;
$dataTag =~ s/\.root//;
system( "$bin/split_mass $dataFile $dataTag $lowMass $highMass $nBins $maxEvts" );

@accMCParts = split /\//, $accMCFile;
$accMCTag = pop @accMCParts;
$accMCTag =~ s/\.root//;
system( "$bin/split_mass $accMCFile $accMCTag $lowMass $highMass $nBins" );

@genMCParts = split /\//, $genMCFile;
$genMCTag = pop @genMCParts;
$genMCTag =~ s/\.root//;
system( "$bin/split_mass $genMCFile $genMCTag $lowMass $highMass $nBins" );

for( $i = 0; $i < $nBins; ++$i ){

  system( "mv *\_$i.root bin_$i" );

  chdir "bin_$i";

  open( CFGOUT, ">bin_$i.cfg" );
  open( CFGIN, $cfgTempl );

  while( <CFGIN> ){

    s/DATAFILE/$dataTag\_$i.root/;
    s/ACCMCFILE/$accMCTag\_$i.root/;
    s/GENMCFILE/$genMCTag\_$i.root/;
    s/NIFILE/bin_$i.ni/;
    s/FITNAME/bin_$i/;

    print CFGOUT $_;
  }

  close CFGOUT;
  close CFGIN;

  system( "touch param_init.cfg" );

  chdir $fitDir;
}

