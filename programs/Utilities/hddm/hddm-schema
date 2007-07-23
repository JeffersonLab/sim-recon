#!/usr/bin/env perl
#
#  hddm-schema: create a xml schema based on a hddm template
#
#  Both templates and schemas are themselves well-formed xml documents,
#  so the transformation can be done automatically using XSLT tools.
#  This implementation uses Xalan from xml.apache.org to carry out the
#  transformation.  The output from the transform is then validated by
#  the XercesDOMParser before printing it on standard output.
#
#  Richard T. Jones
#  September 15, 2003
#  ------------------

use XML::Xerces;

sub Usage()
{
   print <<eoi;

Usage: hddm-schema myfile.hddm

extracts metadata from the header of myfile.hddm and writes to stdout
a schema describing the contents.

eoi
}

if (@ARGV != 1) {
   Usage();
   exit 1;
}
$infile = $ARGV[0];
if (open(HDDM,"<$infile") == 0) {
   die "Unable to open input file $infile\n";
}
$base = `basename $infile`;
chop $base;
$tmpxml = ".tmp-$base.xml";
$tmpxsd = ".tmp-$base.xsd";
if (open(TMP,">$tmpxml") == 0) {
   die "Unable to open temp file $tmpxml\n";
}
while ($line = <HDDM>) {
   if ($line =~ s/<\/HDDM>.*$/<\/HDDM>/) {
      print TMP $line;
      last;
   }
   print TMP $line;
}
close(HDDM);
close(TMP);

# generate the basic schema using XSL transform:
#
#  two translators: xalan-j (java) or xalan-c (c++)
#  chose one of the two following, comment the other
#$cmd = "$ENV{XALANCROOT}/bin/Xalan -o $tmpxsd $tmpxml hddm-schema.xsl";
$cmd = "$ENV{JAVAROOT}/bin/java org.apache.xalan.xslt.Process" .
       " -IN $tmpxml -OUT $tmpxsd -XSL hddm-schema.xsl -L";

if (system($cmd)) {
   if (! -x "$ENV{JAVAROOT}/bin/java") {
      die "Please check that $JAVAROOT/bin/java exists and try again.\n";
   }
   else {
      print "command was:$cmd\n";
      die "hddm-schema: errors returned by Xalan, quitting\n";
   }
}

# check the result for correctness and pretty-print it

$parser = XML::Xerces::XercesDOMParser->new();
$parser->setValidationScheme (0);
$parser->setDoNamespaces (1);
$parser->setCreateEntityReferenceNodes(1);
$parser->setDoSchema (1);

$ERROR_HANDLER = XML::Xerces::PerlErrorHandler->new();
$parser->setErrorHandler($ERROR_HANDLER);
eval {$parser->parse ($tmpxsd)};
XML::Xerces::error($@) if ($@);

$doc = $parser->getDocument();

$impl = XML::Xerces::DOMImplementationRegistry::getDOMImplementation('LS');
$writer = $impl->createDOMWriter();
if ($writer->canSetFeature('format-pretty-print',1)) {
  $writer->setFeature('format-pretty-print',1);
}
$target = XML::Xerces::StdOutFormatTarget->new();
$writer->writeNode($target,$doc);

unlink $tmpxml;
unlink $tmpxsd;
exit 0;
