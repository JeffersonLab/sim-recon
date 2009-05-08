#!/usr/bin/env perl
#
#  schema-hddm: create a hddm template based on a xml schema
#
#  Both templates and schemas are themselves well-formed xml documents,
#  so the transformation can be done automatically using XSLT tools.
#  This implementation uses Xalan from xml.apache.org to carry out the
#  transformation.  Some result is passed through the XercesDOMParser
#  just to check for well-formed xml and improve the readability of the
#  of the template.
#
#  Richard T. Jones
#  September 15, 2003
#  ------------------

use XML::Xerces;

sub Usage()
{
   print <<eoi;

Usage: schema-hddm  myfile.xsd

reads the schema contained in myfile.xsd and writes to stdout a hddm
template, which is just a xml file that shows the data structure.

eoi
}

if (@ARGV != 1) {
   Usage();
   exit 1;
}
$infile = $ARGV[0];
if (open(XSD,"<$infile") == 0) {
   die "Unable to open input file $infile\n";
}
close(XSD);
$base = `basename $infile`;
chop $base;
$tmpxml = ".tmp-$base.xml";

# generate the basic schema using XSL transform:
#
#  two translators: xalan-j (java) or xalan-c (c++)
#  chose one of the two following, comment the other
#$cmd = "$ENV{XALANCROOT}/bin/Xalan -o $tmpxml $infile schema-hddm.xsl";
$cmd = "$ENV{JAVAROOT}/bin/java org.apache.xalan.xslt.Process" .
       " -IN $infile -OUT $tmpxml -XSL schema-hddm.xsl";

if (system($cmd)) {
   print "command was:$cmd\n";
   die "hddm-schema: errors returned by Xalan, quitting\n";
}

# check the result for correctness and pretty-print it

$parser = XML::Xerces::XercesDOMParser->new();
$parser->setValidationScheme (0);
$parser->setDoNamespaces (1);
$parser->setCreateEntityReferenceNodes(1);
$parser->setDoSchema (1);

$ERROR_HANDLER = XML::Xerces::PerlErrorHandler->new();
$parser->setErrorHandler($ERROR_HANDLER);
eval {$parser->parse ($tmpxml)};
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
exit 0;
