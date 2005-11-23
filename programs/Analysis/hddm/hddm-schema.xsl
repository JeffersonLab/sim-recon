<?xml version="1.0" encoding="UTF-8"?> 
<xsl:transform xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
               xmlns:xs="http://www.w3.org/2001/XMLSchema"
               xmlns:hddm="http://www.gluex.org/hddm"
               version="1.0">
<xsl:output method="xml" version="1.0" encoding="iso-8859-1" indent="no"/>

<!-- Top-level processing: make sure document tag is HDDM -->

<xsl:template match="/">
  <xsl:for-each select="*">
    <xsl:if test="name() != 'HDDM'">
      <xsl:message terminate="yes">
  Format error: expected root tag HDDM, found <xsl:value-of select="name()"/>
      </xsl:message>
    </xsl:if>
  </xsl:for-each>

  <!-- create schema header -->

  <xs:schema targetNamespace="http://www.gluex.org/hddm"
             elementFormDefault="qualified">

  <!-- declare hddm data types -->

    <xs:simpleType name="int">
      <xs:restriction base="xs:int"/>
    </xs:simpleType>

    <xs:simpleType name="long">
      <xs:restriction base="xs:long"/>
    </xs:simpleType>

    <xs:simpleType name="float">
      <xs:restriction base="xs:float"/>
    </xs:simpleType>

    <xs:simpleType name="double">
      <xs:restriction base="xs:double"/>
    </xs:simpleType>

    <xs:simpleType name="boolean">
      <xs:restriction base="xs:boolean"/>
    </xs:simpleType>

    <xs:simpleType name="string">
      <xs:restriction base="xs:string"/>
    </xs:simpleType>

    <xs:simpleType name="anyURI">
      <xs:restriction base="xs:anyURI"/>
    </xs:simpleType>

    <xs:simpleType name="Particle_t">
      <xs:restriction base="xs:string">
        <xs:enumeration value="unknown"/>
        <xs:enumeration value="gamma"/>
        <xs:enumeration value="positron"/>
        <xs:enumeration value="electron"/>
        <xs:enumeration value="neutrino"/>
        <xs:enumeration value="mu+"/>
        <xs:enumeration value="mu-"/>
        <xs:enumeration value="pi0"/>
        <xs:enumeration value="pi+"/>
        <xs:enumeration value="pi-"/>
        <xs:enumeration value="kL"/>
        <xs:enumeration value="k+"/>
        <xs:enumeration value="k-"/>
        <xs:enumeration value="neutron"/>
        <xs:enumeration value="proton"/>
        <xs:enumeration value="antiProton"/>
        <xs:enumeration value="kS"/>
        <xs:enumeration value="eta"/>
        <xs:enumeration value="Lambda"/>
        <xs:enumeration value="Sigma+"/>
        <xs:enumeration value="Sigma0"/>
        <xs:enumeration value="Sigma-"/>
        <xs:enumeration value="Xi0"/>
        <xs:enumeration value="Xi-"/>
        <xs:enumeration value="Omega-"/>
        <xs:enumeration value="antiNeutron"/>
        <xs:enumeration value="antiLambda"/>
        <xs:enumeration value="antiSigma-"/>
        <xs:enumeration value="antiSigma0"/>
        <xs:enumeration value="antiSigma+"/>
        <xs:enumeration value="antiXi0"/>
        <xs:enumeration value="antiXi+"/>
        <xs:enumeration value="antiOmega+"/>
        <xs:enumeration value="rho0"/>
        <xs:enumeration value="rho+"/>
        <xs:enumeration value="rho-"/>
        <xs:enumeration value="omega"/>
        <xs:enumeration value="etaPrime"/>
        <xs:enumeration value="phi"/>
        <xs:enumeration value="a0(980)"/>
        <xs:enumeration value="f0(980)"/>
      </xs:restriction>
    </xs:simpleType>

  <!-- declare tags in order from primitive to complex -->
    
    <xsl:apply-templates select="hddm:HDDM"/> 
  </xs:schema>
</xsl:template>

<xsl:template match="//*">
  <xsl:if test="count(preceding::node()[name()=name(current())]) = 0">
    <xsl:apply-templates/>
    <xsl:element name="xs:element">
      <xsl:attribute name="name">
        <xsl:value-of select="name()"/>
      </xsl:attribute>
      <xsl:element name="xs:complexType">
        <xsl:if test="count(./*) > 0">
          <xs:sequence>
          <xsl:for-each select="./*">
            <xsl:element name="xs:element">
              <xsl:attribute name="ref">
                <xsl:text>hddm:</xsl:text>
                <xsl:value-of select="name()"/>
              </xsl:attribute>
              <xsl:if test="@minOccurs">
                <xsl:if test="@minOccurs &lt; 0">
                  <xsl:message terminate="yes">
  Format error: negative value not allowed for minOccurs attribute
                  </xsl:message>
                </xsl:if>
                <xsl:attribute name="minOccurs">
                  <xsl:value-of select="@minOccurs"/>
                </xsl:attribute>
              </xsl:if>
              <xsl:if test="@maxOccurs">
                <xsl:if test="@minOccurs">
                  <xsl:if test="@maxOccurs &lt; minOccurs">
                    <xsl:message terminate="yes">
  Format error: value for maxOccurs must not be less than minOccurs
                    </xsl:message>
                  </xsl:if>
                  <xsl:if test="@maxOccurs &lt; 1">
                    <xsl:message terminate="yes">
  Format error: maxOccurs must not be less than one
                    </xsl:message>
                  </xsl:if>
                </xsl:if>
                <xsl:attribute name="maxOccurs">
                  <xsl:value-of select="@maxOccurs"/>
                </xsl:attribute>
              </xsl:if>
            </xsl:element>
          </xsl:for-each>
          </xs:sequence>
        </xsl:if>

        <xsl:for-each select="@*[name() != 'minOccurs' and name() != 'maxOccurs']">
          <xsl:element name="xs:attribute">
            <xsl:attribute name="name">
              <xsl:value-of select="name()"/>
            </xsl:attribute>
            <xsl:attribute name="use">required</xsl:attribute>
            <xsl:choose>
              <xsl:when test=". = 'int'">
                <xsl:attribute name="type">hddm:int</xsl:attribute>
              </xsl:when>
              <xsl:when test=". = 'long'">
                <xsl:attribute name="type">hddm:long</xsl:attribute>
              </xsl:when>
              <xsl:when test=". = 'float'">
                <xsl:attribute name="type">hddm:float</xsl:attribute>
              </xsl:when>
              <xsl:when test=". = 'double'">
                <xsl:attribute name="type">hddm:double</xsl:attribute>
              </xsl:when>
              <xsl:when test=". = 'boolean'">
                <xsl:attribute name="type">hddm:boolean</xsl:attribute>
              </xsl:when>
              <xsl:when test=". = 'string'">
                <xsl:attribute name="type">hddm:string</xsl:attribute>
              </xsl:when>
              <xsl:when test=". = 'anyURI'">
                <xsl:attribute name="type">hddm:anyURI</xsl:attribute>
              </xsl:when>
              <xsl:when test=". = 'Particle_t'">
                <xsl:attribute name="type">hddm:Particle_t</xsl:attribute>
              </xsl:when>
              <xsl:otherwise>
                <xsl:attribute name="fixed">
                  <xsl:value-of select="."/>
                </xsl:attribute>
                <xsl:attribute name="type">xs:string</xsl:attribute>
              </xsl:otherwise>
            </xsl:choose>
          </xsl:element>
        </xsl:for-each>
      </xsl:element>
    </xsl:element>
  </xsl:if>
</xsl:template>

<!-- ignore all plain text -->

<xsl:template match="text()">
</xsl:template>

</xsl:transform>
